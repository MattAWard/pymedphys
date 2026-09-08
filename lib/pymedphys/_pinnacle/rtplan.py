# Copyright (C) 2019 South Western Sydney Local Health District,
# University of New South Wales

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This work is derived from:
# https://github.com/AndrewWAlexander/Pinnacle-tar-DICOM
# which is released under the following license:

# Copyright (c) [2017] [Colleen Henschel, Andrew Alexander]

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import os
import re
import time

from pymedphys._imports import pydicom
from pymedphys._pinnacle.pinnacle_exceptions import (
    IsocenterNotFoundError,
    MachineDataNotFoundError,
    MissingCTImageError,
    MissingTrialBeamsError,
)

from .constants import (
    GImplementationClassUID,
    GTransferSyntaxUID,
    RTPLANModality,
    RTPlanSOPClassUID,
    RTStructSOPClassUID,
)
from .pinnacle_metadata import (
    append_pinnacle_metadata_for_plan,
    apply_approval_status,
    apply_equipment_stamps,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _new_dataset():
    """Shorthand for creating a new empty DICOM Dataset."""
    return pydicom.dataset.Dataset()


def _new_sequence():
    """Shorthand for creating a new empty DICOM Sequence."""
    return pydicom.sequence.Sequence()


def _sanitize_for_filename(name):
    """Make a trial name safe for use in a DICOM file name."""
    return re.sub(r"[^\w\-.]", "_", str(name)) if name else "trial"


# DICOM LO (Long String) maximum length — used for SeriesDescription.
_DICOM_LO_MAX = 64

# DICOM SH (Short String) maximum length — used for RTPlanLabel.
_DICOM_SH_MAX = 16


def _truncate_sh(value, logger=None, tag=""):
    """Truncate a value to the DICOM SH (Short String) limit of 16 chars.

    Pinnacle plan names regularly exceed 16 characters (e.g.
    "CopyOf_1_LtBreast" -> RTPlanLabel "CopyOf_1_LtBreast.0" = 19), which
    pydicom warns about and which some PACS reject outright.  The full,
    untruncated name is always still available in the LO-VR fields
    (RTPlanName / StructureSetName / SeriesDescription), so nothing is
    lost — only the short label is clipped.
    """
    text = str(value or "")
    if len(text) <= _DICOM_SH_MAX:
        return text
    truncated = text[:_DICOM_SH_MAX]
    if logger is not None:
        logger.warning(
            "%s value %r is %d chars, exceeding the DICOM SH limit of %d; "
            "truncated to %r (the full name is retained in the "
            "corresponding LO-VR attribute).",
            tag or "SH", text, len(text), _DICOM_SH_MAX, truncated,
        )
    return truncated


def _format_ds(value):
    """Format a float for a DICOM DS (Decimal String, max 16 chars)."""
    text = f"{float(value):.6f}".rstrip("0").rstrip(".")
    if text in ("", "-"):
        text = "0"
    return text[:16]


def _append_trial_series_description(ds, trial_name, prefix):
    """Append 'prefix: trial' to SeriesDescription, respecting the LO VR.

    Exported objects carry their trial name so that a reviewer
    can distinguish trials within a plan at the PACS/TPS end.
    """
    label = f"{prefix}: {trial_name}"
    base = (getattr(ds, "SeriesDescription", "") or "").strip()
    combined = f"{base} - {label}" if base else label
    if len(combined) > _DICOM_LO_MAX:
        combined = combined[: _DICOM_LO_MAX - 3].rstrip() + "..."
    ds.SeriesDescription = combined


def _iter_machine_dicts(node, _depth=0):
    """Yield every dict in *node* that looks like a Pinnacle machine entry.

    ``pinn_to_dict`` nests top-level Pinnacle objects under their type name
    — ``plan.Trial`` parses to ``{"Trial": {...}}``, and
    ``plan.Pinnacle.Machines`` likewise wraps each machine (single machine,
    a list of machines, or a keyed container, depending on the file).  The
    previous implementation only looked at the top level, so
    ``machine.get("Name")`` returned ``None`` and every lookup failed —
    which in turn blocked SAD, DosePerMuAtCalibration and (since the MLC
    fallback was removed) the entire RTPLAN export.

    Rather than hard-coding one nesting shape, this walks the parsed
    structure and yields any dict carrying a "Name" key alongside at least
    one machine-ish attribute.  Depth is bounded so a pathological file
    cannot cause runaway recursion.
    """
    if _depth > 6:
        return

    if isinstance(node, dict):
        keys = set(node)
        if "Name" in keys and keys & {
            "PhotonEnergyList",
            "ElectronEnergyList",
            "MultiLeaf",
            "MultiLeafLayout",
            "VersionTimestamp",
            "SourceToAxisDistance",
            "SourceAxisDistance",
            "SAD",
        }:
            yield node
        for value in node.values():
            yield from _iter_machine_dicts(value, _depth + 1)

    elif isinstance(node, list):
        for item in node:
            yield from _iter_machine_dicts(item, _depth + 1)


def _select_machine(machine_info, machinename, machineversion, logger=None):
    """Return the machine dict matching *machinename*/*machineversion*.

    Searches the parsed ``plan.Pinnacle.Machines`` structure at any nesting
    depth (see :func:`_iter_machine_dicts`).  Returns ``None`` when no match
    is found — callers must fall back to safe behaviour (and, for the MLC
    boundary table, refuse the export rather than assume geometry).
    """
    if machine_info is None:
        return None

    candidates = list(_iter_machine_dicts(machine_info))
    if not candidates:
        if logger is not None:
            logger.warning(
                "No machine entries could be located in the parsed "
                "plan.Pinnacle.Machines structure (top-level keys: %s). The "
                "file may use an unexpected layout.",
                sorted(machine_info)[:10] if isinstance(machine_info, dict)
                else type(machine_info).__name__,
            )
        return None

    # 1. Exact match on both name and version timestamp.
    for machine in candidates:
        if machine.get("Name") == machinename and (
            not machineversion
            or machine.get("VersionTimestamp") == machineversion
        ):
            return machine

    # 2. Name-only match (version timestamps occasionally differ between the
    #    Trial reference and the Machines file).
    for machine in candidates:
        if machine.get("Name") == machinename:
            if logger is not None:
                logger.debug(
                    "Machine '%s' matched by name only (version '%s' not "
                    "matched exactly; file has '%s').",
                    machinename, machineversion,
                    machine.get("VersionTimestamp"),
                )
            return machine

    # No match — report what IS in the file so the mismatch is diagnosable
    # without hand-parsing the archive.
    if logger is not None:
        available = [
            f"{m.get('Name')!r} (version {m.get('VersionTimestamp')!r})"
            for m in candidates[:8]
        ]
        logger.warning(
            "Machine '%s' (version '%s') not found among the %d machine "
            "entry/entries in plan.Pinnacle.Machines. Available: %s",
            machinename, machineversion, len(candidates),
            "; ".join(available) or "(none)",
        )
    return None


def _get_machine_sad_mm(machine, logger, beam_name):
    """Read the source-axis distance (mm) from Pinnacle machine data.

    Pinnacle stores geometry in cm; several key spellings are probed.
    Returns ``None`` when the value is absent or fails a
    sanity check, in which case the caller falls back to 1000 mm with a
    logged assumption.
    """
    if not isinstance(machine, dict):
        return None

    for key in ("SourceToAxisDistance", "SourceAxisDistance", "SAD"):
        raw = machine.get(key)
        if raw in (None, ""):
            continue
        try:
            sad_mm = float(raw) * 10  # Pinnacle cm → DICOM mm
        except (TypeError, ValueError):
            continue
        # Sanity window: clinical linac SADs sit comfortably within
        # 500–3000 mm.  Anything outside suggests a unit/parse problem.
        if 500 <= sad_mm <= 3000:
            return sad_mm
        logger.warning(
            "Beam '%s': machine SAD value %s (key '%s') is outside the "
            "plausible 500–3000 mm range after cm→mm conversion; ignoring "
            "and falling back to 1000 mm.",
            beam_name,
            sad_mm,
            key,
        )
        return None

    return None


def _leaf_boundaries_from_machine(machine, expected_pairs, logger):
    """Build LeafPositionBoundaries (mm strings) from Pinnacle machine data.

    Derives the boundary table from the machine's MultiLeaf layout
    (leaf-pair centre positions and widths, stored in cm) instead
    of assuming the Varian-Millennium pattern.

    Returns a list of ``expected_pairs + 1`` boundary strings, or ``None``
    when the machine data is absent/inconsistent — in which case the
    trial's RTPLAN export fails (MachineDataNotFoundError) rather than
    exporting with an assumed boundary table.
    """
    if not isinstance(machine, dict):
        return None

    multileaf = machine.get("MultiLeaf") or machine.get("MultiLeafLayout")
    if not isinstance(multileaf, dict):
        return None

    pair_container = (
        multileaf.get("LeafPairList")
        or multileaf.get("LeafPairArray")
        or multileaf.get("LeafPairs")
    )
    if pair_container is None:
        return None

    # Normalise the parser output into a flat list of leaf-pair dicts.
    if isinstance(pair_container, dict):
        pairs = [v for v in pair_container.values() if isinstance(v, dict)]
        # Some parses nest each pair under a repeated "LeafPair" key that
        # collapses to a single dict/list — handle a list value too.
        if not pairs:
            inner = pair_container.get("LeafPair")
            if isinstance(inner, list):
                pairs = [p for p in inner if isinstance(p, dict)]
            elif isinstance(inner, dict):
                pairs = [inner]
    elif isinstance(pair_container, list):
        pairs = [p for p in pair_container if isinstance(p, dict)]
    else:
        return None

    geometry = []
    for pair in pairs:
        center = pair.get("YCenterPosition", pair.get("CenterPosition"))
        width = pair.get("Width", pair.get("LeafWidth"))
        if center is None or width is None:
            return None
        try:
            geometry.append((float(center) * 10, float(width) * 10))  # cm → mm
        except (TypeError, ValueError):
            return None

    if len(geometry) != expected_pairs:
        logger.warning(
            "Machine MultiLeaf data describes %d leaf pairs but the plan "
            "data implies %d; the machine boundary table cannot be used and "
            "this trial's RTPLAN export will fail.",
            len(geometry),
            expected_pairs,
        )
        return None

    geometry.sort(key=lambda cw: cw[0])

    boundaries = [geometry[0][0] - geometry[0][1] / 2]
    for center, width in geometry:
        boundaries.append(center + width / 2)

    # The table must be strictly increasing to be a valid boundary set.
    for a, b in zip(boundaries, boundaries[1:]):
        if b <= a:
            logger.warning(
                "Machine MultiLeaf boundaries are not strictly increasing "
                "(%s then %s); the machine boundary table cannot be used and "
                "this trial's RTPLAN export will fail.",
                a,
                b,
            )
            return None

    return [_format_ds(b) for b in boundaries]


def _resolve_beam_isocenter(plan, beam):
    """Resolve the isocenter for a specific beam.

    Priority order:

    1. The beam's ``IsocenterName`` from ``plan.Trial``, looked up in
       ``plan.Points`` (case-insensitive).  A named-but-missing point is
       an error: assuming another point would be dangerous.
    2. The plan-level heuristic (``find_iso_center``: PoiInterpretedType,
       iso-like names, CT centre) for archives whose trial data carries
       no isocenter name.

    Raises
    ------
    IsocenterNotFoundError
        When no isocenter can be resolved.
    """
    iso_name = str(beam.get("IsocenterName", "") or "").strip()
    if iso_name:
        for point in plan.points:
            if str(point.get("Name", "")).strip().lower() == iso_name.lower():
                iso = plan.convert_point(point)
                plan.logger.debug(
                    "Beam '%s': isocenter '%s' resolved from plan.Points: %s",
                    beam.get("Name"),
                    iso_name,
                    iso,
                )
                return iso
        raise IsocenterNotFoundError(
            f"Beam '{beam.get('Name')}' references isocenter point "
            f"'{iso_name}' which was not found in plan.Points."
        )

    # No IsocenterName in this archive — fall back to the plan-level
    # heuristic.  plan.iso_center lazily runs find_iso_center, which no
    # longer defaults to an arbitrary first point (returns [] instead).
    iso = plan.iso_center
    if iso is not None and len(iso) >= 3:
        plan.logger.debug(
            "Beam '%s': no IsocenterName in trial data; using plan-level "
            "isocenter heuristic: %s",
            beam.get("Name"),
            iso,
        )
        return iso

    raise IsocenterNotFoundError(
        f"No isocenter could be determined for beam '{beam.get('Name')}': "
        f"the trial data carries no IsocenterName and no isocenter-like "
        f"point exists in plan.Points."
    )


def _gantry_direction_between(prev_angle, next_angle):
    """Return the DICOM rotation direction from *prev_angle* to *next_angle*.

    Angles are in degrees; the shorter arc decides the
    direction, with the delta normalised into (-180, 180].
    """
    try:
        delta = (float(next_angle) - float(prev_angle) + 180.0) % 360.0 - 180.0
    except (TypeError, ValueError):
        return "NONE"
    if delta > 1e-6:
        return "CW"
    if delta < -1e-6:
        return "CC"
    return "NONE"


def _parse_mlc_leaf_positions(control_point):
    """Parse MLC leaf positions from a Pinnacle control point dict.

    Returns (leafpositions, p_count) where leafpositions is the interleaved
    list ready for DICOM and p_count is the total number of raw leaf values.
    """
    points_str = control_point["MLCLeafPositions"]["RawData"]["Points[]"]
    raw_points = points_str.split(",")
    p_count = len(raw_points)

    bank_a = []  # left bank
    bank_b = []  # right bank
    for i, p in enumerate(raw_points):
        leafpoint = float(p.strip())
        if i % 2 == 0:
            bank_a.append(-leafpoint * 10)
        else:
            bank_b.append(leafpoint * 10)

    # Reverse both banks and concatenate
    bank_a = list(reversed(bank_a))
    bank_b = list(reversed(bank_b))
    return bank_a + bank_b, p_count


def _parse_wedge_info(cp_data, plan_logger):
    """Parse wedge information from a Pinnacle control point.

    Returns a dict with wedge details, or None if no wedge is present.
    Keys: type, angle, name, orientation, count.
    """
    wedge_context = cp_data["WedgeContext"]
    wedge_name_raw = wedge_context["WedgeName"]

    if wedge_name_raw in ("No Wedge", ""):
        plan_logger.debug("No wedge present")
        return None

    info = {"count": 1, "angle": wedge_context["Angle"]}
    orientation_raw = wedge_context["Orientation"]

    if "edw" in wedge_name_raw.lower():
        # Enhanced Dynamic Wedge
        plan_logger.debug("EDW wedge present")
        info["type"] = "DYNAMIC"
        if orientation_raw == "WedgeBottomToTop":
            info["name"] = f"{wedge_name_raw.upper()}{info['angle']}IN"
            info["orientation"] = "0"  # TODO confirm orientation mapping
        elif orientation_raw == "WedgeTopToBottom":
            info["name"] = f"{wedge_name_raw.upper()}{info['angle']}OUT"
            info["orientation"] = "180"
        plan_logger.debug("EDW wedge name = %s", info.get("name"))

    elif "UP" in wedge_name_raw:
        # Standard (Universal/Physical) wedge
        plan_logger.debug("Standard wedge present")
        info["type"] = "STANDARD"
        angle_int = int(info["angle"])

        # Map wedge angle to the machine-specific number suffix
        angle_to_suffix = {15: "30", 30: "30", 45: "20", 60: "15"}
        number_suffix = angle_to_suffix.get(angle_int, "")

        orientation_to_label = {
            "WedgeRightToLeft": ("R", "90"),
            "WedgeLeftToRight": ("L", "270"),
            "WedgeTopToBottom": ("OUT", "180"),
            "WedgeBottomToTop": ("IN", "0"),
        }
        label, dicom_orientation = orientation_to_label.get(orientation_raw, ("", "0"))
        info["name"] = f"W{angle_int}{label}{number_suffix}"
        info["orientation"] = dicom_orientation  # TODO: confirm orientation values
        plan_logger.debug("Standard wedge name = %s", info["name"])
    else:
        # Unknown wedge type — treat as no wedge
        plan_logger.warning("Unrecognised wedge name: %s", wedge_name_raw)
        return None

    return info


def _populate_beam_limiting_device_seq(
    beam_ds, p_count, logger=None, machine_boundaries=None
):
    """Populate the BeamLimitingDeviceSequence for a beam.

    Creates entries for ASYMX and ASYMY, plus MLCX when the beam has MLC
    leaf data.  The MLCX boundary table MUST come from the Pinnacle
    machine data (*machine_boundaries*): exporting with an
    assumed table (the old hardcoded Varian-Millennium fallback) is
    dangerous — a wrong table silently shifts every leaf pair — so a
    beam that uses an MLC without derivable machine geometry fails the
    trial's RTPLAN export instead.

    Raises
    ------
    MachineDataNotFoundError
        When p_count > 0 but no consistent boundary table was derived
        from the machine data.
    """
    beam_ds.BeamLimitingDeviceSequence = _new_sequence()

    asymx = _new_dataset()
    asymx.RTBeamLimitingDeviceType = "ASYMX"
    asymx.NumberOfLeafJawPairs = "1"
    beam_ds.BeamLimitingDeviceSequence.append(asymx)

    asymy = _new_dataset()
    asymy.RTBeamLimitingDeviceType = "ASYMY"
    asymy.NumberOfLeafJawPairs = "1"
    beam_ds.BeamLimitingDeviceSequence.append(asymy)

    # NumberOfLeafJawPairs has VR=IS (integer); use integer division so we
    # never emit a fractional value such as "60.0".
    num_pairs = p_count // 2

    if num_pairs == 0:
        # Jaw-only beam: no MLCX device entry (valid DICOM without one).
        if logger is not None:
            logger.debug(
                "Beam has no MLC leaf data; MLCX omitted from "
                "BeamLimitingDeviceSequence."
            )
        return

    if machine_boundaries is None or len(machine_boundaries) != num_pairs + 1:
        raise MachineDataNotFoundError(
            f"MLC LeafPositionBoundaries for {num_pairs} leaf pairs could "
            f"not be derived from the Pinnacle machine data "
            f"(plan.Pinnacle.Machines). Exporting with an assumed boundary "
            f"table is unsafe, so this trial's RTPLAN export is refused. "
            f"Verify the machine file is present in the archive and its "
            f"MultiLeaf layout matches the plan's leaf count."
        )

    mlcx = _new_dataset()
    mlcx.RTBeamLimitingDeviceType = "MLCX"
    mlcx.NumberOfLeafJawPairs = num_pairs
    # IS-725: boundaries derived from the Pinnacle machine MultiLeaf data.
    mlcx.LeafPositionBoundaries = machine_boundaries
    beam_ds.BeamLimitingDeviceSequence.append(mlcx)


def _create_bld_position_entries(x1, x2, y1, y2, leafpositions):
    """Create the BeamLimitingDevicePositionSequence items for a control point.

    Returns a Sequence containing ASYMX and ASYMY entries, plus an MLCX
    entry when the beam has MLC leaf data (jaw-only beams legitimately
    have no MLCX device).
    """
    bld_seq = _new_sequence()

    asymx = _new_dataset()
    asymx.RTBeamLimitingDeviceType = "ASYMX"
    asymx.LeafJawPositions = [x1, x2]
    bld_seq.append(asymx)

    asymy = _new_dataset()
    asymy.RTBeamLimitingDeviceType = "ASYMY"
    asymy.LeafJawPositions = [y1, y2]
    bld_seq.append(asymy)

    if leafpositions:
        mlcx = _new_dataset()
        mlcx.RTBeamLimitingDeviceType = "MLCX"
        mlcx.LeafJawPositions = leafpositions
        bld_seq.append(mlcx)

    return bld_seq


def _create_wedge_position_seq():
    """Create a WedgePositionSequence for a control point with wedge IN."""
    seq = _new_sequence()
    wp = _new_dataset()
    wp.WedgePosition = "IN"
    wp.ReferencedWedgeNumber = "1"
    seq.append(wp)
    return seq


def _create_wedge_sequence(wedge_info):
    """Create the beam-level WedgeSequence from parsed wedge info."""
    seq = _new_sequence()
    wedge = _new_dataset()
    wedge.WedgeNumber = 1
    wedge.WedgeType = wedge_info["type"]
    wedge.WedgeAngle = wedge_info["angle"]
    wedge.WedgeID = wedge_info["name"]
    wedge.WedgeOrientation = wedge_info["orientation"]
    # WedgeFactor (Type 3, VR=DS) is omitted rather than written as "": an
    # empty string is not a valid Decimal String. Populate with the real
    # factor from Pinnacle data when it becomes available.
    seq.append(wedge)
    return seq


def _populate_first_control_point(
    cp,
    beam_ds,
    beam,
    plan,
    beam_energy,
    doserate,
    gantryrotdir,
    numwedges,
    cp_entry,
    iso_center,
):
    """Populate all attributes required by DICOM for the first control point.

    Per DICOM C.8.8.14.5: at the first control point, ALL applicable
    attributes must be present. This includes 1C and 2C attributes.

    All geometric values now come from *cp_entry* — the parsed data of
    the beam's first Pinnacle control point — rather than values
    frozen from whichever control point was parsed last.
    """
    # --- Required energy and dose rate (Type 3 but universally expected) ---
    cp.NominalBeamEnergy = beam_energy
    cp.DoseRateSet = doserate

    # --- Gantry (1C — required at first CP) ---
    cp.GantryAngle = cp_entry["gantry"]
    cp.GantryRotationDirection = gantryrotdir

    # --- Collimator (1C — required at first CP) ---
    cp.BeamLimitingDeviceAngle = cp_entry["collimator"]
    cp.BeamLimitingDeviceRotationDirection = "NONE"

    # --- Patient Support / Couch (1C — required at first CP) ---
    cp.PatientSupportAngle = cp_entry["couch"]
    cp.PatientSupportRotationDirection = "NONE"

    # --- Table Top Eccentric (1C — required at first CP) ---
    cp.TableTopEccentricAngle = "0"
    cp.TableTopEccentricRotationDirection = "NONE"

    # --- Table Top Position (2C — required at first CP, may be empty) ---
    cp.TableTopVerticalPosition = ""
    cp.TableTopLongitudinalPosition = ""
    cp.TableTopLateralPosition = ""

    # --- Isocenter (2C — required at first CP when isocentric) ---
    # Resolved per beam from the trial's IsocenterName.
    cp.IsocenterPosition = iso_center

    # --- Source to Surface Distance (Type 3) ---
    cp.SourceToSurfaceDistance = beam["SSD"] * 10

    # --- Wedge position (1C — required when wedges present) ---
    if numwedges > 0:
        cp.WedgePositionSequence = _create_wedge_position_seq()

    # --- Beam Limiting Device positions (1C) ---
    cp.BeamLimitingDevicePositionSequence = _create_bld_position_entries(
        cp_entry["x1"],
        cp_entry["x2"],
        cp_entry["y1"],
        cp_entry["y2"],
        cp_entry["leafpositions"],
    )

    # --- Beam-level counts (Type 1 — placed here for locality but belong to beam) ---
    beam_ds.NumberOfWedges = numwedges
    beam_ds.NumberOfCompensators = "0"
    beam_ds.NumberOfBoli = "0"
    beam_ds.NumberOfBlocks = "0"


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------


def convert_plan(plan, export_path):
    """Export RTPLAN files for every trial in the plan.

    For each trial: switch the plan's active trial, generate fresh per-trial
    UIDs, and call ``convert_plan_for_trial`` to write a single RTPLAN file.
    The matching RTSTRUCT UID for that trial is also generated here so the
    written RTPLAN's ``ReferencedStructureSetSequence`` lines up with the
    RTSTRUCT that ``convert_struct`` will emit for the same trial.
    """
    if not plan.primary_image:
        plan.logger.error("No primary image found for plan. Unable to generate RTPLAN.")
        raise MissingCTImageError("Plan has no primary image associated with it.")

    # TODO Test the RTPLAN export functionality and remove this warning
    plan.logger.warning(
        "RTPLAN export functionality is currently not validated and not stable. "
        "Use with caution."
    )

    for trial_info in plan.trials:
        plan.active_trial = trial_info["Name"]
        plan.logger.info("Exporting RTPLAN for trial: %s", trial_info["Name"])

        uids = plan.generate_uids_for_trial(trial_info)
        try:
            convert_plan_for_trial(
                plan,
                trial_info,
                plan_instance_uid=uids["plan"],
                struct_instance_uid=uids["struct"],
                series_instance_uid=uids["series_plan"],
                export_path=export_path,
            )
        except MissingTrialBeamsError as exc:
            plan.logger.warning(
                "Skipping RTPLAN for trial '%s': %s", trial_info["Name"], exc
            )
            continue
        except (IsocenterNotFoundError, MachineDataNotFoundError) as exc:
            # Refusing to guess geometry; surface
            # loudly so the operator knows this trial's RTPLAN was refused.
            plan.logger.error(
                "RTPLAN export failed for trial '%s': %s", trial_info["Name"], exc
            )
            continue


# ---------------------------------------------------------------------------
# Per-trial RTPLAN generation
# ---------------------------------------------------------------------------


def convert_plan_for_trial(
    plan,
    trial_info,
    plan_instance_uid,
    struct_instance_uid,
    series_instance_uid,
    export_path,
):
    """Write a single RTPLAN DICOM file for one specific trial."""

    patient_info = plan.pinnacle.patient_info
    plan_info = plan.plan_info
    image_info = plan.primary_image.image_info[0]
    machine_info = plan.machine_info
    patient_position = plan.patient_position

    # --- File meta ---
    file_meta = _new_dataset()
    file_meta.MediaStorageSOPClassUID = RTPlanSOPClassUID
    file_meta.TransferSyntaxUID = GTransferSyntaxUID
    file_meta.MediaStorageSOPInstanceUID = plan_instance_uid
    file_meta.ImplementationClassUID = GImplementationClassUID

    safe_trial = _sanitize_for_filename(trial_info.get("Name"))
    rp_filename = f"RP.{safe_trial}.{plan_instance_uid}.dcm"
    ds = pydicom.dataset.FileDataset(
        rp_filename, {}, file_meta=file_meta, preamble=b"\x00" * 128
    )

    # --- Study / Patient level ---
    ds.SpecificCharacterSet = "ISO_IR 100"
    ds.InstanceCreationDate = time.strftime("%Y%m%d")
    ds.InstanceCreationTime = time.strftime("%H%M%S")
    ds.SOPClassUID = RTPlanSOPClassUID
    ds.SOPInstanceUID = plan_instance_uid

    datetimesplit = plan_info["ObjectVersion"]["WriteTimeStamp"].split()
    if trial_info:
        datetimesplit = trial_info["ObjectVersion"]["WriteTimeStamp"].split()

    ds.StudyDate = datetimesplit[0].replace("-", "")
    ds.StudyTime = datetimesplit[1].replace(":", "")
    ds.AccessionNumber = ""
    ds.Modality = RTPLANModality
    ds.Manufacturer = ""  # Type 2; overwritten by apply_equipment_stamps
    ds.OperatorsName = ""
    ds.ManufacturerModelName = plan_info.get("ToolType", "")
    ds.SoftwareVersions = [plan_info["PinnacleVersionDescription"]]

    # Apply site-specific equipment identification stamps from config
    apply_equipment_stamps(
        ds,
        plan.pinnacle.equipment_cfg,
        pinnacle_model=plan_info.get("ToolType", ""),
        pinnacle_sw=plan_info.get("PinnacleVersionDescription", ""),
    )
    ds.PhysiciansOfRecord = patient_info["RadiationOncologist"]
    ds.PatientName = patient_info["FullName"]
    ds.PatientBirthDate = patient_info["DOB"]
    ds.PatientID = patient_info["MedicalRecordNumber"]
    ds.PatientSex = patient_info.get("Gender", "")[:1]
    ds.StudyInstanceUID = image_info["StudyInstanceUID"]
    ds.SeriesInstanceUID = series_instance_uid
    ds.StudyID = plan.primary_image.image["StudyID"]
    ds.FrameOfReferenceUID = image_info["FrameUID"]
    ds.PositionReferenceIndicator = ""

    # Carry the trial name so reviewers can distinguish
    # trials within a plan at the PACS/TPS end.
    _append_trial_series_description(ds, trial_info.get("Name", ""), "Plan")

    # --- Plan identification ---
    # RTPlanLabel is VR SH (max 16 chars); RTPlanName is LO (max 64) and
    # keeps the full untruncated Pinnacle plan name.
    ds.RTPlanLabel = _truncate_sh(
        f"{plan_info['PlanName']}.0", plan.logger, "RTPlanLabel")
    ds.RTPlanName = plan_info["PlanName"]
    ds.RTPlanDescription = append_pinnacle_metadata_for_plan(
        None,
        plan,
        trial_info,
        max_length=1024,
    )
    ds.RTPlanDate = ds.StudyDate
    ds.RTPlanTime = ds.StudyTime
    ds.PlanIntent = ""  # Type 3 — no curative/palliative source in Pinnacle data
    # PATIENT is correct because this RTPLAN always carries a
    # ReferencedStructureSetSequence pointing at an image-based RTSTRUCT
    # (the exporter refuses to run without a primary CT image).
    ds.RTPlanGeometry = "PATIENT"

    # --- Referenced Structure Set ---
    ds.ReferencedStructureSetSequence = _new_sequence()
    ref_struct = _new_dataset()
    ref_struct.ReferencedSOPClassUID = RTStructSOPClassUID
    ref_struct.ReferencedSOPInstanceUID = struct_instance_uid
    ds.ReferencedStructureSetSequence.append(ref_struct)

    # Derived from Pinnacle PlanLockStatus (locked → APPROVED
    # with reviewer/timestamp audit fields, unlocked → UNAPPROVED).
    apply_approval_status(ds, plan)

    # --- Fraction Group ---
    ds.FractionGroupSequence = _new_sequence()
    fraction_group = _new_dataset()
    fraction_group.ReferencedBeamSequence = _new_sequence()
    ds.FractionGroupSequence.append(fraction_group)

    # --- Sequences that are populated per-beam ---
    ds.BeamSequence = _new_sequence()
    ds.PatientSetupSequence = _new_sequence()

    num_fractions = 0
    beam_count = 0

    beam_list = trial_info["BeamList"] if trial_info["BeamList"] else []
    if len(beam_list) == 0:
        plan.logger.warning("No Beams found in Trial. Unable to generate RTPLAN.")
        raise MissingTrialBeamsError("No Beams found in Trial.")

    # =======================================================================
    # BEAM LOOP
    # =======================================================================
    for beam in beam_list:
        beam_count += 1
        plan.logger.info("Exporting Plan for beam: %s", beam["Name"])

        # Meterset weights are per-beam: reset here so beam N never inherits
        # cumulative weights parsed from an earlier beam.
        metersetweight = ["0"]

        # --- Patient Setup (one per beam) ---
        patient_setup = _new_dataset()
        patient_setup.PatientPosition = patient_position
        patient_setup.PatientSetupNumber = beam_count
        ds.PatientSetupSequence.append(patient_setup)

        # --- Referenced Beam in Fraction Group ---
        ref_beam = _new_dataset()
        ref_beam.ReferencedBeamNumber = beam_count
        fraction_group.ReferencedBeamSequence.append(ref_beam)

        # --- Beam dataset ---
        beam_ds = _new_dataset()
        ds.BeamSequence.append(beam_ds)

        beam_ds.Manufacturer = ds.Manufacturer  # Consistent with plan-level stamp
        beam_ds.BeamNumber = beam_count
        beam_ds.TreatmentDeliveryType = "TREATMENT"
        beam_ds.ReferencedPatientSetupNumber = beam_count
        # SourceAxisDistance is overridden below from the
        # Pinnacle machine data once the beam's machine has been resolved;
        # 1000 mm remains only as a logged fallback.
        beam_ds.SourceAxisDistance = "1000"
        beam_ds.PrimaryDosimeterUnit = "MU"

        # Primary Fluence Mode
        beam_ds.PrimaryFluenceModeSequence = _new_sequence()
        fluence = _new_dataset()
        fluence.FluenceMode = "STANDARD"
        beam_ds.PrimaryFluenceModeSequence.append(fluence)

        beam_ds.BeamName = beam["FieldID"]
        beam_ds.BeamDescription = beam["Name"]

        # Radiation type
        modality = beam["Modality"]
        if "Photons" in modality:
            beam_ds.RadiationType = "PHOTON"
        elif "Electrons" in modality:
            beam_ds.RadiationType = "ELECTRON"
        else:
            plan.logger.warning(
                "Beam '%s': unrecognised modality '%s'; RadiationType left "
                "empty (proton/ion plans are not yet supported).",
                beam["Name"],
                modality,
            )
            beam_ds.RadiationType = ""

        # Beam type
        if "STATIC" in beam["SetBeamType"].upper():
            beam_ds.BeamType = beam["SetBeamType"].upper()
        else:
            beam_ds.BeamType = "DYNAMIC"

        beam_ds.TreatmentMachineName = beam["MachineNameAndVersion"].partition(":")[0]

        # --- Isocenter (per beam) ---
        # Resolved from the trial's IsocenterName where available; raises
        # IsocenterNotFoundError (failing this trial's RTPLAN) rather than
        # assuming an arbitrary point.
        beam_iso_center = _resolve_beam_isocenter(plan, beam)

        # --- Dose Reference Point ---
        doserefpt = None
        for point in plan.points:
            if point["Name"] == beam["PrescriptionPointName"]:
                doserefpt = plan.convert_point(point)
                plan.logger.debug("Dose reference point found: %s", point["Name"])

        if not doserefpt:
            plan.logger.debug("No dose reference point, setting to isocenter")
            doserefpt = beam_iso_center

        plan.logger.debug("Dose reference point: %s", doserefpt)
        ref_beam.BeamDoseSpecificationPoint = doserefpt

        # --- Control Point Manager ---
        beam_ds.ControlPointSequence = _new_sequence()

        cp_manager = beam["CPManager"]
        if "CPManagerObject" in cp_manager:
            cp_manager = cp_manager["CPManagerObject"]

        numctrlpts = cp_manager["NumberOfControlPoints"]
        plan.logger.debug("Number of control points: %s", numctrlpts)

        # --- Parse control point data from Pinnacle ---
        # Every Pinnacle control point is parsed into its own entry so the
        # builders can give each DICOM control point its own jaw and MLC
        # positions and mechanical angles — previously only the last CP's
        # values survived the loop and were reused for every DICOM CP.
        cp_data_list = []
        for cp_data in cp_manager["ControlPointList"]:
            metersetweight.append(cp_data["Weight"])

            leafpositions, p_count = _parse_mlc_leaf_positions(cp_data)

            cp_data_list.append(
                {
                    "x1": -cp_data["LeftJawPosition"] * 10,
                    "x2": cp_data["RightJawPosition"] * 10,
                    "y1": -cp_data["BottomJawPosition"] * 10,
                    "y2": cp_data["TopJawPosition"] * 10,
                    "leafpositions": leafpositions,
                    "p_count": p_count,
                    "gantry": cp_data["Gantry"],
                    "collimator": cp_data["Collimator"],
                    "couch": cp_data["Couch"],
                }
            )

        if not cp_data_list:
            plan.logger.warning(
                "Beam '%s' has no control points; skipping beam.", beam["Name"]
            )
            raise MissingTrialBeamsError(
                f"Beam '{beam['Name']}' has no control points."
            )

        # Wedge context is constant across a beam's control points in
        # Pinnacle; read it from the first CP.
        wedge_info = _parse_wedge_info(
            cp_manager["ControlPointList"][0], plan.logger
        )
        p_count = cp_data_list[0]["p_count"]

        numwedges = wedge_info["count"] if wedge_info else 0

        # --- Prescription and energy ---
        prescription = [
            p
            for p in trial_info["PrescriptionList"]
            if p["Name"] == beam["PrescriptionName"]
        ][0]

        mnv = beam["MachineNameAndVersion"]
        if ": " in mnv:
            machinename, machineversion = mnv.split(": ", 1)
        else:
            plan.logger.warning(
                "Beam '%s': MachineNameAndVersion '%s' is not in the expected "
                "'name: version' form; version treated as empty.",
                beam["Name"],
                mnv,
            )
            machinename, machineversion = mnv, ""
        machineenergyname = beam["MachineEnergyName"]

        # --- Machine data ---
        # plan.machine_info may be None when plan.Pinnacle.Machines is
        # missing or unparseable; _select_machine handles that (previously
        # machine_info["Name"] would raise TypeError here).
        machine = _select_machine(
            machine_info, machinename, machineversion, plan.logger
        )
        if machine is None:
            plan.logger.warning(
                "Beam '%s': no machine entry matching '%s' (version '%s'); "
                "SAD falls back to 1000 mm and, if this beam uses an MLC, "
                "the trial's RTPLAN export will fail (no assumed "
                "leaf-boundary table is used). See the preceding message "
                "for the machine names present in the file.",
                beam["Name"],
                machinename,
                machineversion,
            )

        # SourceAxisDistance from the machine geometry (fallback 1000 mm).
        sad_mm = _get_machine_sad_mm(machine, plan.logger, beam["Name"])
        if sad_mm is not None:
            beam_ds.SourceAxisDistance = _format_ds(sad_mm)
            plan.logger.debug(
                "Beam '%s': SourceAxisDistance %s mm read from machine data.",
                beam["Name"],
                beam_ds.SourceAxisDistance,
            )
        else:
            plan.logger.warning(
                "Beam '%s': SourceAxisDistance not found in machine data; "
                "assuming standard 1000 mm.",
                beam["Name"],
            )

        # MLC LeafPositionBoundaries from the machine MultiLeaf layout.
        machine_boundaries = _leaf_boundaries_from_machine(
            machine, p_count // 2, plan.logger
        )
        if machine_boundaries is not None:
            plan.logger.debug(
                "Beam '%s': LeafPositionBoundaries derived from machine data "
                "(%d boundaries).",
                beam["Name"],
                len(machine_boundaries),
            )

        energy_matches = re.findall(r"[-+]?\d*\.\d+|\d+", machineenergyname)
        if energy_matches:
            beam_energy = energy_matches[0]
        else:
            plan.logger.warning(
                "Beam '%s': could not parse a numeric energy from '%s'; "
                "defaulting NominalBeamEnergy to 0.",
                beam["Name"],
                machineenergyname,
            )
            beam_energy = "0"

        # Find DosePerMuAtCalibration from machine data
        dose_per_mu_at_cal = -1
        if machine is not None:
            for energy in machine.get("PhotonEnergyList", []) or []:
                if energy.get("Name") == machineenergyname:
                    dose_per_mu_at_cal = energy["PhysicsData"]["OutputFactor"][
                        "DosePerMuAtCalibration"
                    ]
                    plan.logger.debug(
                        "Using DosePerMuAtCalibration of: %s", dose_per_mu_at_cal
                    )

        prescripdose = beam["MonitorUnitInfo"]["PrescriptionDose"]
        normdose = beam["MonitorUnitInfo"]["NormalizedDose"]

        if normdose == 0:
            ref_beam.BeamMeterset = 0
        elif dose_per_mu_at_cal <= 0:
            # No valid calibration was located (machine/energy mismatch, or a
            # non-positive value). Computing prescripdose / (normdose *
            # dose_per_mu_at_cal) here would produce a negative meterset or a
            # division by zero, so leave the Type-3 BeamMeterset/BeamDose unset
            # and warn instead of emitting an invalid value.
            plan.logger.warning(
                "Beam '%s': no valid DosePerMuAtCalibration found (machine "
                "'%s', version '%s', energy '%s'); BeamMeterset and BeamDose "
                "left unset to avoid an invalid value.",
                beam["Name"],
                machinename,
                machineversion,
                machineenergyname,
            )
        else:
            ref_beam.BeamDose = prescripdose / 100
            ref_beam.BeamMeterset = prescripdose / (normdose * dose_per_mu_at_cal)

        # Gantry rotation direction
        is_ccw = cp_manager.get("GantryIsCCW") == 1
        is_cw = cp_manager.get("GantryIsCW") == 1
        if is_ccw and is_cw:
            plan.logger.warning(
                "Beam '%s': both GantryIsCCW and GantryIsCW are set; "
                "defaulting GantryRotationDirection to CW.",
                beam["Name"],
            )
        gantryrotdir = "NONE"
        if is_ccw:
            gantryrotdir = "CC"
        if is_cw:
            gantryrotdir = "CW"

        plan.logger.debug("Beam MU: %s", getattr(ref_beam, "BeamMeterset", None))

        doserate = beam.get("DoseRate", 0)

        # ===================================================================
        # Branch: Step & Shoot vs. Non-Step-and-Shoot
        # ===================================================================
        is_step_and_shoot = (
            "STEP" in beam["SetBeamType"].upper()
            and "SHOOT" in beam["SetBeamType"].upper()
        )

        if is_step_and_shoot:
            _build_step_and_shoot_control_points(
                beam_ds,
                beam,
                plan,
                numctrlpts,
                metersetweight,
                beam_energy,
                doserate,
                gantryrotdir,
                numwedges,
                wedge_info,
                cp_data_list,
                p_count,
                machine_boundaries,
                beam_iso_center,
            )
        else:
            _build_non_ss_control_points(
                beam_ds,
                beam,
                plan,
                numctrlpts,
                metersetweight,
                beam_energy,
                doserate,
                gantryrotdir,
                numwedges,
                wedge_info,
                cp_data_list,
                p_count,
                machine_boundaries,
                beam_iso_center,
            )

        num_fractions = prescription["NumberOfFractions"]
        numwedges = 0  # Reset for next beam

    # --- Fraction Group summary ---
    fraction_group.FractionGroupNumber = 1
    fraction_group.NumberOfFractionsPlanned = num_fractions
    fraction_group.NumberOfBeams = beam_count
    fraction_group.NumberOfBrachyApplicationSetups = "0"

    # --- Save ---
    output_file = os.path.join(export_path, rp_filename)
    plan.logger.info("Creating Plan file: %s", output_file)
    ds.save_as(output_file, enforce_file_format=True)


# ---------------------------------------------------------------------------
# Step & Shoot control point builder
# ---------------------------------------------------------------------------


def _build_step_and_shoot_control_points(
    beam_ds,
    beam,
    plan,
    numctrlpts,
    metersetweight,
    beam_energy,
    doserate,
    gantryrotdir,
    numwedges,
    wedge_info,
    cp_data_list,
    p_count,
    machine_boundaries,
    iso_center,
):
    """Build control points for a Step & Shoot beam.

    Pinnacle control point *i* maps to DICOM control points
    ``2i`` and ``2i+1`` (the segment's dose is delivered between the
    pair), and each pair carries that segment's own jaw and MLC
    positions — previously every control point reused the last
    segment's aperture.

    FinalCumulativeMetersetWeight is set to the last control
    point's cumulative weight instead of a hardcoded "1".
    """
    plan.logger.debug("Using Step & Shoot")

    total_cps = numctrlpts * 2
    beam_ds.NumberOfControlPoints = total_cps
    beam_ds.SourceToSurfaceDistance = beam["SSD"] * 10

    if numwedges > 0:
        beam_ds.WedgeSequence = _create_wedge_sequence(wedge_info)

    # --- Pass 1: cumulative meterset weights per DICOM CP ---
    # Odd control points carry the segment's meterset weight increment.
    cumulative_weights = []
    currentmeterset = 0.0
    metercount = 1
    for j in range(total_cps):
        if j % 2 == 1:
            increment = float(metersetweight[metercount])
            if increment < 0:
                plan.logger.warning(
                    "Beam '%s': negative meterset weight (%s) at Pinnacle "
                    "control point %d; cumulative weights will not be "
                    "monotonic. Verify the plan data.",
                    beam["Name"],
                    increment,
                    metercount,
                )
            currentmeterset += increment
            metercount += 1
        cumulative_weights.append(currentmeterset)

    final_weight = cumulative_weights[-1] if cumulative_weights else 0.0
    beam_ds.FinalCumulativeMetersetWeight = _format_ds(final_weight)
    if abs(final_weight - 1.0) > 1e-3:
        plan.logger.warning(
            "Beam '%s': cumulative meterset weights sum to %s (expected "
            "~1.0). FinalCumulativeMetersetWeight is set to the actual "
            "total, which keeps the plan DICOM-conformant, but verify the "
            "Pinnacle weights.",
            beam["Name"],
            final_weight,
        )

    # --- Pass 2: build the control points ---
    for j in range(total_cps):
        cp = _new_dataset()
        beam_ds.ControlPointSequence.append(cp)

        cp.ControlPointIndex = j
        cp.ReferencedDoseReferenceSequence = _new_sequence()

        dose_ref = _new_dataset()
        cp.ReferencedDoseReferenceSequence.append(dose_ref)

        cmw = cumulative_weights[j]
        cp.CumulativeMetersetWeight = _format_ds(cmw)
        # Coefficient is the delivered fraction: CMW / FCMW.
        coefficient = cmw / final_weight if final_weight else 0.0
        dose_ref.CumulativeDoseReferenceCoefficient = _format_ds(coefficient)
        dose_ref.ReferencedDoseReferenceNumber = "1"

        # DICOM CP j belongs to Pinnacle segment j // 2.
        cp_entry = cp_data_list[min(j // 2, len(cp_data_list) - 1)]

        if j == 0:
            # First control point: all attributes must be present (DICOM C.8.8.14.5)
            _populate_first_control_point(
                cp,
                beam_ds,
                beam,
                plan,
                beam_energy,
                doserate,
                gantryrotdir,
                numwedges,
                cp_entry,
                iso_center,
            )
        else:
            # Subsequent control points: this segment's own jaw and MLC
            # positions (jaws can change between segments too).
            cp.BeamLimitingDevicePositionSequence = _create_bld_position_entries(
                cp_entry["x1"],
                cp_entry["x2"],
                cp_entry["y1"],
                cp_entry["y2"],
                cp_entry["leafpositions"],
            )

    # Beam Limiting Device Sequence (beam level)
    _populate_beam_limiting_device_seq(
        beam_ds, p_count, plan.logger, machine_boundaries
    )


# ---------------------------------------------------------------------------
# Non-Step-and-Shoot (conformal arc / dynamic) control point builder
# ---------------------------------------------------------------------------


def _build_non_ss_control_points(
    beam_ds,
    beam,
    plan,
    numctrlpts,
    metersetweight,
    beam_energy,
    doserate,
    gantryrotdir,
    numwedges,
    wedge_info,
    cp_data_list,
    p_count,
    machine_boundaries,
    iso_center,
):
    """Build control points for a non-Step-and-Shoot beam (e.g. conformal arc).

    Each DICOM control point carries its own jaw, MLC and gantry
    values from the matching Pinnacle control point (the final, appended
    control point repeats the last Pinnacle aperture).

    FinalCumulativeMetersetWeight equals the last control point's
    cumulative weight instead of a hardcoded "1".

    GantryAngle is written on every control point and
    GantryRotationDirection is derived per control point from the angle
    deltas, so beams that reverse direction mid-delivery are represented.
    """
    plan.logger.debug("Not using Step & Shoot")

    total_cps = numctrlpts + 1
    beam_ds.NumberOfControlPoints = total_cps
    beam_ds.SourceToSurfaceDistance = beam["SSD"] * 10

    if numwedges > 0:
        beam_ds.WedgeSequence = _create_wedge_sequence(wedge_info)

    # --- Cumulative meterset weights --------
    cumulative_weights = []
    running = 0.0
    for j in range(total_cps):
        if j > 0:
            running += float(metersetweight[j])
        cumulative_weights.append(running)

    final_weight = cumulative_weights[-1] if cumulative_weights else 0.0
    beam_ds.FinalCumulativeMetersetWeight = _format_ds(final_weight)
    if abs(final_weight - 1.0) > 1e-3:
        plan.logger.warning(
            "Beam '%s': final cumulative meterset weight is %s (expected "
            "~1.0). FinalCumulativeMetersetWeight is set to the actual "
            "value, which keeps the plan DICOM-conformant, but verify the "
            "Pinnacle weights.",
            beam["Name"],
            final_weight,
        )

    for j in range(total_cps):
        cp = _new_dataset()
        beam_ds.ControlPointSequence.append(cp)

        cp.ControlPointIndex = j
        cp.ReferencedDoseReferenceSequence = _new_sequence()

        dose_ref = _new_dataset()
        cp.ReferencedDoseReferenceSequence.append(dose_ref)

        cmw = cumulative_weights[j]
        cp.CumulativeMetersetWeight = _format_ds(cmw)
        coefficient = cmw / final_weight if final_weight else 0.0
        dose_ref.CumulativeDoseReferenceCoefficient = _format_ds(coefficient)
        dose_ref.ReferencedDoseReferenceNumber = "1"

        # appended final CP repeats the last Pinnacle aperture.
        pinn_idx = min(j, len(cp_data_list) - 1)
        cp_entry = cp_data_list[pinn_idx]

        if j == 0:
            # First control point: all attributes must be present
            _populate_first_control_point(
                cp,
                beam_ds,
                beam,
                plan,
                beam_energy,
                doserate,
                gantryrotdir,
                numwedges,
                cp_entry,
                iso_center,
            )
        else:
            # Subsequent control points: this control point's own jaw and
            # MLC positions.
            cp.BeamLimitingDevicePositionSequence = _create_bld_position_entries(
                cp_entry["x1"],
                cp_entry["x2"],
                cp_entry["y1"],
                cp_entry["y2"],
                cp_entry["leafpositions"],
            )

            # per-control-point gantry angle and rotation
            # direction so direction reversals mid-delivery are captured.
            prev_entry = cp_data_list[min(j - 1, len(cp_data_list) - 1)]
            cp.GantryAngle = cp_entry["gantry"]
            cp.GantryRotationDirection = _gantry_direction_between(
                prev_entry["gantry"], cp_entry["gantry"]
            )

    # Beam Limiting Device Sequence (beam level)
    _populate_beam_limiting_device_seq(
        beam_ds, p_count, plan.logger, machine_boundaries
    )