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
from .pinnacle_metadata import append_pinnacle_metadata_for_plan, apply_equipment_stamps


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


def _generate_leaf_position_boundaries():
    """Generate the standard 51-element Varian-style MLC leaf position boundaries.

    Boundaries go from -200 to +200 mm in a pattern of:
      -200 to -100 in 10 mm steps (11 values)
      -95  to  +95 in  5 mm steps (39 values, but -100 and +100 already counted)
      +100 to +200 in 10 mm steps
    Total: 51 boundaries for 50 leaf pairs.
    """
    boundaries = []
    # Outer leaves: -200 to -100 in steps of 10
    for v in range(-200, -100, 10):
        boundaries.append(str(v))
    # Inner leaves: -100 to +100 in steps of 5
    for v in range(-100, 105, 5):
        boundaries.append(str(v))
    # Outer leaves: +110 to +200 in steps of 10
    for v in range(110, 210, 10):
        boundaries.append(str(v))
    return boundaries


# Pre-compute the leaf boundaries so we don't rebuild the list on every beam.
_LEAF_POSITION_BOUNDARIES = _generate_leaf_position_boundaries()


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
            info["orientation"] = "0"  # TODO: confirm orientation mapping
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


def _populate_beam_limiting_device_seq(beam_ds, p_count, logger=None):
    """Populate the BeamLimitingDeviceSequence for a beam.

    Creates entries for ASYMX, ASYMY, and MLCX.
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

    mlcx = _new_dataset()
    mlcx.RTBeamLimitingDeviceType = "MLCX"
    # NumberOfLeafJawPairs has VR=IS (integer); use integer division so we
    # never emit a fractional value such as "60.0".
    num_pairs = p_count // 2
    mlcx.NumberOfLeafJawPairs = num_pairs
    mlcx.LeafPositionBoundaries = _LEAF_POSITION_BOUNDARIES
    # DICOM requires len(LeafPositionBoundaries) == NumberOfLeafJawPairs + 1.
    # The boundary table is a fixed Varian-Millennium pattern, so warn when the
    # data implies a different MLC (e.g. Elekta 80-pair, Varian HD120).
    if logger is not None and num_pairs + 1 != len(_LEAF_POSITION_BOUNDARIES):
        logger.warning(
            "MLC leaf-pair count from the plan data (%d) does not match the "
            "hardcoded LeafPositionBoundaries table (%d boundaries for %d "
            "pairs); the exported MLC geometry may be incorrect for this "
            "machine.",
            num_pairs,
            len(_LEAF_POSITION_BOUNDARIES),
            len(_LEAF_POSITION_BOUNDARIES) - 1,
        )
    beam_ds.BeamLimitingDeviceSequence.append(mlcx)


def _create_bld_position_entries(x1, x2, y1, y2, leafpositions):
    """Create the three BeamLimitingDevicePositionSequence items for CP 0.

    Returns a Sequence containing ASYMX, ASYMY, and MLCX entries.
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

    mlcx = _new_dataset()
    mlcx.RTBeamLimitingDeviceType = "MLCX"
    mlcx.LeafJawPositions = leafpositions
    bld_seq.append(mlcx)

    return bld_seq


def _create_mlc_only_position_entry(leafpositions):
    """Create a single-item BeamLimitingDevicePositionSequence for MLC only.

    Used for control points after the first, where only MLC positions change.
    """
    bld_seq = _new_sequence()
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
    gantryangle,
    colangle,
    psupportangle,
    gantryrotdir,
    numwedges,
    x1,
    x2,
    y1,
    y2,
    leafpositions,
):
    """Populate all attributes required by DICOM for the first control point.

    Per DICOM C.8.8.14.5: at the first control point, ALL applicable
    attributes must be present. This includes 1C and 2C attributes.
    """
    # --- Required energy and dose rate (Type 3 but universally expected) ---
    cp.NominalBeamEnergy = beam_energy
    cp.DoseRateSet = doserate

    # --- Gantry (1C — required at first CP) ---
    cp.GantryAngle = gantryangle
    cp.GantryRotationDirection = gantryrotdir

    # --- Collimator (1C — required at first CP) ---
    cp.BeamLimitingDeviceAngle = colangle
    cp.BeamLimitingDeviceRotationDirection = "NONE"

    # --- Patient Support / Couch (1C — required at first CP) ---
    cp.PatientSupportAngle = psupportangle
    cp.PatientSupportRotationDirection = "NONE"

    # --- Table Top Eccentric (1C — required at first CP) ---
    cp.TableTopEccentricAngle = "0"
    cp.TableTopEccentricRotationDirection = "NONE"

    # --- Table Top Position (2C — required at first CP, may be empty) ---
    cp.TableTopVerticalPosition = ""
    cp.TableTopLongitudinalPosition = ""
    cp.TableTopLateralPosition = ""

    # --- Isocenter (2C — required at first CP when isocentric) ---
    cp.IsocenterPosition = plan.iso_center

    # --- Source to Surface Distance (Type 3) ---
    cp.SourceToSurfaceDistance = beam["SSD"] * 10

    # --- Wedge position (1C — required when wedges present) ---
    if numwedges > 0:
        cp.WedgePositionSequence = _create_wedge_position_seq()

    # --- Beam Limiting Device positions (1C) ---
    cp.BeamLimitingDevicePositionSequence = _create_bld_position_entries(
        x1, x2, y1, y2, leafpositions
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

    # TODO Fix the RTPLAN export functionality and remove this warning
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

    # --- Plan identification ---
    ds.RTPlanLabel = f"{plan_info['PlanName']}.0"
    ds.RTPlanName = plan_info["PlanName"]
    ds.RTPlanDescription = append_pinnacle_metadata_for_plan(
        None,
        plan,
        trial_info,
        max_length=1024,
    )
    ds.RTPlanDate = ds.StudyDate
    ds.RTPlanTime = ds.StudyTime
    ds.PlanIntent = ""  # TODO: palliative/curative — find source
    ds.RTPlanGeometry = "PATIENT"

    # --- Referenced Structure Set ---
    ds.ReferencedStructureSetSequence = _new_sequence()
    ref_struct = _new_dataset()
    ref_struct.ReferencedSOPClassUID = RTStructSOPClassUID
    ref_struct.ReferencedSOPInstanceUID = struct_instance_uid
    ds.ReferencedStructureSetSequence.append(ref_struct)

    ds.ApprovalStatus = "UNAPPROVED"  # TODO: derive from trial file

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
        beam_ds.SourceAxisDistance = "1000"
        beam_ds.FinalCumulativeMetersetWeight = "1"
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

        # --- Dose Reference Point ---
        doserefpt = None
        for point in plan.points:
            if point["Name"] == beam["PrescriptionPointName"]:
                doserefpt = plan.convert_point(point)
                plan.logger.debug("Dose reference point found: %s", point["Name"])

        if not doserefpt:
            plan.logger.debug("No dose reference point, setting to isocenter")
            doserefpt = plan.iso_center

        plan.logger.debug("Dose reference point: %s", doserefpt)
        ref_beam.BeamDoseSpecificationPoint = doserefpt

        # --- Control Point Manager ---
        beam_ds.ControlPointSequence = _new_sequence()

        cp_manager = beam["CPManager"]
        if "CPManagerObject" in cp_manager:
            cp_manager = cp_manager["CPManagerObject"]

        numctrlpts = cp_manager["NumberOfControlPoints"]
        currentmeterset = 0.0
        plan.logger.debug("Number of control points: %s", numctrlpts)

        # --- Parse control point data from Pinnacle ---
        # Extract jaw positions (from first CP that has them) and all leaf positions
        x1 = x2 = y1 = y2 = None
        leafpositions = []
        p_count = 0
        gantryangle = colangle = psupportangle = 0
        wedge_info = None

        for cp_data in cp_manager["ControlPointList"]:
            metersetweight.append(cp_data["Weight"])

            # Jaw positions — keep only the first values encountered
            if x1 is None:
                x1 = -cp_data["LeftJawPosition"] * 10
            if x2 is None:
                x2 = cp_data["RightJawPosition"] * 10
            if y2 is None:
                y2 = cp_data["TopJawPosition"] * 10
            if y1 is None:
                y1 = -cp_data["BottomJawPosition"] * 10

            # MLC leaf positions
            leafpositions, p_count = _parse_mlc_leaf_positions(cp_data)

            # Mechanical angles
            gantryangle = cp_data["Gantry"]
            colangle = cp_data["Collimator"]
            psupportangle = cp_data["Couch"]

            # Wedge
            wedge_info = _parse_wedge_info(cp_data, plan.logger)

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
        if (
            machine_info["Name"] == machinename
            and machine_info["VersionTimestamp"] == machineversion
        ):
            for energy in machine_info["PhotonEnergyList"]:
                if energy["Name"] == machineenergyname:
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
                gantryangle,
                colangle,
                psupportangle,
                gantryrotdir,
                numwedges,
                wedge_info,
                x1,
                x2,
                y1,
                y2,
                leafpositions,
                p_count,
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
                gantryangle,
                colangle,
                psupportangle,
                gantryrotdir,
                numwedges,
                wedge_info,
                x1,
                x2,
                y1,
                y2,
                leafpositions,
                p_count,
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
    gantryangle,
    colangle,
    psupportangle,
    gantryrotdir,
    numwedges,
    wedge_info,
    x1,
    x2,
    y1,
    y2,
    leafpositions,
    p_count,
):
    """Build control points for a Step & Shoot beam."""
    plan.logger.debug("Using Step & Shoot")

    total_cps = numctrlpts * 2
    beam_ds.NumberOfControlPoints = total_cps
    beam_ds.SourceToSurfaceDistance = beam["SSD"] * 10

    if numwedges > 0:
        beam_ds.WedgeSequence = _create_wedge_sequence(wedge_info)

    metercount = 1
    currentmeterset = 0.0

    for j in range(total_cps):
        cp = _new_dataset()
        beam_ds.ControlPointSequence.append(cp)

        cp.ControlPointIndex = j
        cp.BeamLimitingDevicePositionSequence = _new_sequence()
        cp.ReferencedDoseReferenceSequence = _new_sequence()

        dose_ref = _new_dataset()
        cp.ReferencedDoseReferenceSequence.append(dose_ref)

        if j % 2 == 1:  # odd control points carry the meterset weight
            currentmeterset += float(metersetweight[metercount])
            metercount += 1

        cp.CumulativeMetersetWeight = currentmeterset
        dose_ref.CumulativeDoseReferenceCoefficient = currentmeterset
        dose_ref.ReferencedDoseReferenceNumber = "1"

        if j == 0:
            # First control point: all attributes must be present (DICOM C.8.8.14.5)
            _populate_first_control_point(
                cp,
                beam_ds,
                beam,
                plan,
                beam_energy,
                doserate,
                gantryangle,
                colangle,
                psupportangle,
                gantryrotdir,
                numwedges,
                x1,
                x2,
                y1,
                y2,
                leafpositions,
            )
        else:
            # Subsequent control points: only MLC changes
            cp.BeamLimitingDevicePositionSequence = _create_mlc_only_position_entry(
                leafpositions
            )

    # Beam Limiting Device Sequence (beam level)
    _populate_beam_limiting_device_seq(beam_ds, p_count, plan.logger)


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
    gantryangle,
    colangle,
    psupportangle,
    gantryrotdir,
    numwedges,
    wedge_info,
    x1,
    x2,
    y1,
    y2,
    leafpositions,
    p_count,
):
    """Build control points for a non-Step-and-Shoot beam (e.g. conformal arc)."""
    plan.logger.debug("Not using Step & Shoot")

    total_cps = numctrlpts + 1
    beam_ds.NumberOfControlPoints = total_cps
    beam_ds.SourceToSurfaceDistance = beam["SSD"] * 10

    if numwedges > 0:
        beam_ds.WedgeSequence = _create_wedge_sequence(wedge_info)

    for j in range(total_cps):
        cp = _new_dataset()
        beam_ds.ControlPointSequence.append(cp)

        cp.ControlPointIndex = j
        cp.BeamLimitingDevicePositionSequence = _new_sequence()
        cp.ReferencedDoseReferenceSequence = _new_sequence()

        dose_ref = _new_dataset()
        cp.ReferencedDoseReferenceSequence.append(dose_ref)

        cp.CumulativeMetersetWeight = metersetweight[j]

        if j == 0:
            # First control point: all attributes must be present
            dose_ref.CumulativeDoseReferenceCoefficient = "0"
            dose_ref.ReferencedDoseReferenceNumber = "1"

            _populate_first_control_point(
                cp,
                beam_ds,
                beam,
                plan,
                beam_energy,
                doserate,
                gantryangle,
                colangle,
                psupportangle,
                gantryrotdir,
                numwedges,
                x1,
                x2,
                y1,
                y2,
                leafpositions,
            )
        else:
            # Subsequent control points: only MLC
            cp.BeamLimitingDevicePositionSequence = _create_mlc_only_position_entry(
                leafpositions
            )
            dose_ref.CumulativeDoseReferenceCoefficient = "1"
            dose_ref.ReferencedDoseReferenceNumber = "1"

    # Beam Limiting Device Sequence (beam level)
    _populate_beam_limiting_device_seq(beam_ds, p_count, plan.logger)
