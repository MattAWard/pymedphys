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

from pymedphys._imports import pydicom

from .pinn_yaml import pinn_to_dict
from .rtstruct import find_iso_center

# Medical Connections offers free valid UIDs
# (http://www.medicalconnections.co.uk/FreeUID.html)
# Their service was used to obtain the following root UID for this tool:
UID_PREFIX = "1.2.826.0.1.3680043.10.202."


class PinnaclePlan:
    """Represents a plan within the Pinnacle data.

    This class manages the data specific to a plan within a Pinnacle dataset.

    Parameters
    ----------
        pinnacle : PinnacleExport
            PinnacleExport object representing the dataset.
        path : str
            Path to raw Pinnacle data (directoy containing 'Patient' file).
        plan : dict
            Plan info dict from 'Patient' file.
    """

    def __init__(self, pinnacle, path, plan):
        self._pinnacle = pinnacle
        self._path = path
        self._plan_info = plan

        self._machine_info = None  # Data of the machines used for this plan
        self._trials = None  # Data found in plan.Trial
        self._trial_info = None  # 'Active' trial for this plan
        self._points = None  # Data found in plan.Points
        self._patient_setup = None  # Data found in PatientSetup file
        self._primary_image = None  # Primary image for this plan
        self._plan_info_file = None  # Data found in plan.PlanInfo (per-plan)
        self._uid_prefix = UID_PREFIX  # The prefix for UIDs generated

        self._roi_count = 0  # Store the total number of ROIs
        self._iso_center = []  # Store the iso center of this plan
        self._ct_center = []  # Store the ct center of the primary image of this plan
        self._dose_ref_pt = []  # Store the dose ref point of this plan

        self._plan_inst_uid = None  # UID for RTPlan instance
        self._dose_inst_uid = None  # UID for RTDose instance
        self._struct_inst_uid = None  # UID for RTStruct instance
        self._trial_uid_cache = {}  # Cache: trial key → uid dict

        for image in pinnacle.images:
            if image.image["ImageSetID"] == self.plan_info["PrimaryCTImageSetID"]:
                self._primary_image = image

        if not self._primary_image:
            self.logger.warning("Primary Image Not Available")

    @property
    def logger(self):
        """Gets the configured logger.

        Returns
        -------
        logger : Logger
            Logger configured.
        """
        return self._pinnacle.logger

    @property
    def pinnacle(self):
        """Gets the PinnacleExport object.

        Returns
        -------
        pinnacle : PinnacleExport
            PinnacleExport object for this dataset.
        """
        return self._pinnacle

    @property
    def path(self):
        """Gets the path of the Pinnacle data.

        Returns
        -------
        path : str
            Path containing the Pinnacle data.
        """
        return self._path

    @property
    def primary_image(self):
        """Gets the primary image for this plan.

        Returns
        -------
        primary_image : PinnacleImage
            PinnacleImage representing the primary image for this plan.
        """
        return self._primary_image

    @property
    def machine_info(self):
        """Gets the machine info for this plan.

        Returns
        -------
        machine_info : dict
            Machine info read from 'plan.Pinnacle.Machines' file.
        """

        if not self._machine_info:
            path_machine = os.path.join(self._path, "plan.Pinnacle.Machines")
            self.logger.debug("Reading machine data from: %s", path_machine)
            self._machine_info = pinn_to_dict(path_machine)

        return self._machine_info

    @property
    def trials(self):
        """Gets the trials within this plan.

        Returns
        -------
        trials : list
            List of all trials found within this plan.
        """

        if not self._trials:
            path_trial = os.path.join(self._path, "plan.Trial")
            self.logger.debug("Reading trial data from: %s", path_trial)
            self._trials = pinn_to_dict(path_trial)
            if isinstance(self._trials, dict):
                self._trials = [
                    self._trials["Trial"]
                ]  # make sure trials is always a list of dicts

            # trial info is always a dict (we can ignore this)
            if not self._trial_info:
                self._trial_info = self._trials[0]

            self.logger.debug("Number of trials read: %s", len(self._trials))
            self.logger.debug("Active Trial: %s", self._trial_info["Name"])

        return self._trials

    @property
    def active_trial(self):
        """Get and set the active trial for this plan.

        When DICOM objects are exported, data from the active trial is
        used to generate the output.
        """

        return self._trial_info

    @active_trial.setter
    def active_trial(self, trial_name):
        if isinstance(trial_name, str):
            for trial in self.trials:
                if trial["Name"] == trial_name:
                    self._trial_info = trial
                    self.logger.info("Active Trial set: %s", trial_name)
                    return

        raise KeyError

    @property
    def plan_info(self):
        """Gets the plan information for this plan.

        Returns
        -------
        plan_info : dict
            Plan info as found in the Pinnacle 'Patient' file.
        """
        return self._plan_info

    @property
    def plan_info_file(self):
        """Gets the per-plan ``plan.PlanInfo`` file contents.

        This is distinct from :attr:`plan_info` (which is the ``PlanList``
        entry from the top-level ``Patient`` file). The per-plan
        ``Plan_N/plan.PlanInfo`` file carries richer metadata that the
        ``Patient`` file doesn't include — most notably the lock/approval
        audit fields (e.g. ``PlanLockStatus``).

        Returns an empty dict if the file is absent (older Pinnacle
        versions, or partial archives), so callers can use ``.get(...)``
        without further None-checks.

        Returns
        -------
        plan_info_file : dict
            Contents of ``Plan_N/plan.PlanInfo``, or ``{}`` if missing.
        """
        if self._plan_info_file is None:
            path = os.path.join(self._path, "plan.PlanInfo")
            if os.path.exists(path):
                self.logger.debug("Reading plan info file from: %s", path)
                self._plan_info_file = pinn_to_dict(path) or {}
            else:
                self.logger.debug("plan.PlanInfo not found at: %s", path)
                self._plan_info_file = {}
        return self._plan_info_file

    @property
    def trial_info(self):
        """Gets the trial information of the active trial.

        Returns
        -------
        trial_info : dict
            Trial info from the 'plan.Trial' file.
        """

        if not self._trial_info:
            # Ensures that the trials are read and a default trial_info is set
            _ = self.trials

        return self._trial_info

    @property
    def points(self):
        """Gets the points defined within the plan.

        Returns
        -------
        points : list
            List of points read from the 'plan.Points' file.
        """

        if not self._points:
            path_points = os.path.join(self._path, "plan.Points")
            self.logger.debug("Reading points data from: %s", path_points)
            self._points = pinn_to_dict(path_points)

            if isinstance(self._points, dict):
                self._points = [self._points["Poi"]]

            if self._points is None:
                self._points = []

        return self._points

    @property
    def patient_position(self):
        """Gets the patient position

        Returns
        -------
        patient_position : str
            The patient position for this plan.
        """

        if not self._patient_setup:
            self._patient_setup = pinn_to_dict(
                os.path.join(self._path, "plan.PatientSetup")
            )

        pat_pos = ""

        if "Head First" in self._patient_setup["Orientation"]:
            pat_pos = "HF"
        elif "Feet First" in self._patient_setup["Orientation"]:
            pat_pos = "FF"

        if "supine" in self._patient_setup["Position"]:
            pat_pos = f"{pat_pos}S"
        elif "prone" in self._patient_setup["Position"]:
            pat_pos = f"{pat_pos}P"
        elif (
            "decubitus right" in self._patient_setup["Position"]
            or "Decuibitus Right" in self._patient_setup["Position"]
        ):
            pat_pos = f"{pat_pos}DR"
        elif (
            "decubitus left" in self._patient_setup["Position"]
            or "Decuibitus Left" in self._patient_setup["Position"]
        ):
            pat_pos = f"{pat_pos}DL"

        return pat_pos

    @property
    def iso_center(self):
        """Gets and sets the iso center for this plan."""

        if len(self._iso_center) == 0:
            find_iso_center(self)

        return self._iso_center

    @iso_center.setter
    def iso_center(self, iso_center):
        self._iso_center = iso_center

    @staticmethod
    def is_prefix_valid(prefix):
        """Check if a UID prefix is valid.

        Parameters
        ----------
            prefix :  str
                The UID prefix to check.

        Returns:
            True if valid, False otherwise.
        """

        if re.match(pydicom.uid.RE_VALID_UID_PREFIX, prefix):
            return True

        return False

    def generate_uids_for_trial(self, trial_info, uid_type="RANDOM"):
        """Generate (or retrieve cached) RTPLAN, RTDOSE and RTSTRUCT UIDs for
        a specific trial.

        This is the canonical UID source for multi-trial export. It does NOT
        mutate ``self._plan_inst_uid`` / ``self._struct_inst_uid`` so it is
        safe to call repeatedly inside a loop over trials without UID values
        from one trial leaking into another.

        Results are cached per trial (keyed by trial name + write timestamp)
        so that ``convert_struct``, ``convert_plan`` and ``convert_dose`` all
        receive the same UIDs when they independently call this method for the
        same trial.  This ensures cross-references between RT objects remain
        consistent.  Call ``clear_uid_cache()`` to force fresh UIDs on the
        next call.

        Each call produces both SOP Instance UIDs and Series Instance UIDs.
        The Series UIDs are always random (never deterministic) because
        DICOM requires SeriesInstanceUID to be distinct from SOPInstanceUID
        — they live at different levels of the DICOM hierarchy.

        Parameters
        ----------
            trial_info : dict
                The trial dictionary (one entry from ``self.trials``) to
                generate UIDs for. Trial name and write timestamp are mixed
                into the entropy so the same trial always hashes to the
                same UIDs across runs (HASH mode only).
            uid_type : str, optional
                If 'HASH', deterministic entropy-based SOP UIDs are
                generated.  Anything else falls back to pydicom's random
                UIDs.  Default: 'RANDOM'.

                .. note::
                    The default was changed from 'HASH' to 'RANDOM' because
                    deterministic UIDs cause "inconsistent link" rejections
                    when a PACS/Conquest destination already holds objects
                    from a previous export of the same patient.  Callers
                    that genuinely need idempotent UIDs (e.g. offline
                    scripting where the output folder is wiped each time)
                    can still pass ``uid_type='HASH'`` explicitly.

        Returns
        -------
        uids : dict
            Mapping with keys ``"plan"``, ``"dose"``, ``"struct"``,
            ``"series_plan"``, ``"series_dose"`` and ``"series_struct"``.
            The first three are SOP Instance UIDs; the last three are
            Series Instance UIDs (always random).
        """

        # Build a cache key from trial identity so that every caller
        # (convert_struct, convert_plan, convert_dose) that asks for UIDs
        # for the *same* trial gets back the *same* set.  Without this,
        # each independent call generates fresh random UIDs and the
        # cross-references between RT objects are broken.
        cache_key = (
            trial_info["Name"],
            trial_info["ObjectVersion"]["WriteTimeStamp"],
        )
        if cache_key in self._trial_uid_cache:
            self.logger.debug(
                "Trial '%s' UIDs (cached) - plan: %s, dose: %s, struct: %s",
                trial_info["Name"],
                self._trial_uid_cache[cache_key]["plan"],
                self._trial_uid_cache[cache_key]["dose"],
                self._trial_uid_cache[cache_key]["struct"],
            )
            return self._trial_uid_cache[cache_key]

        entropy_srcs = None
        if uid_type == "HASH":
            entropy_srcs = [
                self._pinnacle.patient_info["MedicalRecordNumber"],
                self.plan_info["PlanName"],
                trial_info["Name"],
                trial_info["ObjectVersion"]["WriteTimeStamp"],
            ]

        plan_uid = pydicom.uid.generate_uid(
            prefix=f"{self._uid_prefix}1.", entropy_srcs=entropy_srcs
        )
        dose_uid = pydicom.uid.generate_uid(
            prefix=f"{self._uid_prefix}2.", entropy_srcs=entropy_srcs
        )
        struct_uid = pydicom.uid.generate_uid(
            prefix=f"{self._uid_prefix}3.", entropy_srcs=entropy_srcs
        )

        # Series Instance UIDs are always random — they must never equal
        # the SOP Instance UID (different DICOM hierarchy levels) and
        # deterministic series UIDs offer no benefit.
        series_plan_uid = pydicom.uid.generate_uid()
        series_dose_uid = pydicom.uid.generate_uid()
        series_struct_uid = pydicom.uid.generate_uid()

        self.logger.debug(
            "Trial '%s' UIDs - plan: %s, dose: %s, struct: %s",
            trial_info["Name"],
            plan_uid,
            dose_uid,
            struct_uid,
        )

        uids = {
            "plan": plan_uid,
            "dose": dose_uid,
            "struct": struct_uid,
            "series_plan": series_plan_uid,
            "series_dose": series_dose_uid,
            "series_struct": series_struct_uid,
        }
        self._trial_uid_cache[cache_key] = uids
        return uids

    def clear_uid_cache(self):
        """Discard all cached trial UIDs.

        Call this before re-exporting the same plan if fresh UIDs are
        required (e.g. when writing to a new output directory).
        """
        self._trial_uid_cache.clear()

    def generate_uids(self, uid_type="RANDOM"):
        """Generates UIDs for the *active* trial and caches them on the plan.

        Retained for backwards compatibility with callers (notably the CLI
        and the legacy ``plan_inst_uid`` / ``struct_inst_uid`` properties)
        that expect a single set of UIDs per plan instance.

        Multi-trial export paths should call ``generate_uids_for_trial``
        directly with each trial in turn instead of relying on this method.

        Parameters
        ----------
            uid_type : str, optional
                If 'HASH', the entropy will be generated
                to hash to consistent UIDs. If not then random UIDs will be
                generated. Default: 'RANDOM'
        """

        uids = self.generate_uids_for_trial(self.trial_info, uid_type=uid_type)
        self._plan_inst_uid = uids["plan"]
        self._dose_inst_uid = uids["dose"]
        self._struct_inst_uid = uids["struct"]

    @property
    def plan_inst_uid(self):
        """Gets the instance UID for RTPLAN.

        Returns
        -------
        uid : str
            The UID to use for the plan.
        """

        if not self._plan_inst_uid:
            self.generate_uids()

        return self._plan_inst_uid

    @property
    def struct_inst_uid(self):
        """Gets the instance UID for RTSTRUCT.

        Returns
        -------
        uid : str
            The UID to use for the struct.
        """

        if not self._struct_inst_uid:
            self.generate_uids()

        return self._struct_inst_uid

    # Convert the point from the pinnacle plan format to dicom
    def convert_point(self, point):
        """Convert a point from Pinnacle coordinates to DICOM coordinates.

        Parameters
        ----------
            point : list
                The point to convert.

        Returns
        -------
            The converted point.

        """

        image_header = self.primary_image.image_header

        refpoint = [point["XCoord"] * 10, point["YCoord"] * 10, point["ZCoord"] * 10]
        if (
            image_header["patient_position"] == "HFP"
            or image_header["patient_position"] == "FFS"
        ):
            refpoint[0] = -refpoint[0]
        if (
            image_header["patient_position"] == "HFS"
            or image_header["patient_position"] == "FFS"
        ):
            refpoint[1] = -(refpoint[1])
        if (
            image_header["patient_position"] == "HFS"
            or image_header["patient_position"] == "HFP"
        ):
            refpoint[2] = -(refpoint[2])

        point["refpoint"] = refpoint

        refpoint[0] = round(refpoint[0], 5)
        refpoint[1] = round(refpoint[1], 5)
        refpoint[2] = round(refpoint[2], 5)

        return refpoint
