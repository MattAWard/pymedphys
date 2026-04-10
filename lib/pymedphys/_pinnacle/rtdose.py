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


import math
import os
import re
import struct
import time

from pymedphys._imports import numpy as np
from pymedphys._imports import pydicom
from pymedphys._pinnacle.pinnacle_exceptions import MissingBeamDoseError, MissingCTImageError, MissingTrialBeamsError

from pymedphys._dicom.orientation import IMAGE_ORIENTATION_MAP

from .constants import (
    GImplementationClassUID,
    GTransferSyntaxUID,
    Manufacturer,
    RTDOSEModality,
    RTDoseSOPClassUID,
    RTPlanSOPClassUID,
)

def construct_dose_from_binary(binary_data, array):
    """
    Read binary data into empty dose array
    """
    X, Y, Z = array.shape
    idx=0
    for z in range(Z - 1, -1, -1):
        for y in range(Y):
            for x in range(X):
                data_element = binary_data[idx:idx+4]
                value = struct.unpack(">f", data_element)[0]
                array[x, y, z] = value
                idx += 4
    return array


def read_binary_data(binary_file):
    """
    Check if the supplied binary file is non-empty and return the data if so
    """
    if os.path.isfile(binary_file):
        size = os.path.getsize(binary_file)
        with open(binary_file, "rb") as b:
            data = b.read()
            if all(byte == 0 for byte in data):
                return False
            else:
                return data


def trilinear_interpolation(idx, grid):
    """
    Return trilinear interpolated value for a voxel with index idx within the grid
    """

    int_idx = [math.floor(f) for f in idx]
    frac_idx = [f % 1 for f in idx]

    l1 = [[[0 for x in range(2)] for x in range(2)] for x in range(2)]
    for x in range(0, 2):
        for y in range(0, 2):
            for z in range(0, 2):
                l1[x][y][z] = grid[int_idx[0] + x, int_idx[1] + y, int_idx[2] + z]

    l2 = [[0 for x in range(2)] for x in range(2)]
    for y in range(0, 2):
        for z in range(0, 2):
            l2[y][z] = l1[0][y][z] * (1 - frac_idx[0]) + l1[1][y][z] * frac_idx[0]

    l3 = [0 for x in range(2)]
    for z in range(0, 2):
        l3[z] = l2[0][z] * (1 - frac_idx[1]) + l2[1][z] * frac_idx[1]

    return l3[0] * (1 - frac_idx[2]) + l3[1] * frac_idx[2]


def convert_dose(plan, export_path):
    """Export RTDose files for all trials in the plan.

    For each trial, generates unique UIDs and creates an RTDOSE DICOM file.

    Parameters
    ----------
    plan : PinnaclePlan
        The plan object containing trials and plan information.
    export_path : str
        Directory where DICOM files will be saved.
    """
    # Check that the plan has a primary image
    if not plan.primary_image:
        plan.logger.error("No primary image found for plan. Unable to generate RTDOSE.")
        raise MissingCTImageError("Plan has no primary image associated with it.")

    supported_orientations = ("HFS", "HFP", "FFS", "FFP")

    patient_info = plan.pinnacle.patient_info
    plan_info = plan.plan_info
    image_info = plan.primary_image.image_info[0]
    patient_position = plan.patient_position

    if patient_position not in supported_orientations:
        raise NotImplementedError(
            f"{patient_position} orientation not supported. Only: "
            f"{supported_orientations}"
        )

    # Calculate dose origin shifts (plan-level, same for all trials)
    if patient_position in ("HFP", "FFS"):
        dose_origin_x_sign = -1
    elif patient_position in ("HFS", "FFP"):
        dose_origin_x_sign = 1

    if patient_position in ("HFS", "FFS"):
        dose_origin_y_sign = -1
    elif patient_position in ("HFP", "FFP"):
        dose_origin_y_sign = 1

    if patient_position in ("HFS", "HFP"):
        dose_origin_z_sign = -1
    elif patient_position in ("FFS", "FFP"):
        dose_origin_z_sign = 1

    # Create base DICOM dataset with plan-level metadata (used for all trials)
    file_meta = pydicom.dataset.Dataset()
    file_meta.MediaStorageSOPClassUID = RTDoseSOPClassUID
    file_meta.TransferSyntaxUID = GTransferSyntaxUID
    file_meta.ImplementationClassUID = GImplementationClassUID

    ds = pydicom.dataset.FileDataset(
        "RD.dcm", {}, file_meta=file_meta, preamble=b"\x00" * 128
    )
    ds.SpecificCharacterSet = "ISO_IR 100"
    ds.InstanceCreationDate = time.strftime("%Y%m%d")
    ds.InstanceCreationTime = time.strftime("%H%M%S")

    ds.SOPClassUID = RTDoseSOPClassUID  # RT Dose Storage

    ds.AccessionNumber = ""
    ds.Modality = RTDOSEModality
    ds.Manufacturer = Manufacturer
    ds.OperatorsName = ""
    ds.ManufacturerModelName = plan_info.get("ToolType", "")
    ds.SoftwareVersions = [plan_info["PinnacleVersionDescription"]]
    ds.PhysiciansOfRecord = patient_info["RadiationOncologist"]
    ds.PatientName = patient_info["FullName"]
    ds.PatientBirthDate = patient_info["DOB"]
    ds.PatientID = patient_info["MedicalRecordNumber"]
    ds.PatientSex = patient_info["Gender"][0]

    ds.StudyInstanceUID = image_info["StudyInstanceUID"]
    ds.FrameOfReferenceUID = image_info["FrameUID"]
    ds.StudyID = plan.primary_image.image["StudyID"]

    ds.ImageOrientationPatient = IMAGE_ORIENTATION_MAP[patient_position]
    ds.PositionReferenceIndicator = ""
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"

    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 0
    ds.DoseUnits = "GY"
    ds.DoseType = "PHYSICAL"
    ds.DoseSummationType = "PLAN"

    ds.ReferencedRTPlanSequence = pydicom.sequence.Sequence()
    ds.ReferencedRTPlanSequence.append(pydicom.dataset.Dataset())
    ds.ReferencedRTPlanSequence[0].ReferencedSOPClassUID = RTPlanSOPClassUID

    ds.TissueHeterogeneityCorrection = "IMAGE"

    # Use the same plan UID for all trials, and generate a fresh dose UID each time.
    planInstanceUID = plan.plan_inst_uid
    ds.ReferencedRTPlanSequence[0].ReferencedSOPInstanceUID = planInstanceUID

    # Process each trial
    for trial_info in plan.trials:
        plan.active_trial = trial_info["Name"]
        plan.logger.info("Exporting Dose for trial: %s", trial_info["Name"])

        # Separate dose UID for each trial
        doseInstanceUID = pydicom.uid.generate_uid(
            prefix=f"{plan._uid_prefix}2.",
            entropy_srcs=[
                plan.pinnacle.patient_info["MedicalRecordNumber"],
                plan_info["PlanName"],
                trial_info["Name"],
                trial_info["ObjectVersion"]["WriteTimeStamp"],
            ],
        )

        # Calculate dose origin for this trial
        dose_origin = [
            dose_origin_x_sign * trial_info["DoseGrid .Origin .X"] * 10,
            dose_origin_y_sign * trial_info["DoseGrid .Origin .Y"] * 10,
            dose_origin_z_sign * trial_info["DoseGrid .Origin .Z"] * 10,
        ]

        # Call the trial-specific dose conversion function
        convert_dose_for_trial(
            plan,
            trial_info,
            doseInstanceUID,
            planInstanceUID,
            dose_origin,
            patient_position,
            ds,
            export_path
        )


def convert_dose_for_trial(plan, trial_info, doseInstanceUID, planInstanceUID,
                           dose_origin, patient_position, ds, export_path):
    """Convert dose for a specific trial.

    Parameters
    ----------
    plan : PinnaclePlan
        The plan object.
    trial_info : dict
        The trial dictionary containing dose and beam information.
    doseInstanceUID : str
        The DICOM UID for this dose instance.
    planInstanceUID : str
        The DICOM UID for the plan instance.
    dose_origin : list
        The dose origin [x, y, z] in mm.
    patient_position : str
        The patient position code (e.g., "HFS").
    export_path : str
        Directory where the DICOM file will be saved.
    """

    # Unpack dose origin
    dose_origin_x, dose_origin_y, dose_origin_z = dose_origin

    # Calculate trial-specific dose grid shifts
    ydoseshift = (
        trial_info["DoseGrid .VoxelSize .Y"] * 10 * trial_info["DoseGrid .Dimension .Y"]
        - trial_info["DoseGrid .VoxelSize .Y"] * 10
    )
    zdoseshift = (
        trial_info["DoseGrid .VoxelSize .Z"] * 10 * trial_info["DoseGrid .Dimension .Z"]
        - trial_info["DoseGrid .VoxelSize .Z"] * 10
    )

    # Calculate ImagePositionPatient based on patient position and trial-specific shifts
    if patient_position == "HFS":
        image_position_patient = [
            dose_origin_x,
            dose_origin_y - ydoseshift,
            dose_origin_z - zdoseshift,
        ]
    elif patient_position == "HFP":
        image_position_patient = [
            dose_origin_x,
            dose_origin_y + ydoseshift,
            dose_origin_z - zdoseshift,
        ]
    elif patient_position == "FFS":
        image_position_patient = [
            dose_origin_x,
            dose_origin_y - ydoseshift,
            dose_origin_z + zdoseshift,
        ]
    elif patient_position == "FFP":
        image_position_patient = [
            dose_origin_x,
            dose_origin_y + ydoseshift,
            dose_origin_z + zdoseshift,
        ]

    # Update trial-specific DICOM fields
    trial_name = trial_info.get("Name", "Unknown")

    # Update the SOP Instance UID for this trial
    ds.SOPInstanceUID = doseInstanceUID
    ds.file_meta.MediaStorageSOPInstanceUID = doseInstanceUID

    # Update filename with trial name and UID
    RDfilename = f"RD.{trial_name}.{doseInstanceUID}.dcm"
    ds.filename = RDfilename

    # Update trial-specific image position and grid parameters
    ds.ImagePositionPatient = image_position_patient
    ds.NumberOfFrames = int(trial_info["DoseGrid .Dimension .Z"])
    ds.Rows = int(trial_info["DoseGrid .Dimension .Y"])
    ds.Columns = int(trial_info["DoseGrid .Dimension .X"])
    ds.PixelSpacing = [
        trial_info["DoseGrid .VoxelSize .X"] * 10,
        trial_info["DoseGrid .VoxelSize .Y"] * 10,
    ]
    ds.SliceThickness = trial_info["DoseGrid .VoxelSize .Z"] * 10

    # Update trial-specific reference IDs
    ds.SeriesInstanceUID = doseInstanceUID
    ds.ReferencedRTPlanSequence[0].ReferencedSOPInstanceUID = planInstanceUID

    # Update grid frame offset vector for this trial
    grid_frame_offset_vector = []
    for p in range(0, int(trial_info["DoseGrid .Dimension .Z"])):
        grid_frame_offset_vector.append(
            p * float(trial_info["DoseGrid .VoxelSize .Z"] * 10)
        )
    ds.GridFrameOffsetVector = grid_frame_offset_vector

    # Array in which to sum the dose values of all beams
    summed_pixel_values = []

    # For each beam in the trial, convert the dose from the Pinnacle binary file and sum
    beam_list = trial_info["BeamList"] if trial_info["BeamList"] else []
    if len(beam_list) == 0:
        plan.logger.warning("No Beams found in Trial: %s. Unable to generate RTDOSE.", trial_name)
        raise MissingTrialBeamsError(f"No Beams found in Trial: {trial_name}")

    empty_beam_count = 0
    for beam in beam_list:

        plan.logger.info("Exporting Dose for beam: %s", beam["Name"])

        # Get the binary file for this beam
        binary_id = re.findall("\\d+", beam["DoseVolume"])[0]
        filled_binary_id = str(binary_id).zfill(3)
        binary_file = os.path.join(plan.path, f"plan.Trial.binary.{filled_binary_id}")

        # check whether the binary file is non-empty
        binary_data = read_binary_data(binary_file)
        if binary_data is False:
            plan.logger.warning("No Dose found for beam: %s. Skipping beam.", beam['Name'])
            empty_beam_count += 1
            if empty_beam_count == len(beam_list):
                plan.logger.error("All beams in plan are missing dose. Unable to generate RTDOSE.")
                raise MissingBeamDoseError("All beams in plan are missing dose.")
            continue

        # Get the prescription for this beam (need this for number of fractions)
        prescription = [
            p
            for p in trial_info["PrescriptionList"]
            if p["Name"] == beam["PrescriptionName"]
        ][0]

        # Get the prescription point
        plan.logger.debug("PrescriptionPointName: %s", beam["PrescriptionPointName"])
        points = plan.points
        prescription_point = []
        for p in points:
            if p["Name"] == beam["PrescriptionPointName"]:
                plan.logger.debug(
                    "Presc Point: %s %s %s %s",
                    p["Name"],
                    p["XCoord"],
                    p["YCoord"],
                    p["ZCoord"],
                )
                prescription_point = plan.convert_point(p)
                break

        if len(prescription_point) < 3:
            plan.logger.warning(
                "No valid prescription point found for beam! Beam will be ignored for "
                "Dose conversion. Dose will most likely be incorrect"
            )
            continue

        plan.logger.debug("Presc Point Dicom: %s, %s", p["Name"], prescription_point)
        total_prescription = (
            beam["MonitorUnitInfo"]["PrescriptionDose"]
            * prescription["NumberOfFractions"]
        )
        plan.logger.debug("Total Prescription %s", total_prescription)

        # Read the dose into a grid, so that we can interpolate for the prescription
        # point and determine the MU for the grid
        dose_grid = np.zeros(
            (
                trial_info["DoseGrid .Dimension .X"],
                trial_info["DoseGrid .Dimension .Y"],
                trial_info["DoseGrid .Dimension .Z"],
            )
        )
        spacing = [
            trial_info["DoseGrid .VoxelSize .X"] * 10,
            trial_info["DoseGrid .VoxelSize .Y"] * 10,
            trial_info["DoseGrid .VoxelSize .Z"] * 10,
        ]
        origin = [
            ds.ImagePositionPatient[0],
            ds.ImagePositionPatient[1],
            ds.ImagePositionPatient[2],
        ]
        dose_grid = construct_dose_from_binary(binary_data, dose_grid)

        # Get the index within that grid of the dose reference point
        idx = [0.0, 0.0, 0.0]
        orientation_matrix = np.zeros((3, 3))
        orientation_matrix[0, :] = IMAGE_ORIENTATION_MAP[patient_position][:3]
        orientation_matrix[1, :] = IMAGE_ORIENTATION_MAP[patient_position][3:]
        orientation_matrix[2, :] = np.cross(
            orientation_matrix[0, :], orientation_matrix[1, :]
        )

        for i in range(3):
            idx[i] = -(origin[i] - prescription_point[i]) / spacing[i]
            idx[i] *= orientation_matrix[i, i]

        plan.logger.debug("Index of prescription point within grid: %s", idx)

        # Trilinear interpolation of that point within the dose grid
        cgy_mu = trilinear_interpolation(idx, dose_grid)
        plan.logger.debug("cgy_mu: %s", cgy_mu)

        # Now that we have the cgy/mu value of the dose reference point, we can
        # extract an accurate value for MU
        beam_mu = (total_prescription / cgy_mu) / prescription["NumberOfFractions"]
        plan.logger.debug("Beam MU: %s", beam_mu)

        pixel_data_list = []
        for z in range(trial_info["DoseGrid .Dimension .Z"] - 1, -1, -1):
            for y in range(0, trial_info["DoseGrid .Dimension .Y"]):
                for x in range(0, trial_info["DoseGrid .Dimension .X"]):
                    value = (
                        float(prescription["NumberOfFractions"])
                        * dose_grid[x, y, z]
                        * beam_mu
                        / 100
                    )
                    pixel_data_list.append(value)

        ds.FrameIncrementPointer = ds.data_element("GridFrameOffsetVector").tag

        main_pix_array = []
        for h in range(0, trial_info["DoseGrid .Dimension .Z"]):
            pixelsforframe = []
            for k in range(
                0,
                trial_info["DoseGrid .Dimension .X"]
                * trial_info["DoseGrid .Dimension .Y"],
            ):
                pixelsforframe.append(
                    float(
                        pixel_data_list[
                            h
                            * trial_info["DoseGrid .Dimension .Y"]
                            * trial_info["DoseGrid .Dimension .X"]
                            + k
                        ]
                    )
                )

            main_pix_array = main_pix_array + list(reversed(pixelsforframe))

        main_pix_array = list(reversed(main_pix_array))

        # Add the values from this beam to the summed values
        if len(summed_pixel_values) == 0:
            summed_pixel_values = main_pix_array
        else:
            for i, values in enumerate(summed_pixel_values):
                summed_pixel_values[i] = values + main_pix_array[i]

    # Compute the scaling factor
    scale = max(summed_pixel_values) / 16384
    ds.DoseGridScaling = scale
    plan.logger.debug("Dose Grid Scaling: %s", ds.DoseGridScaling)

    # Scale by the scaling factor
    pixelvaluelist = []
    for _, element in enumerate(summed_pixel_values, 0):
        if scale != 0:
            element = round(element / scale)
        else:
            element = 0
        pixelvaluelist.append(int(element))

    # Set the PixelData
    pixel_binary_block = struct.pack("%sh" % len(pixelvaluelist), *pixelvaluelist)
    ds.PixelData = pixel_binary_block

    # If Feet first, flip the dose grid
    if patient_position in ("FFS", "FFP"):
        arr = ds.pixel_array
        ds.PixelData = np.flip(arr, axis=0).tostring()

    # Save the RTDose Dicom File
    output_file = os.path.join(export_path, RDfilename)
    plan.logger.info("Creating Dose file: %s", output_file)
    ds.save_as(output_file)
