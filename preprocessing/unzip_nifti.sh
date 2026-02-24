#!/bin/bash
read -p "Enter the subject number without leading zeros: " subject_number
#convert to 4 digits
subject_number=$(printf "%04d" "$subject_number")

data_dir="/home/data/eccolab/MNS/data/fmri/nifti/derivatives/fmriprep-25.0.0/sub-${subject_number}/func"

for gz_file in "$data_dir"/sub-*_task-*_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz; do
  [[ -e "$gz_file" ]] || continue

  nii_file="${gz_file%.gz}"

  if [[ -f "$nii_file" ]]; then
    echo "Skipping already unzipped file: $nii_file"
  else
    echo "Unzipping: $gz_file"
    gunzip -k "$gz_file"
  fi
done
