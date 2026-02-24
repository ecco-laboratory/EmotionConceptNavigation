#!/bin/bash
#SBATCH --account=default
#SBATCH --exclude=node3
#SBATCH --mem=48G
#SBATCH --partition week-long  # Queue names you can submit to
# Outputs ----------------------------------
#SBATCH -o ./slurm/R-%x.%A_%a.out
#SBATCH -e ./slurm/R-%x.%A_%a.err
#SBATCH --mail-user=yma355@emory.edu
#SBATCH --mail-type=ALL
# ------------------------------------------

read -p "Enter subj num to warp with no leading 0s (but with spaces between each subj): " SUBJ_NUM
SUBJ_NUMS_4D=()
for num in $SUBJ_NUM; do
    SUBJ_NUMS_4D+=("sub-$(printf '%04d' "$num")")
done

# First, let's define some paths to the fMRIPrep and tedana outputs
fmriprep_dir="/home/data/eccolab/MNS/data/fmri/derivatives/fmriprep-25.0.0"
tedana_dir="/home/data/eccolab/MNS/data/fmri/derivatives/fmriprep-25.0.0/tedana"

standard_space="MNI152NLin2009cAsym"

# find all denoised bold files in the tedana output
find "${tedana_dir}" -maxdepth 1 -type d -name "sub-*" | while IFS= read -r subject_dir; do
  subject=$(basename "${subject_dir}")

  # check if this subject is in the input list
  if [[ ! " ${SUBJ_NUMS_4D[*]} " =~ " ${subject} " ]]; then
    continue
  fi

  echo "Processing ${subject}..."

  # find all denoised bold files for the current subject
  find "${subject_dir}" -type f -name "sub-*_desc-denoised_bold.nii.gz" | while IFS= read -r file_to_warp; do
    echo "  Found denoised file: ${file_to_warp}"

  # get task/run identifier
  identifier=$(basename "${file_to_warp}" | sed -e "s/${subject}_//g" -e "s/_desc-denoised_bold.nii.gz//g")
  echo "  Identifier: ${identifier}"

  # construct the output filename
  out_file="${fmriprep_dir}/func/${subject}/${subject}_${identifier}_space-${standard_space}_desc-optcomDenoised_bold.nii.gz"
  echo "  Output file: ${out_file}"

  # construct the path to the corresponding standard space bold reference image from fMRIPrep
  standard_space_file="${fmriprep_dir}/func/${subject}/${subject}_${identifier}_space-${standard_space}_boldref.nii.gz"

  # construct the paths to the necessary transforms from fMRIPrep
  xform_native_to_t1w="${fmriprep_dir}/func/${subject}/${subject}_${identifier}_from-scanner_to-T1w_mode-image_xfm.txt"
  xform_t1w_to_std="${fmriprep_dir}/anat/${subject}/${subject}_from-T1w_to-${standard_space}_mode-image_xfm.h5"

  # check if the necessary files exist before running ANTs
  if [ -f "${standard_space_file}" ] && [ -f "${xform_native_to_t1w}" ] && [ -f "${xform_t1w_to_std}" ]; then
    echo "  applying transforms..."
    antsApplyTransforms \
      -e 3 \
      -i "${file_to_warp}" \
      -r "${standard_space_file}" \
      -o "${out_file}" \
      -n LanczosWindowedSinc \
      -t "${xform_t1w_to_std}" \
      -t "${xform_native_to_t1w}" # apply native to T1w first
    echo "  Registration complete for ${file_to_warp}"
  else
    echo "  Warning: Missing necessary files for ${file_to_warp}. Skipping registration."
    if [ ! -f "${standard_space_file}" ]; then
      echo "    Missing: ${standard_space_file}"
    fi
    if [ ! -f "${xform_native_to_t1w}" ]; then
      echo "    Missing: ${xform_native_to_t1w}"
    fi
    if [ ! -f "${xform_t1w_to_std}" ]; then
      echo "    Missing: ${xform_t1w_to_std}"
    fi
  fi
  echo "" 
done

echo "Finished processing all denoised files."