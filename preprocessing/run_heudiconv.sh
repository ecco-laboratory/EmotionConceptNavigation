#!/bin/bash
# This study uses a single heuristic file that should be consistent across all subjects
# because they all get scanned with the exact same sequence
# This also only runs in an interactive bash instance, not in slurm
# Monica tried setting it up to run through slurm and it appeared to convert the data okay
# but then didn't update the participants.tsv file with the finished subjects??

#singularity run -B /archival/projects/MNS/data/fmri:/base \
# /home/data/shared/SingularityImages/heudiconv_0.5.4.sif \
# -d /base/dicom/mns_{subject}*/*/*.dcm \
# -o /base/nifti/dicominfo/ \
# -f convertall -s 9999 -c none 

read -p "Enter subj num to heuristic with no leading 0s: " SUBJ_NUM

printf -v SUBJ_NUM_4D '%04d' $SUBJ_NUM 
if (( SUBJ_NUM % 2 == 0 )); then
    HEURISTIC="/base/heuristic_evensubj.py"
else
    HEURISTIC="/base/heuristic_oddsubj.py"
fi

singularity run -B /archival/projects/MNS/data/fmri:/base \
/home/data/shared/SingularityImages/heudiconv_1.3.3.sif \
-d /base/dicom/mns_{subject}*/*/*.dcm \
-o /base/nifti \
-f "$HEURISTIC" \
-s $SUBJ_NUM_4D \
-c dcm2niix \
-b \
--overwrite

