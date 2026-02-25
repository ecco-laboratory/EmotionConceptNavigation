# A grid-like basis for affective space in human ventromedial prefrontal cortex

In this fMRI study, we tested whether abstract emotion concepts are represented in the human brain using grid-like codes.

## Dependencies 
This code uses [Canlab Core Tools](https://github.com/canlab/CanlabCore/tree/master), which is an object oriented toolbox that uses the [SPM software](https://www.fil.ion.ucl.ac.uk/spm/) to process fMRI data. Instructions for installing CANlab Core Tools which provided many of the functions and Neuroimaging Pattern Masks used in the analyses can be found [here](https://canlab.github.io/_pages/canlab_help_1_installing_tools/canlab_help_1_installing_tools.html).

Instructions for installing SPM12 can be found [here](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/).

Analyses were performed using MATLAB R2024a, Python 3.9.18, and R 4.4.1 on Ubuntu 20.04.6 LTS (Focal Fossa). No non-standard or specialized hardware was required.

## Summary of analyses

**Behavioral analyses**
- Get trajectory angles from ratings 
  (`get_theta_table.ipynb`) 
   - Output: `./data/beh/behTables/subAvg_table_face.csv` & `./data/beh/behTables/subAvg_table_word.csv`
- Visualize behavioral data 
  (`plot_beh.ipynb`). 
- Test relationship between distance and choice behavior
  (`beh_consistency_stats.R`). 

**Brain analyses**
1. GLM
    - Get timing of stimuli  
  (`get_spm_condition_timing_faceRatingButton.m`  
  `get_spm_condition_timing_wordRatingButton.m`)
    - Get motion regressors from fMRIPrep output  
    (`get_motion_regressors.m`)
    - Fit GLM to obtain single-trial betas for each subject  
    (`estimate_spm_level1model_singleTrial.m`)

2. Hexadirectional modulation analysis
    - Prepare trajectory angle file  
    (`get_theta_vaDiff_distance_category.m`)
    - Estimate grid orientation (phi)  
    (`get_phi_singleTrialBeta.m`)
    - Compute aligned–misaligned pattern similarity differences  
    (`compare_multivariateOnOffGrid_singleTrialBeta_wPhi.m`)
        - Output:  
    `./outputs/singleTrialBetaAnalysis/noICA/incl_all_subs_trials/onoffGridcontrast_multivariate/unsmoothed/Subavg/periodicity6/xModalityRun/phi_voxelComponentAverage/singleTrialBeta/csv/contrast_values.csv`
    - Visualize results  
    (`plot_brain_results.ipynb`)

3. Distance encoding & choice behavior
    - Fit encoding model and compute pattern expression score  
  (`encode_distance_singleTrialBeta_loso_multivariate.m`)
    - Save results to CSV  
    (`save_beh_Yproj_to_csv.m`)
        - Output:
    `./outputs/singleTrialBetaAnalysis/noICA/incl_all_subs_trials/distanceCentered_multivariate_encoding_yhat_ytest/unsmoothed/both_modalities/Subavg/va/csv/beh_Yproj_raw_wDistance_OFC2016ConstantinescuR5.csv`
    - Test relationship between distance expression and choice behavior  
    (`beh_consistency_stats_brain.R`)
    - Visualize results  
    (`plot_brain_results.ipynb`)

