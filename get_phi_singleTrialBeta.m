%Compute phi from single trial betas
%compute for each voxel in each region and save as nifti file 
% and also average sin and cos betas across all voxels to get component average phi
[project_dir, fmri_data_dir, ~, theta_dir, ~, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();


ICA_type = 'noICA'; % 'noICA' or 'wICA'
spm1stlevel_dir = [spm1stlevel_dir, '_', ICA_type, '_wButton'];
phi_dir_voxelwise = fullfile(phi_dir, 'univariate', 'singleTrialBeta_voxelwise', ICA_type); 
phi_dir_regionwise = fullfile(phi_dir, 'univariate', 'singleTrialBeta_regionwise', ICA_type);

non_norm_theta_trials = 'includeNonNorm';%'excludeNonNorm';%'';
angle_sources = {'Subavg','Subspec'};
periodicity = {6, 4, 5, 7, 8};
runs = [1, 2]; modalities = {'word', 'face'};
brain_atlas = load_atlas('canlab2018');
smooth = 'unsmoothed'; phi_dir_regionwise = fullfile(phi_dir_regionwise, smooth); phi_dir_voxelwise = fullfile(phi_dir_voxelwise, smooth); spm1stlevel_dir = fullfile(spm1stlevel_dir, smooth);


region_names = {'OFC2016ConstantinescuR5','vmPFC2019BaoR5','vmPFC2016ConstantinescuR5','vmPFC2019Bao2R5','HC','ERC'};
region_masks = {fullfile(project_dir, 'masks', 'OFC_2016Constantinescu_r5.nii'),...
                fullfile(project_dir, 'masks', 'vmPFC_2019Bao_r5.nii'),...
                fullfile(project_dir, 'masks', 'vmPFC_2016Constantinescu_r5.nii'),...
                fullfile(project_dir, 'masks', 'vmPFC_2019Bao2_r5.nii'),...
                fullfile(project_dir, 'masks', 'vmPFC_2016Constantinescu_r5.nii'),...
                fullfile(project_dir, 'masks', 'HC_Julich.nii'),...
                fullfile(project_dir, 'masks', 'ERC_Julich.nii')};


modality_to_use = 'both_modalities';
cross_validations = {'xModalityRun', 'xModality', 'xRun'};%{'LORO'};%{'xModalityRun', 'xModality', 'xRun', 'avgAllRuns'};
phi_averaging_methods = 'component_average';
save_voxelwise_phi = true;
save_regionwise_phi = true;


[brain_data_all, brain_modality_run_idx, brain_sub_ids] = load_and_prepare_brain_data(subjects, spm1stlevel_dir, region_masks, region_names, modality_to_use); %a struct with concatenated brain data (across runs and subjects) for each region
template_nifti_object = fmri_data(fullfile(spm1stlevel_dir, ['sub', subjects{1}], 'singleTrial', 'beta_0001.nii'));

for p = 1:length(periodicity)
    current_periodicity = periodicity{p};
    fprintf('Processing periodicity %d...\n', current_periodicity);

    
    for t = 1:length(angle_sources)
        current_angle_source = angle_sources{t};

        [angle_data_all, subject_ids, modality_run_ids] = load_and_prepare_angle_data(subjects, runs, spm1stlevel_dir, theta_dir, current_angle_source, modality_to_use);
        %make sure brain_modality_run_idx and angles_modality_run_idx are the same
        if ~isequal(brain_modality_run_idx, modality_run_ids)
            disp(brain_modality_run_idx); disp(modality_run_ids); error('brain_modality_run_idx and modality_run_ids are not the same');
        end
        
        for s = 1:length(subjects)
            current_subject = subjects{s};

            
            for region_idx = 1:length(region_names)
                phiRadDivByPeriod_struct = struct();
                current_region_name = region_names{region_idx}; current_region_mask = region_masks{region_idx};
                brain_data_region = brain_data_all.(current_region_name);
                tmp_obj_region = apply_mask(template_nifti_object, current_region_mask);

                for cv = 1:length(cross_validations)
                    current_cross_validation = cross_validations{cv};
        
                    if strcmp(current_cross_validation, 'xModalityRun')
                        training_idx = {{'face1', 'word1'}, {'face1', 'word2'}, {'face2', 'word1'}, {'face2', 'word2'}};
                        test_idx = {{'face2', 'word2'}, {'face2', 'word1'}, {'face1', 'word2'}, {'face1', 'word1'}};
                    elseif strcmp(current_cross_validation, 'xModality')
                        training_idx = {{'face1','face2'}, {'word1', 'word2'}};
                        test_idx = {{'word1', 'word2'}, {'face1', 'face2'}};
                    elseif strcmp(current_cross_validation, 'xRun')
                        training_idx = {{'face1'}, {'face2'}, {'word1'}, {'word2'}};
                        test_idx = {{'face2'}, {'face1'}, {'word2'}, {'word1'}};
                    elseif strcmp(current_cross_validation, 'avgAllRuns')
                        training_idx = {{'face1','face2', 'word1', 'word2'}};
                        test_idx = {{'face1','face2', 'word1', 'word2'}};
                    elseif strcmp(current_cross_validation, 'LORO')
                        training_idx = {{'face2', 'word1','word2'}, {'face1', 'word1','word2'}, {'face1', 'face2','word2'}, {'face1', 'face2','word1'}};
                        test_idx = {{'face1'}, {'face2'}, {'word1'}, {'word2'}};
                    end

                
                    if strcmp(phi_averaging_methods, 'component_average')
                        beta_sin_allvoxels = nan(size(brain_data_region, 2), 1);
                        beta_cos_allvoxels = nan(size(brain_data_region, 2), 1);
                        
                        for t_idx = 1:length(training_idx)
                            current_training_idx = training_idx{t_idx}; current_test_idx = test_idx{t_idx};
                            mask = strcmp(subject_ids, current_subject) & ismember(modality_run_ids, current_training_idx);
                            X_train = angle_data_all(mask, :);
                            Y_train = brain_data_region(mask, :);
                            
                            X_train_noPhi = [cos(current_periodicity * deg2rad(X_train)), sin(current_periodicity * deg2rad(X_train))];
                            X_design_noPhi = [ones(size(X_train_noPhi, 1), 1), X_train_noPhi];

                            for v = 1:size(Y_train, 2)
                                y = Y_train(:, v);
                                b_noPhi = X_design_noPhi \ y;
                                beta_cos_allvoxels(v) = b_noPhi(2);
                                beta_sin_allvoxels(v) = b_noPhi(3);
                            end
                            if save_voxelwise_phi
                                tmp_obj_region.dat = atan2(beta_sin_allvoxels, beta_cos_allvoxels) / current_periodicity;
                                path_to_save = fullfile(phi_dir_voxelwise, ['sub', current_subject], non_norm_theta_trials, ['periodicity', num2str(current_periodicity)], ['theta', current_angle_source],'nifti', current_cross_validation);
                                if ~exist(path_to_save, 'dir'), mkdir(path_to_save); end
                                tmp_obj_region.fullpath = fullfile(path_to_save, ['voxelwisePhiRadDivByPeriod_',strjoin(current_training_idx, '_'), current_region_name, '.nii']);
                                tmp_obj_region.write;

                                %save to csv file
                                csv_path = fullfile(phi_dir_voxelwise, ['sub', current_subject], non_norm_theta_trials, ['periodicity', num2str(current_periodicity)], ['theta', current_angle_source],'csv', current_cross_validation);
                                if ~exist(csv_path, 'dir'), mkdir(csv_path); end
                                csv_file = fullfile(csv_path, ['voxelwisePhiRadDivByPeriod_',strjoin(current_training_idx, '_'), current_region_name, '.csv']);
                                writetable(array2table(tmp_obj_region.dat), csv_file);
                            end
                            if save_regionwise_phi
                                phiRadDivByPeriod_struct.crossValidations.(current_cross_validation).(strjoin(current_training_idx, '_')) = atan2(mean(beta_sin_allvoxels), mean(beta_cos_allvoxels)) / current_periodicity;
                            end
                        end
                    end
                end
                if save_regionwise_phi
                    path_to_save = fullfile(phi_dir_regionwise, ['sub', current_subject], non_norm_theta_trials, ['periodicity', num2str(current_periodicity)], ['theta', current_angle_source]);
                    if ~exist(path_to_save, 'dir'), mkdir(path_to_save); end
                    phi_file_name = ['PhiRadDivByPeriod_', region_names{region_idx}, '.mat'];
                    save_and_merge_phi_data(path_to_save, phi_file_name, phiRadDivByPeriod_struct);
                end
            end 
        end
    end
end


function [angle_data_all, subject_ids, modality_run_ids] = load_and_prepare_angle_data(subjects, runs, spm1stlevel_dir, theta_dir, current_angle_source, modality_to_use)

    angle_data_all = [];
    subject_ids = {};
    modality_run_ids = {};

    for s = 1:length(subjects)
        current_subject = subjects{s};
        
        singleTrial_dir = fullfile(spm1stlevel_dir, ['sub', current_subject], 'singleTrial');
        load(fullfile(singleTrial_dir, 'SPM.mat'));
        spm_regressor_names = SPM.xX.name;
        stim_trials = contains(spm_regressor_names, 'trial');
        stim_trials_idx = find(stim_trials);
        
        first_stim_trial = spm_regressor_names(find(stim_trials, 1, 'first'));
        last_stim_trial = spm_regressor_names(find(stim_trials, 1, 'last'));
        
        if strcmp(modality_to_use, 'both_modalities') && contains(first_stim_trial, 'face') && contains(last_stim_trial, 'word')
            modalities = {'face', 'word'};
        elseif strcmp(modality_to_use, 'both_modalities') && contains(first_stim_trial, 'word') && contains(last_stim_trial, 'face')
            modalities = {'word', 'face'};
        elseif strcmp(modality_to_use, 'face') || strcmp(modality_to_use, 'word')
            modalities = {modality_to_use};
        else
            error('Modality order unclear.');
        end
        
        % Load angle data
        angles_all_modalities_runs = [];
        angles_modality_run_idx = {};
        for m = 1:length(modalities)
            current_modality = modalities{m};
            for r = 1:length(runs)
                current_run = runs(r);
                angle_file = fullfile(theta_dir, ['sub', current_subject], ['sub', current_subject, '_', current_modality, '_thetas_run', num2str(current_run), '.mat']);
                load(angle_file);
                
                angle_data = thetas.(current_angle_source)';
                angle_avg_data = thetas.Subavg';
                if strcmp(current_angle_source, 'Subspec')
                    nan_idx = isnan(angle_data);
                    angle_data(nan_idx) = angle_avg_data(nan_idx);
                end
                
                angles_all_modalities_runs = [angles_all_modalities_runs; angle_data];
                current_modality_run_idx = [current_modality, num2str(current_run)];
                angles_modality_run_idx = [angles_modality_run_idx; repmat({current_modality_run_idx}, size(angle_data, 1), 1)];
            end
        end
        angle_data_all = [angle_data_all; angles_all_modalities_runs];
        subject_ids = [subject_ids; repmat({current_subject}, size(angles_all_modalities_runs, 1), 1)];
        modality_run_ids = [modality_run_ids; angles_modality_run_idx];
    end
end

function [brain_data_all, brain_modality_run_idx, brain_sub_ids] = load_and_prepare_brain_data(subjects, spm1stlevel_dir, region_masks, region_names, modality_to_use)
    brain_data_all = struct();
    beta_images = {};
    brain_modality_run_idx = {};
    brain_sub_ids = {};

    for s = 1:length(subjects)
        current_subject = subjects{s};
        
        singleTrial_dir = fullfile(spm1stlevel_dir, ['sub', current_subject], 'singleTrial');
        load(fullfile(singleTrial_dir, 'SPM.mat'));
        spm_regressor_names = SPM.xX.name;
        stim_trials = contains(spm_regressor_names, 'trial');
        stim_trials_idx = find(stim_trials);
        
        % Prepare beta images
        beta_images_current_subject = {};
        brain_modality_run_idx_current_subject = {};
        brain_sub_ids_current_subject = {};
        for i = 1:length(stim_trials_idx)
            current_trial = stim_trials_idx(i);
            current_beta_image = fullfile(singleTrial_dir, ['beta_', num2str(current_trial, '%04d'), '.nii']);
            beta_images_current_subject{i, 1} = current_beta_image;
            stim_name = spm_regressor_names{current_trial};
            stim_info = regexp(stim_name, 'Sn\((\d+)\)\s+(\w+)_run(\d+)_trial(\d+)\*bf\((\d+)\)', 'tokens');
            brain_modality_run_idx_current_subject{i, 1} = [stim_info{1}{2}, stim_info{1}{3}];
            brain_sub_ids_current_subject{i, 1} = current_subject;
        end
        beta_images = [beta_images; beta_images_current_subject];
        brain_modality_run_idx = [brain_modality_run_idx; brain_modality_run_idx_current_subject];
        brain_sub_ids = [brain_sub_ids; brain_sub_ids_current_subject];
    end
    if strcmp(modality_to_use, 'face') || strcmp(modality_to_use, 'word')
        modality_subset = contains(brain_modality_run_idx, modality_to_use);
        brain_modality_run_idx = brain_modality_run_idx(modality_subset);
        brain_sub_ids = brain_sub_ids(modality_subset);
        beta_images = beta_images(modality_subset);
        fprintf('Using %s modality for subject %s...\n', modality_to_use, current_subject);
    end
    % Load brain data
    brain_data = fmri_data(beta_images);
    for r = 1:length(region_masks)
        current_region_mask = region_masks{r}; current_region_name = region_names{r};
        brain_data_region_obj = apply_mask(brain_data, current_region_mask);
        brain_data_all.(current_region_name) = brain_data_region_obj.dat';
    end
end

function save_and_merge_phi_data(phi_file_dir, phi_file_name, new_data)
    full_file_path = fullfile(phi_file_dir, phi_file_name);
    
    if exist(full_file_path, 'file')
        % load and merge with existing data
        existing_data = load(full_file_path); existing_data = existing_data.phiRadDivByPeriod_struct;
        fields = fieldnames(new_data.crossValidations);
        for i = 1:length(fields)
            existing_data.crossValidations.(fields{i}) = new_data.crossValidations.(fields{i}); %this will overwrite the existing data if the same fieldname exists in both new and existing data
        end
        phiRadDivByPeriod_struct = existing_data;
        save(full_file_path, 'phiRadDivByPeriod_struct');
    else
        phiRadDivByPeriod_struct = new_data;
        save(full_file_path, 'phiRadDivByPeriod_struct');
    end
end
