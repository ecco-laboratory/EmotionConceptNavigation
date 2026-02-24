[project_dir, fmri_data_dir, ~, theta_dir, category_dir, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();

ICA_type = 'noICA'; % 'noICA' or 'wICA'
spm1stlevel_dir = [spm1stlevel_dir, '_', ICA_type, '_wButton'];
beh_sources = {'Subavg','Subspec'};
smooth = 'unsmoothed'; spm1stlevel_dir = fullfile(spm1stlevel_dir, smooth);

runs = [1, 2]; 
modality_to_use = {'both_modalities', 'face', 'word'};
brain_atlas = load_atlas('canlab2018');
region_group = 'hcecvmpfc';
region_names = {'OFC2016ConstantinescuR5','HC', 'ERC'};
region_masks = {fullfile(project_dir, 'masks', 'OFC_2016Constantinescu_r5.nii'),...
                fullfile(project_dir, 'masks', 'HC_Julich.nii'),...
                fullfile(project_dir, 'masks', 'ERC_Julich.nii')};

save_nifti_files_yhat_ytest = false;       
save_nifti_files_performance = false;    
save_original_results = false;
save_Yproj = true;
if_process_noise_ceiling = false;
run_permutation = false;   
n_permutations = 3;     
shuffle_type = 'same_shuffle_across_conditions';

va_dir = fullfile(project_dir, 'data', 'beh', 'VA');
inscan_judgment_dir = fullfile(project_dir, 'data', 'beh', 'inscan_judgment');
distance_types = {'catDominantAvg','catDominantEuclidean','catMaxEuclidean','vSigned','aSigned','judgedSigned','va', 'judged', 'catCosine', 'catMax', 'v', 'a'};

job_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
[distance_type_idx, modality_to_use_idx] = ind2sub([numel(distance_types), numel(modality_to_use)], job_id);
distance_types = distance_types(distance_type_idx);
modality_to_use = modality_to_use{modality_to_use_idx};

yhat_output_dir = fullfile(project_dir, 'outputs', 'singleTrialBetaAnalysis', ICA_type, 'incl_all_subs_trials', 'distanceCentered_multivariate_encoding_yhat_ytest', smooth, modality_to_use); if ~exist(yhat_output_dir, 'dir'), mkdir(yhat_output_dir); end  
performance_output_dir = fullfile(project_dir, 'outputs', 'singleTrialBetaAnalysis', ICA_type, 'incl_all_subs_trials', 'distanceCentered_multivariate_encoding_performance', smooth, modality_to_use); if ~exist(performance_output_dir, 'dir'), mkdir(performance_output_dir); end
voxelperformance_output_dir = fullfile(project_dir, 'outputs', 'singleTrialBetaAnalysis', ICA_type, 'incl_all_subs_trials', 'distanceCentered_multivariate_encoding_performance_voxelwise', smooth, modality_to_use); if ~exist(voxelperformance_output_dir, 'dir'), mkdir(voxelperformance_output_dir); end


[brain_data_all, brain_modality_run_idx, brain_sub_ids] = load_and_prepare_brain_data(subjects, spm1stlevel_dir, region_masks, region_names, modality_to_use); %a struct with concatenated brain data (across runs and subjects) for each region
template_nifti_object = fmri_data(fullfile(spm1stlevel_dir, ['sub', subjects{1}], 'singleTrial', 'beta_0001.nii'));

category_formats = {'startContinuous', 'endContinuous'};
[category_data_concatSubs_allFormatsSources, category_info] = load_and_prepare_category_data(subjects, runs, spm1stlevel_dir, category_dir, beh_sources, category_formats, modality_to_use); 
[va_data_allSources, subject_ids, modality_run_ids] = load_and_prepare_va_data(subjects, runs, spm1stlevel_dir, va_dir, beh_sources, modality_to_use); 

if ~isequal(brain_modality_run_idx, modality_run_ids), disp(brain_modality_run_idx); disp(modality_run_ids); error('brain_modality_run_idx and modality_run_ids are not the same'); end
if ~isequal(brain_sub_ids, subject_ids), disp(brain_sub_ids); disp(subject_ids); error('brain_sub_ids and subject_ids are not the same'); end
if ~isequal(brain_modality_run_idx, category_info.(beh_sources{1}).(category_formats{1}).modality_run), disp(brain_modality_run_idx); disp(category_info.(beh_sources{1}).(category_formats{1}).modality_run);
    error('brain_modality_run_idx and category_info.(beh_sources{1}).(category_formats{1}).modality_run are not the same');
end
if ~isequal(brain_sub_ids, category_info.(beh_sources{1}).(category_formats{1}).sub_ids), disp(brain_sub_ids); disp(category_info.(beh_sources{1}).(category_formats{1}).sub_ids);
    error('brain_sub_ids and category_info.(beh_sources{1}).(category_formats{1}).sub_ids are not the same');
end
clear category_info;

if run_permutation
    if strcmp(shuffle_type, 'same_shuffle_across_conditions')
        permutation_indices = generate_permutation_indices_loso(subject_ids, modality_run_ids, n_permutations, subjects);
    end
end

%% LOSO encoding with permutation testing
for t = 1:length(beh_sources)
    current_beh_source = beh_sources{t};
    category_data_all = category_data_concatSubs_allFormatsSources.(current_beh_source); 
    va_data_all = va_data_allSources.(current_beh_source); 
    [~, comparison_data_all,~,~] = load_and_prepare_scanbeh_data(subjects, runs, spm1stlevel_dir, inscan_judgment_dir, current_beh_source, modality_to_use); 
    for c = 1:length(distance_types)
        current_distance_type = distance_types{c};
        
        
        current_obs_output_dir = fullfile(performance_output_dir, 'csv', current_beh_source, current_distance_type);
        current_perm_output_dir = fullfile(performance_output_dir, shuffle_type, current_beh_source, current_distance_type);
        current_yhat_output_dir = fullfile(yhat_output_dir, current_beh_source, current_distance_type);
        current_voxelperformance_output_dir = fullfile(voxelperformance_output_dir, current_beh_source, current_distance_type);
        if ~exist(current_obs_output_dir, 'dir'), mkdir(current_obs_output_dir); end
        if ~exist(current_perm_output_dir, 'dir'), mkdir(current_perm_output_dir); end
        if ~exist(current_yhat_output_dir, 'dir'), mkdir(current_yhat_output_dir); end
        if ~exist(current_voxelperformance_output_dir, 'dir'), mkdir(current_voxelperformance_output_dir); end

        X_original = prepare_X_original(category_data_all, va_data_all, comparison_data_all, current_distance_type);
        
        if run_permutation || save_original_results || save_Yproj || save_nifti_files_yhat_ytest || save_nifti_files_performance
            tic;
            region_results = perform_loso_encoding(brain_data_all, X_original, ...
                                                subject_ids, subjects, region_masks, region_names, ...
                                                    save_nifti_files_yhat_ytest, save_nifti_files_performance, save_Yproj, ...
                                                current_yhat_output_dir, current_voxelperformance_output_dir,...
                                                    template_nifti_object, project_dir);
            csv_file_performance = fullfile(current_obs_output_dir, ['loso_performance_', region_group, '.csv']);
            if save_original_results
                writetable(region_results, csv_file_performance);
            end
            display('Time taken to perform LOSO encoding:');
            toc;
        end
        if if_process_noise_ceiling
            region_results = perform_loso_encoding_noise_ceiling(brain_data_all, X_original, subject_ids, subjects, region_names);
            csv_file_performance = fullfile(current_obs_output_dir, ['loso_performance_noiseceilingTrainAvgTestCorr_', region_group, '.csv']);
            writetable(region_results, csv_file_performance);
        end
        if run_permutation
            tic;
            fprintf('  Running permutation test...\n');
            
            [p_values, null_distribution, permutation_fold_results] = run_permutation_test(brain_data_all, X_original, region_results, ...
                                          subject_ids, modality_run_ids, subjects, ...
                                          region_masks, region_names, ...
                                          permutation_indices, shuffle_type, project_dir);
            display('Time taken to run permutation test:');
            toc;
            %save_performance_results_wperm(region_results, p_values, null_distribution, csv_file_performance_wperm);
            save_permutation_details(region_results, permutation_fold_results, current_perm_output_dir);
        end
       
    end
end

fprintf('LOSO encoding analysis completed!\n');

function [X_original] = prepare_X_original(category_data_all, va_data_all, comparison_data_all, current_distance_type)
    vSignedDist = va_data_all.endValence - va_data_all.startValence;
    aSignedDist = va_data_all.endArousal - va_data_all.startArousal;
    vDist = abs(vSignedDist);
    aDist = abs(aSignedDist);
    vaDist = sqrt(vDist.^2 + aDist.^2);
    
    if strcmp(current_distance_type, 'va')
        X_original = vaDist - mean(vaDist);
    elseif strcmp(current_distance_type, 'v')
        X_original = vDist - mean(vDist);
    elseif strcmp(current_distance_type, 'a')
        X_original = aDist - mean(aDist);
    elseif strcmp(current_distance_type, 'vSigned')
        X_original = vSignedDist - mean(vSignedDist);
    elseif strcmp(current_distance_type, 'aSigned')
        X_original = aSignedDist - mean(aSignedDist);
    elseif strcmp(current_distance_type, 'judged') || strcmp(current_distance_type, 'judgedSigned')
        judgedDist = zeros(size(comparison_data_all, 1), 1); signed_judgedDist = zeros(size(comparison_data_all, 1), 1);
        for t = 1:size(comparison_data_all, 1)
            if strcmpi(comparison_data_all{t}, 'PLEASANTNESS')
                judgedDist(t) = vDist(t); signed_judgedDist(t) = vSignedDist(t);
            elseif strcmpi(comparison_data_all{t}, 'ACTIVATION')
                judgedDist(t) = aDist(t); signed_judgedDist(t) = aSignedDist(t);
            else
                error('Unexpected comparison_data value at trial %d', t);
            end
        end
        if strcmp(current_distance_type, 'judged'), X_original = judgedDist - mean(judgedDist); elseif strcmp(current_distance_type, 'judgedSigned'), X_original = signed_judgedDist - mean(signed_judgedDist); end
    elseif strcmp(current_distance_type, 'catCosine') || strcmp(current_distance_type, 'catMax') || strcmp(current_distance_type, 'catMaxEuclidean') || strcmp(current_distance_type, 'catDominantAvg') || strcmp(current_distance_type, 'catDominantEuclidean')
        startMat = table2array(category_data_all.('startContinuous')); endMat = table2array(category_data_all.('endContinuous')); 
        nTrials  = size(startMat,1);
        cosineDist = zeros(nTrials,1); catMaxDist = zeros(nTrials,1); catDominantEuclideanDist = zeros(nTrials,1); catDominantAvgDist = zeros(nTrials,1); catMaxEuclideanDist = zeros(nTrials,1);
        [~, startMaxIdx] = max(startMat, [], 2); [~, endMaxIdx]   = max(endMat, [], 2);
        for t = 1:nTrials
            s = startMat(t,:); e = endMat(t,:);
            cosineDist(t) = 1 - (dot(s,e) / (norm(s)*norm(e) + eps));
            catMaxDist(t) = mean([abs(s(startMaxIdx(t)) - e(startMaxIdx(t))), abs(s(endMaxIdx(t))   - e(endMaxIdx(t)))]);
            if startMaxIdx(t) == endMaxIdx(t)
                catMaxEuclideanDist(t) = catMaxDist(t);
                catDominantEuclideanDist(t) = 0; 
                catDominantAvgDist(t) = 0;
            else
                catMaxEuclideanDist(t) = sqrt((s(startMaxIdx(t)) - e(startMaxIdx(t)))^2 + (s(endMaxIdx(t)) - e(endMaxIdx(t)))^2);
                catDominantEuclideanDist(t) = sqrt((s(startMaxIdx(t)))^2 + (e(endMaxIdx(t)))^2);
                catDominantAvgDist(t) = mean([s(startMaxIdx(t)),e(endMaxIdx(t))]);
            end
        end
        if strcmp(current_distance_type, 'catCosine')
            X_original = cosineDist - mean(cosineDist);
        elseif strcmp(current_distance_type, 'catMax')
            X_original = catMaxDist - mean(catMaxDist);
        elseif strcmp(current_distance_type, 'catMaxEuclidean')
            X_original = catMaxEuclideanDist - mean(catMaxEuclideanDist);
        elseif strcmp(current_distance_type, 'catDominantEuclidean')
            X_original = catDominantEuclideanDist - mean(catDominantEuclideanDist);
        elseif strcmp(current_distance_type, 'catDominantAvg')
            X_original = catDominantAvgDist - mean(catDominantAvgDist);
        end
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

function modalities = get_modality_order(spm1stlevel_dir, subject, modality_to_use)
    singleTrial_dir = fullfile(spm1stlevel_dir, ['sub', subject], 'singleTrial');
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
end

function [consistency_data_all, comparison_data_all, subject_ids, modality_run_ids] = load_and_prepare_scanbeh_data(subjects, runs, spm1stlevel_dir, scanbeh_dir, current_beh_source, modality_to_use)

    consistency_data_all = [];
    comparison_data_all = [];
    subject_ids = {};
    modality_run_ids = {};

    for s = 1:length(subjects)
        current_subject = subjects{s};
        modalities = get_modality_order(spm1stlevel_dir, current_subject, modality_to_use);
        
        % Load angle data
        consistency_all_modalities_runs = []; comparison_all_modalities_runs = [];
        scanbeh_modality_run_idx = {}; 
        for m = 1:length(modalities)
            current_modality = modalities{m};
            for r = 1:length(runs)
                current_run = runs(r);
                scanbeh_file = fullfile(scanbeh_dir, ['sub', current_subject], ['sub', current_subject, '_', current_modality, '_consistency_run', num2str(current_run), '.mat']);
                comparison_file = fullfile(scanbeh_dir, ['sub', current_subject], ['sub', current_subject, '_', current_modality, '_comparison_run', num2str(current_run), '.mat']);
                load(scanbeh_file); load(comparison_file);
                
                consistency_data = inscan_consistency_struct.(current_beh_source)'; comparison_data = inscan_comparison_struct.(current_beh_source)';
                consistency_avg_data = inscan_consistency_struct.Subavg';
                if strcmp(current_beh_source, 'Subspec')
                    nan_idx = isnan(consistency_data);
                    consistency_data(nan_idx) = consistency_avg_data(nan_idx);
                end
                
                consistency_all_modalities_runs = [consistency_all_modalities_runs; consistency_data]; comparison_all_modalities_runs = [comparison_all_modalities_runs; comparison_data];
                current_modality_run_idx = [current_modality, num2str(current_run)];
                scanbeh_modality_run_idx = [scanbeh_modality_run_idx; repmat({current_modality_run_idx}, size(consistency_data, 1), 1)]; 
            end
        end
        consistency_data_all = [consistency_data_all; consistency_all_modalities_runs]; comparison_data_all = [comparison_data_all; comparison_all_modalities_runs];
        subject_ids = [subject_ids; repmat({current_subject}, size(consistency_all_modalities_runs, 1), 1)];
        modality_run_ids = [modality_run_ids; scanbeh_modality_run_idx];
    end
end


function [category_data_concatSubs_allFormatsSources, category_info] = load_and_prepare_category_data(subjects, runs, spm1stlevel_dir, category_dir, category_sources, category_formats, modality_to_use)

    category_data_concatSubs_allFormatsSources = struct();
    category_info = struct();

    for s = 1:length(subjects)
        current_subject = subjects{s};
        modalities = get_modality_order(spm1stlevel_dir, current_subject, modality_to_use);
        
        % Load category data
        for m = 1:length(modalities)
            current_modality = modalities{m};
            for r = 1:length(runs)
                current_run = runs(r);
                category_file = fullfile(category_dir, ['sub', current_subject], ['sub', current_subject, '_', current_modality, '_categories_run', num2str(current_run), '.mat']);
                load(category_file);

                
                for cat_src = 1:length(category_sources)
                    cat_source = category_sources{cat_src};
                    for cat_fmt = 1:length(category_formats)
                        cat_format = category_formats{cat_fmt};
    
                        Y = categories.(cat_format).(cat_source);
                        subavg = categories.(cat_format).Subavg;
                        % fill NaNs in subspec with subavg
                        if strcmp(cat_source, 'Subspec')
                            for col = 1:width(Y)
                                nan_idx = isnan(Y{:,col});
                                Y{nan_idx,col} = subavg{nan_idx,col};
                            end
                        end
    
                        if ~isfield(category_data_concatSubs_allFormatsSources, cat_source)
                            category_data_concatSubs_allFormatsSources.(cat_source) = struct();
                        end
                        if ~isfield(category_data_concatSubs_allFormatsSources.(cat_source), cat_format)
                            category_data_concatSubs_allFormatsSources.(cat_source).(cat_format) = Y;
                            category_info.(cat_source).(cat_format).sub_ids = repmat({current_subject}, size(Y,1), 1);
                            category_info.(cat_source).(cat_format).modality_run = repmat({[current_modality,num2str(current_run)]}, size(Y,1), 1);
                        else
                            category_data_concatSubs_allFormatsSources.(cat_source).(cat_format) = [category_data_concatSubs_allFormatsSources.(cat_source).(cat_format); Y];
                            category_info.(cat_source).(cat_format).sub_ids = [category_info.(cat_source).(cat_format).sub_ids; repmat({current_subject}, size(Y,1), 1)];
                            category_info.(cat_source).(cat_format).modality_run = [category_info.(cat_source).(cat_format).modality_run; repmat({[current_modality,num2str(current_run)]}, size(Y,1), 1)];
                        end
                    end
                end
            end
        end
    end
end

function [va_data_all, subject_ids, modality_run_ids] = load_and_prepare_va_data(subjects, runs, spm1stlevel_dir, va_dir, va_sources, modality_to_use)

    va_data_all = struct();
    subject_ids = {};
    modality_run_ids = {};

    for s = 1:length(subjects)
        current_subject = subjects{s};
        modalities = get_modality_order(spm1stlevel_dir, current_subject, modality_to_use);
        
        % Load va data
        for m = 1:length(modalities)
            current_modality = modalities{m};
            for r = 1:length(runs)
                current_run = runs(r);
                va_file = fullfile(va_dir, ['sub', current_subject], ['sub', current_subject, '_', current_modality, '_va_run', num2str(current_run), '.mat']);
                load(va_file);

                for va_src = 1:length(va_sources)
                    va_source = va_sources{va_src};

                    %get a table of the va data with variable names startValence, endValence, startArousal, endArousal
                    Y = table(va.startValence.(va_source)', va.endValence.(va_source)', va.startArousal.(va_source)', va.endArousal.(va_source)',...
                            'VariableNames', {'startValence', 'endValence', 'startArousal', 'endArousal'});
                    subavg = table(va.startValence.Subavg', va.endValence.Subavg', va.startArousal.Subavg', va.endArousal.Subavg',...
                        'VariableNames', {'startValence', 'endValence', 'startArousal', 'endArousal'});
                    % fill NaNs in subspec with subavg
                    if strcmp(va_source, 'Subspec')
                        for col = 1:width(Y)
                            nan_idx = isnan(Y{:,col});
                            Y{nan_idx,col} = subavg{nan_idx,col};
                        end
                    end
                    if ~isfield(va_data_all, va_source)
                        va_data_all.(va_source) = Y;
                    else
                        va_data_all.(va_source) = [va_data_all.(va_source); Y];
                    end
                    if va_src == 1
                        subject_ids = [subject_ids; repmat({current_subject}, size(Y,1), 1)];
                        modality_run_ids = [modality_run_ids; repmat({[current_modality,num2str(current_run)]}, size(Y,1), 1)];
                    end
                end
            end
        end
    end
end

function region_results = perform_loso_encoding(brain_data_all, X_original, subject_ids, subjects, region_masks, region_names, save_nifti_files_yhat_ytest, save_nifti_files_performance, save_Yproj, output_dir_yhat_ytest, output_dir_performance, template_nifti_object, project_dir)

    region_results = table();
    if save_Yproj, Yproj_allregions = struct(); Ypredproj_allregions = struct(); end
    
    for r = 1:length(region_names)
        current_region_name = region_names{r};
        current_region_mask = region_masks{r};
        
        brain_data_region = brain_data_all.(current_region_name);
        if save_Yproj
            Yproj_allregions.(current_region_name) = table([], {}, 'VariableNames', {'Yproj','test_subject'}); 
            Ypredproj_allregions.(current_region_name) = table([], {}, 'VariableNames', {'Ypredproj','test_subject'});
        end

        for test_sub_idx = 1:length(subjects)
            test_subject = subjects{test_sub_idx};
            
            fprintf('    Testing on subject %s\n', test_subject);
            
            is_test = strcmp(subject_ids, test_subject);
            is_train = ~is_test;
            
            X_train = X_original(is_train, :);
            Y_train = brain_data_region(is_train, :);
            X_test = X_original(is_test, :);
            Y_test = brain_data_region(is_test, :);

            [~, ~, ~, ~, beta_cv] = plsregress(X_train, Y_train, 1);
            yhat_test_voxelwise = [ones(size(X_test, 1), 1), X_test] * beta_cv;

            b_slopes = beta_cv(2:end, :)'; % n_voxels x 1

            if save_Yproj
                new_entries = table(Y_test * b_slopes, repmat({test_subject}, size(Y_test, 1), 1), 'VariableNames', {'Yproj', 'test_subject'});
                Yproj_allregions.(current_region_name) = [Yproj_allregions.(current_region_name); new_entries];
                new_entries = table(yhat_test_voxelwise * b_slopes, repmat({test_subject}, size(yhat_test_voxelwise, 1), 1), 'VariableNames', {'Ypredproj', 'test_subject'});
                Ypredproj_allregions.(current_region_name) = [Ypredproj_allregions.(current_region_name); new_entries];
            end

            corrs_voxelwise_test = diag(corr(yhat_test_voxelwise, Y_test))'; % 1 x n_voxels
            
            region_result_row = table({test_subject}, mean(corrs_voxelwise_test), {current_region_name}, ...
                                        'VariableNames', {'test_subject', 'corr', 'region'});
            region_results = [region_results; region_result_row];
            
            if save_nifti_files_yhat_ytest || save_nifti_files_performance 
                save_nifti_files_func(yhat_test_voxelwise, Y_test, corrs_voxelwise_test, ...
                                    output_dir_yhat_ytest, output_dir_performance, ...
                                    current_region_name, current_region_mask, ...
                                    test_subject,  template_nifti_object, ...
                                    save_nifti_files_yhat_ytest, save_nifti_files_performance);
            end
        end
    end
    if save_Yproj, save(fullfile(output_dir_yhat_ytest, 'Yproj_allregions.mat'), 'Yproj_allregions'); save(fullfile(output_dir_yhat_ytest, 'Ypredproj_allregions.mat'), 'Ypredproj_allregions'); end
end

function region_results = perform_loso_encoding_noise_ceiling_resub(brain_data_all, X_original, subject_ids, subjects,region_names)

    region_results = table();
    
    for r = 1:length(region_names)
        current_region_name = region_names{r};
        
        brain_data_region = brain_data_all.(current_region_name);

        [~, ~, ~, ~, beta_allvoxels] = plsregress(X_original, brain_data_region, 1);
        Yhat = [ones(size(X_original, 1), 1), X_original] * beta_allvoxels;

        for test_sub_idx = 1:length(subjects)
            test_subject = subjects{test_sub_idx};
            
            fprintf('    Testing on subject %s\n', test_subject);
            
            is_test = strcmp(subject_ids, test_subject);

            Y_test = brain_data_region(is_test, :);

            corrs_voxelwise_test = diag(corr(Yhat(is_test, :), Y_test))'; % 1 x n_voxels
            region_result_row = table({test_subject}, mean(corrs_voxelwise_test), {current_region_name}, ...
                                        'VariableNames', {'test_subject', 'corr', 'region'});
            region_results = [region_results; region_result_row];
        end
    end
end

function region_results = perform_loso_encoding_noise_ceiling(brain_data_all, X_original, subject_ids, subjects,region_names)

    region_results = table();
    
    for r = 1:length(region_names)
        current_region_name = region_names{r};
        
        brain_data_region = brain_data_all.(current_region_name);

        for test_sub_idx = 1:length(subjects)
            test_subject = subjects{test_sub_idx};
            
            fprintf('    Testing on subject %s\n', test_subject);
            
            is_test = strcmp(subject_ids, test_subject);
            is_train = ~is_test;


            Y_test = brain_data_region(is_test, :);
            Y_train = brain_data_region(is_train, :);
            
            n_voxels = size(Y_test, 2);
            trainAvg_test_corr = zeros(1, n_voxels);
            
            for v = 1:n_voxels
                Y_train_v = Y_train(:, v);
                train_subject_ids = subject_ids(is_train);
                Y_train_v_subs = [];
                unique_train_subjects = unique(train_subject_ids);
                for s = 1:length(unique_train_subjects)
                    current_train_subject = unique_train_subjects{s};
                    is_train_subject = strcmp(train_subject_ids, current_train_subject);
                    Y_train_v_subs = [Y_train_v_subs, Y_train_v(is_train_subject)];
                end
                Y_train_v_avg_ts= mean(Y_train_v_subs, 2);
                trainAvg_test_corr(v) = corr(Y_train_v_avg_ts, Y_test(:, v));
            end

            region_result_row = table({test_subject}, mean(trainAvg_test_corr), {current_region_name}, ...
                                        'VariableNames', {'test_subject', 'corr', 'region'});
            region_results = [region_results; region_result_row];
        end
    end
end

function save_nifti_files_func(yhat_test, Y_test, corrs_test, output_dir_yhat_ytest, output_dir_performance, current_region_name, current_region_mask, test_subject, template_nifti_object, save_nifti_files_yhat_ytest, save_nifti_files_performance)

    if ~exist(output_dir_yhat_ytest, 'dir'), mkdir(output_dir_yhat_ytest); end
    if ~exist(output_dir_performance, 'dir'), mkdir(output_dir_performance); end
    temp_nifti_object = apply_mask(template_nifti_object, current_region_mask);

    if save_nifti_files_yhat_ytest
        temp_nifti_object.dat = yhat_test';
        temp_nifti_object.fullpath = fullfile(output_dir_yhat_ytest, ['yhat_test_sub', test_subject, '_', current_region_name, '.nii']);
        temp_nifti_object.write;

        temp_nifti_object.dat = Y_test';
        temp_nifti_object.fullpath = fullfile(output_dir_yhat_ytest, ['ytest_sub', test_subject, '_', current_region_name, '.nii']);
        temp_nifti_object.write;
    end

    if save_nifti_files_performance
        temp_nifti_object.dat = corrs_test';
        temp_nifti_object.fullpath = fullfile(output_dir_performance, ['loso_performance_test_sub', test_subject, '_', current_region_name, '.nii']);
        temp_nifti_object.write;
    end
end

function permutation_indices = generate_permutation_indices_loso(subject_ids, modality_run_ids, n_permutations, subjects)
    % Generate permutation indices for LOSO analysis

    permutation_indices = cell(n_permutations, 1);
    rng(42);
    for perm = 1:n_permutations
        perm_idx = 1:size(subject_ids, 1);
        for s = 1:length(subjects)
            subject_mask = strcmp(subject_ids, subjects{s});
            if any(subject_mask)
                subject_runs = unique(modality_run_ids(subject_mask));
                for r = 1:length(subject_runs)
                    run_mask = subject_mask & strcmp(modality_run_ids, subject_runs{r});
                    run_indices = find(run_mask);
                    if length(run_indices) > 1
                        shuffled_within_run = run_indices(randperm(length(run_indices)));
                        perm_idx(run_indices) = shuffled_within_run;
                    end
                end
            end
        end
        permutation_indices{perm} = perm_idx;
    end
end

function [p_values, null_distribution, permutation_fold_results] = run_permutation_test(brain_data_all, X_original, region_results, subject_ids, modality_run_ids,  subjects, region_masks, region_names, permutation_indices, shuffle_type, project_dir)

    n_permutations = length(permutation_indices);
    p_values = struct();

    % Get observed performance for each region
    observed_performance = struct();
    for r = 1:length(region_names)
        current_region_name = region_names{r};
        observed_performance.(current_region_name) = mean([region_results(strcmp(region_results.region, current_region_name), :).corr]); % mean across folds
    end

    % Run permutations
    null_distribution = struct();
    brain_data_all_perm = struct();
    for r = 1:length(region_names)
        current_region_name = region_names{r};
        null_distribution.(current_region_name) = nan(n_permutations, 1);

        if strcmp(shuffle_type, 'independent_shuffle_across_conditions')
            permutation_indices_region = generate_permutation_indices_loso(subject_ids, modality_run_ids, n_permutations, subjects);
        elseif strcmp(shuffle_type, 'same_shuffle_across_conditions')
            permutation_indices_region = permutation_indices;
        end
        for perm = 1:n_permutations
            perm_idx = permutation_indices_region{perm};
            brain_data_perm_region = brain_data_all.(current_region_name)(perm_idx, :);
            brain_data_all_perm.(current_region_name) = brain_data_perm_region;
        end
    end
    permutation_fold_results = table();
    for perm = 1:n_permutations
        if mod(perm, 100) == 0
            fprintf('      Permutation %d/%d\n', perm, n_permutations);
        end
        perm_performance_all_regions = perform_loso_encoding(brain_data_all_perm, X_original, subject_ids, subjects, region_masks, region_names, false, false, false, '', '', '', project_dir);
        perm_performance_all_regions.permutation = repmat(string(perm), size(perm_performance_all_regions, 1), 1);
        null_distribution.(current_region_name)(perm) = mean(perm_performance_all_regions.corr); % mean across folds
        permutation_fold_results = [permutation_fold_results; perm_performance_all_regions];
    end


    % Calculate p-values
    for r = 1:length(region_names)
        current_region_name = region_names{r};
        observed = observed_performance.(current_region_name);
        null_dist = null_distribution.(current_region_name);
        p_values.(current_region_name) = (sum(null_dist >= observed)+1) / (n_permutations+1);
    end 
end

function save_performance_results_wperm(region_results, p_values, null_distribution, csv_file_performance_wperm)
    % Save performance results with observed and permutation stats in one CSV file
    
    performance_table = table();
    unique_regions = unique(region_results.region);
    
    for r = 1:length(unique_regions)
        current_region = unique_regions{r};
        region_mask = strcmp(region_results.region, current_region);
        
        % Calculate observed performance (average across folds)
        observed_corr_avgfolds = mean([region_results(region_mask, :).corr]);
        
        % Get permutation stats
        p_value = p_values.(current_region);
        permutation_mean = mean(null_distribution.(current_region));
        permutation_std = std(null_distribution.(current_region));
        
        % Create row for this region
        region_row = table({current_region}, observed_corr_avgfolds, p_value, permutation_mean, permutation_std, ...
                          'VariableNames', {'region', 'observed_corr_avgfolds', 'p_value', 'permutation_mean', 'permutation_std'});
        
        performance_table = [performance_table; region_row];
    end
    
    writetable(performance_table, csv_file_performance_wperm, 'WriteMode', 'append');
end

function save_permutation_details(region_results, permutation_fold_results, csv_dir_perm)
    % Save detailed permutation results for each region in separate CSV files
    
    unique_regions = unique(region_results.region);
    
    for r = 1:length(unique_regions)
        current_region = unique_regions{r};
        
        csv_file_region = fullfile(csv_dir_perm, sprintf('loso_permutation_details_%s.csv', current_region));
        
        % Add observed data first
        observed_region_results = region_results(strcmp(region_results.region, current_region), :);
        observed_region_results.permutation = repmat({'observed'}, size(observed_region_results, 1), 1);
        
        % Add permutation results (fold-level)
        perm_region_results = permutation_fold_results(strcmp(permutation_fold_results.region, current_region), :);
        details_table = [observed_region_results; perm_region_results];
        
        writetable(details_table, csv_file_region, 'WriteMode', 'append');
    end
end
