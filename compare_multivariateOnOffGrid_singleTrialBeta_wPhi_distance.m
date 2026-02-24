[project_dir, fmri_data_dir, ~, theta_dir, category_dir, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();

ICA_type = 'noICA';
spm1stlevel_dir = [spm1stlevel_dir, '_', ICA_type, '_wButton'];
theta_sources = {'Subavg', 'Subspec'};
periodicity = {6, 4, 5, 7, 8};
cross_validations = {'xModalityRun', 'xModality', 'xRun'};%'phi0''avgAllRuns',
phi_averaging_methods = 'voxelComponentAverage'; 
phi_calculation_methods = 'singleTrialBeta';
phi_dir_voxelcircularmean = fullfile(phi_dir, 'univariate', phi_calculation_methods, ICA_type);

base_output_dir = fullfile(project_dir, 'outputs', 'singleTrialBetaAnalysis', ICA_type, 'incl_all_subs_trials','onoffGridcontrast_multivariate');
phi_dir_voxelcomponentaverage = fullfile(project_dir, 'outputs', 'phi', 'univariate', [phi_calculation_methods, '_regionwise'], ICA_type);


if ~exist(base_output_dir, 'dir'), mkdir(base_output_dir); end
runs = [1, 2]; 
smooth = 'unsmoothed'; base_output_dir = fullfile(base_output_dir, smooth); phi_dir_voxelcomponentaverage = fullfile(phi_dir_voxelcomponentaverage, smooth); spm1stlevel_dir = fullfile(spm1stlevel_dir, smooth);
%phi_source_region = 'OFC2016ConstantinescuR5'; base_output_dir = fullfile(base_output_dir, ['phi_', phi_source_region]);%'current'
phi_source_region= 'current';

brain_atlas = load_atlas('canlab2018');
region_names = {'OFC2016ConstantinescuR5','vmPFC2019BaoR5','vmPFC2016ConstantinescuR5','vmPFC2019Bao2R5','HC','ERC'};
region_masks = {fullfile(project_dir, 'masks', 'OFC_2016Constantinescu_r5.nii'),...
                fullfile(project_dir, 'masks', 'vmPFC_2019Bao_r5.nii'),...
                fullfile(project_dir, 'masks', 'vmPFC_2016Constantinescu_r5.nii'),...
                fullfile(project_dir, 'masks', 'vmPFC_2019Bao2_r5.nii'),...
                fullfile(project_dir, 'masks', 'vmPFC_2016Constantinescu_r5.nii'),...
                fullfile(project_dir, 'masks', 'HC_Julich.nii'),...
                fullfile(project_dir, 'masks', 'ERC_Julich.nii')};

[brain_data_all, brain_modality_run_idx, brain_sub_ids] = load_and_prepare_brain_data(subjects, spm1stlevel_dir, region_masks, region_names);
template_nifti_object = fmri_data(fullfile(spm1stlevel_dir, ['sub', subjects{1}], 'singleTrial', 'beta_0001.nii'));

va_dir = fullfile(project_dir, 'data', 'beh', 'VA');
inscan_judgment_dir = fullfile(project_dir, 'data', 'beh', 'inscan_judgment');
distance_types = {'va', 'judged', 'catCosine', 'catMax', 'v', 'a'};

category_formats = {'startContinuous', 'endContinuous'};
[category_data_concatSubs_allFormatsSources, category_info] = load_and_prepare_category_data(subjects, runs, spm1stlevel_dir, category_dir, theta_sources, category_formats); 
[va_data_allSources, subject_ids, modality_run_ids] = load_and_prepare_va_data(subjects, runs, spm1stlevel_dir, va_dir, theta_sources); 

if ~isequal(brain_modality_run_idx, modality_run_ids), disp(brain_modality_run_idx); disp(modality_run_ids); error('brain_modality_run_idx and modality_run_ids are not the same'); end
if ~isequal(brain_sub_ids, subject_ids), disp(brain_sub_ids); disp(subject_ids); error('brain_sub_ids and subject_ids are not the same'); end
if ~isequal(brain_modality_run_idx, category_info.(theta_sources{1}).(category_formats{1}).modality_run), disp(brain_modality_run_idx); disp(category_info.(theta_sources{1}).(category_formats{1}).modality_run);
    error('brain_modality_run_idx and category_info.(theta_sources{1}).(category_formats{1}).modality_run are not the same');
end
if ~isequal(brain_sub_ids, category_info.(theta_sources{1}).(category_formats{1}).sub_ids), disp(brain_sub_ids); disp(category_info.(theta_sources{1}).(category_formats{1}).sub_ids);
    error('brain_sub_ids and category_info.(theta_sources{1}).(category_formats{1}).sub_ids are not the same');
end
clear category_info;


for t = 1:length(theta_sources)
    current_angle_source = theta_sources{t};
    category_data_all = category_data_concatSubs_allFormatsSources.(current_angle_source); 
    va_data_all = va_data_allSources.(current_angle_source); 
    [consistency_data_all, comparison_data_all,~,~] = load_and_prepare_scanbeh_data(subjects, runs, spm1stlevel_dir, inscan_judgment_dir, current_angle_source); 
    [angle_data_all, subject_ids, modality_run_ids] = load_and_prepare_angle_data(subjects, runs, spm1stlevel_dir, theta_dir, current_angle_source);

    distance_data_alltypes = struct();
    for d = 1:length(distance_types)
        current_distance_type = distance_types{d};
        distance_data_alltypes.(current_distance_type) = prepare_distance_data(category_data_all, va_data_all, comparison_data_all, current_distance_type);
    end

    for p = 1:length(periodicity)
        current_periodicity = periodicity{p};
        fprintf('Processing periodicity %d...\n', current_periodicity);
        for s = 1:length(subjects)
            current_subject = subjects{s}; current_subject_idx = strcmp(subject_ids, current_subject);
            for r = 1:length(region_names)
                current_region_name = region_names{r};
                current_region_mask = region_masks{r};
                brain_data_current_region = brain_data_all.(current_region_name);
                
    
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
                    elseif strcmp(current_cross_validation, 'phi0')
                        training_idx = {'phi0', 'phi0', 'phi0'}; test_idx = {{'face1', 'word1', 'face2', 'word2'}, {'face1', 'face2'}, {'word1', 'word2'}};
                    elseif strcmp(current_cross_validation, 'avgAllRuns')
                        test_idx = {{'face1'},{'face2'},{'word1'},{'word2'},{'face1', 'face2', 'word1', 'word2'}, {'face1', 'face2'}, {'word1', 'word2'}};
                        training_idx = repmat({{'face1', 'face2', 'word1', 'word2'}}, 1, length(test_idx));
                    elseif strcmp(current_cross_validation, 'xModalitySingleRun')
                        training_idx = {{'face1'}, {'face1'}, {'word1'}, {'word1'}, {'face2'}, {'face2'}, {'word2'}, {'word2'}};
                        test_idx = {{'word1'}, {'word2'}, {'face1'}, {'face2'}, {'word1'}, {'word2'}, {'face1'}, {'face2'}};
                    end

                    for t_idx = 1:length(training_idx)
                        current_training_idx = training_idx{t_idx}; current_test_idx = test_idx{t_idx};
                    
                
                        if strcmp(current_cross_validation, 'phi0')
                            phi_value = 0;
                        else
                            if strcmp(phi_source_region, 'current')
                                phi_file = fullfile(phi_dir_voxelcomponentaverage, ['sub', current_subject], 'includeNonNorm', ['periodicity', num2str(current_periodicity)], ['theta', current_angle_source],...
                                    ['PhiRadDivByPeriod_', current_region_name, '.mat']);
                            else
                                phi_file = fullfile(phi_dir_voxelcomponentaverage, ['sub', current_subject], 'includeNonNorm', ['periodicity', num2str(current_periodicity)], ['theta', current_angle_source],...
                                    ['PhiRadDivByPeriod_', phi_source_region, '.mat']);
                            end
                            load(phi_file);
                            if strcmp(current_cross_validation, 'xModalitySingleRun')
                                phi_value = phiRadDivByPeriod_struct.crossValidations.('xRun').(strjoin(current_training_idx, '_'));
                            else
                                phi_value = phiRadDivByPeriod_struct.crossValidations.(current_cross_validation).(strjoin(current_training_idx, '_'));
                            end
                        end
                        
                        data_mask = current_subject_idx & ismember(modality_run_ids, current_test_idx);
                        distance_data_alltypes_masked = struct();
                        for d = 1:length(distance_types)
                            current_distance_type = distance_types{d};
                            distance_data = distance_data_alltypes.(current_distance_type);
                            distance_data_alltypes_masked.(current_distance_type) = distance_data(data_mask, :);
                        end
                        T = get_pattern_similarity_wDistance(brain_data_current_region(data_mask, :), ...
                                                                    angle_data_all(data_mask, :), phi_value, current_periodicity, distance_data_alltypes_masked,...
                                                                    consistency_data_all(data_mask, :), comparison_data_all(data_mask, :));
                        %new table rows
                        %add column subject and region and test runs to T
                        T.subject = repmat(string(current_subject), size(T, 1), 1);
                        T.region = repmat(string(current_region_name), size(T, 1), 1);
                        T.test_runs = repmat(string(strjoin(current_test_idx, '_')), size(T, 1), 1);
                        
                        %save csv files
                        path_to_save = fullfile(base_output_dir, current_angle_source, ['periodicity', num2str(current_periodicity)], current_cross_validation, ['phi_', phi_averaging_methods], phi_calculation_methods, 'csv');
                        if ~exist(path_to_save, 'dir'), mkdir(path_to_save); end
                        csv_file_contrastvalues = fullfile(path_to_save, 'pattern_similarity_wDistance.csv'); 
                        if exist(csv_file_contrastvalues, 'file'), writetable(T, csv_file_contrastvalues, 'WriteMode', 'append'); else writetable(T, csv_file_contrastvalues); end

                    end
                end   
            end
        end
    end
end


function [brain_data_all, brain_modality_run_idx, brain_sub_ids] = load_and_prepare_brain_data(subjects, spm1stlevel_dir, region_masks, region_names)
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
    % Load brain data
    brain_data = fmri_data(beta_images);
    for r = 1:length(region_masks)
        current_region_mask = region_masks{r}; current_region_name = region_names{r};
        brain_data_region = apply_mask(brain_data, current_region_mask).dat';
        brain_data_all.(current_region_name) = brain_data_region;
    end
end

function [angle_data_all, subject_ids, modality_run_ids] = load_and_prepare_angle_data(subjects, runs, spm1stlevel_dir, theta_dir, current_angle_source)

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
        
        if contains(first_stim_trial, 'face') && contains(last_stim_trial, 'word')
            modalities = {'face', 'word'};
        elseif contains(first_stim_trial, 'word') && contains(last_stim_trial, 'face')
            modalities = {'word', 'face'};
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

function [distance_data] = prepare_distance_data(category_data_all, va_data_all, comparison_data_all, current_distance_type)
    vDist = abs(va_data_all.endValence - va_data_all.startValence);
    aDist = abs(va_data_all.endArousal - va_data_all.startArousal);
    vaDist = sqrt(vDist.^2 + aDist.^2);
    if strcmp(current_distance_type, 'va')
        distance_data = vaDist - mean(vaDist);
    elseif strcmp(current_distance_type, 'v')
        distance_data = vDist - mean(vDist);
    elseif strcmp(current_distance_type, 'a')
        distance_data = aDist - mean(aDist);
    elseif strcmp(current_distance_type, 'judged')
        judgedDist = zeros(size(comparison_data_all, 1), 1);
        for t = 1:size(comparison_data_all, 1)
            if strcmpi(comparison_data_all{t}, 'PLEASANTNESS')
                judgedDist(t) = vDist(t);
            elseif strcmpi(comparison_data_all{t}, 'ACTIVATION')
                judgedDist(t) = aDist(t);
            else
                error('Unexpected comparison_data value at trial %d', t);
            end
        end
        distance_data = judgedDist - mean(judgedDist);
    elseif strcmp(current_distance_type, 'catCosine') || strcmp(current_distance_type, 'catMax')
        startMat = table2array(category_data_all.('startContinuous')); endMat = table2array(category_data_all.('endContinuous')); 
        nTrials  = size(startMat,1);
        cosineDist = zeros(nTrials,1); catMaxDist = zeros(nTrials,1);
        [~, startMaxIdx] = max(startMat, [], 2); [~, endMaxIdx]   = max(endMat, [], 2);
        for t = 1:nTrials
            s = startMat(t,:); e = endMat(t,:);
            cosineDist(t) = 1 - (dot(s,e) / (norm(s)*norm(e) + eps));
            catMaxDist(t) = mean([abs(s(startMaxIdx(t)) - e(startMaxIdx(t))), abs(s(endMaxIdx(t))   - e(endMaxIdx(t)))]);
        end
        if strcmp(current_distance_type, 'catCosine')
            distance_data = cosineDist - mean(cosineDist);
        elseif strcmp(current_distance_type, 'catMax')
            distance_data = catMaxDist - mean(catMaxDist);
        end
    end
end

function modalities = get_modality_order(spm1stlevel_dir, subject)
    singleTrial_dir = fullfile(spm1stlevel_dir, ['sub', subject], 'singleTrial');
        load(fullfile(singleTrial_dir, 'SPM.mat'));
        spm_regressor_names = SPM.xX.name;
        stim_trials = contains(spm_regressor_names, 'trial');
        stim_trials_idx = find(stim_trials);
        
        first_stim_trial = spm_regressor_names(find(stim_trials, 1, 'first'));
        last_stim_trial = spm_regressor_names(find(stim_trials, 1, 'last'));
        
        if contains(first_stim_trial, 'face') && contains(last_stim_trial, 'word')
            modalities = {'face', 'word'};
        elseif contains(first_stim_trial, 'word') && contains(last_stim_trial, 'face')
            modalities = {'word', 'face'};
        else
            error('Modality order unclear.');
        end
end

function [consistency_data_all, comparison_data_all, subject_ids, modality_run_ids] = load_and_prepare_scanbeh_data(subjects, runs, spm1stlevel_dir, scanbeh_dir, current_beh_source)

    consistency_data_all = [];
    comparison_data_all = [];
    subject_ids = {};
    modality_run_ids = {};

    for s = 1:length(subjects)
        current_subject = subjects{s};
        modalities = get_modality_order(spm1stlevel_dir, current_subject);
        
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


function [category_data_concatSubs_allFormatsSources, category_info] = load_and_prepare_category_data(subjects, runs, spm1stlevel_dir, category_dir, category_sources, category_formats)

    category_data_concatSubs_allFormatsSources = struct();
    category_info = struct();

    for s = 1:length(subjects)
        current_subject = subjects{s};
        modalities = get_modality_order(spm1stlevel_dir, current_subject);
        
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

function [va_data_all, subject_ids, modality_run_ids] = load_and_prepare_va_data(subjects, runs, spm1stlevel_dir, va_dir, va_sources)

    va_data_all = struct();
    subject_ids = {};
    modality_run_ids = {};

    for s = 1:length(subjects)
        current_subject = subjects{s};
        modalities = get_modality_order(spm1stlevel_dir, current_subject);
        
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


function [bin_centers, bin_types] = getGridBinTypes(aligned_angle_data, periodicity)
    %deal with float-point error
    tol = 10 * eps(max(abs(aligned_angle_data)));

    % classify angles into: 1 = on-grid, -1 = off-grid
    bin_centers = nan(size(aligned_angle_data));
    bin_types = nan(size(aligned_angle_data));

    [on_grid_bins, off_grid_bins, on_centers, off_centers] = get_grid_bins_with_edges(periodicity);
    %the first two on_grid_bins are share the same bin_center, so add one more bin_center
    on_centers = [on_centers(1), on_centers];
    
    % assign on-grid = 1
    for i = 1:size(on_grid_bins, 1)
        lower = on_grid_bins(i, 1);
        upper = on_grid_bins(i, 2);
        idx = aligned_angle_data >= lower - tol & aligned_angle_data < upper + tol & isnan(bin_types);
        bin_types(idx) = 1;
        bin_centers(idx) = on_centers(i);
    end
    
    % assign off-grid = -1
    for i = 1:size(off_grid_bins, 1)
        lower = off_grid_bins(i, 1);
        upper = off_grid_bins(i, 2);
        idx = aligned_angle_data >= lower - tol & aligned_angle_data < upper + tol & isnan(bin_types);
        bin_types(idx) = -1;
        bin_centers(idx) = off_centers(i);
    end
    
end
   

function [on_grid_bins, off_grid_bins, on_centers, off_centers] = get_grid_bins_with_edges(periodicity)
    
        bin_width = 360 / (2 * periodicity);                    
        num_bins = 2 * periodicity;
    
        grid_angles = mod((0:periodicity-1) * 360/periodicity, 360);  
        on_centers = grid_angles;                                        
        off_centers = mod(grid_angles + 360/(2 * periodicity), 360);           
    
        on_grid_bins = expand_bins(on_centers, bin_width);
        off_grid_bins = expand_bins(off_centers, bin_width);
end
    
function bins = expand_bins(centers, width)
    
    half_width = width / 2;
    bins = [];

    for i = 1:length(centers)
        start_angle = mod(centers(i) - half_width, 360);
        end_angle = mod(centers(i) + half_width, 360);

        if start_angle < end_angle
            bins = [bins; start_angle, end_angle];
        else
            bins = [bins; start_angle, 360; 0, end_angle];
        end
    end
end
    

function T = get_pattern_similarity_wDistance(singleTrial_betas, angle_data, phi_value,current_periodicity, distance_data_alltypes_masked, consistency_data_all, comparison_data_all)

        phi_deg = rad2deg(phi_value); % convert phi to degrees

        %define bin width and centers
        period_deg = 360 / current_periodicity; %60 for 6-fold
        bin_width = period_deg/2; %30 for 6-fold
        bin_half = bin_width/2; %15 for 6-fold
        bin_centers = 0:bin_width:(180-bin_width); 
        align_half = bin_half/2; %7.5 for 6-fold

        %find aligned trials
        angle_modP = mod(angle_data - phi_deg, period_deg);
        aligned_trials = (abs(angle_modP) < align_half) | (abs(angle_modP - period_deg) < align_half); %deal with wrap-around (like 59-degree should belong to 0 degree bin for 6 fold)
        aligned_idx = find(aligned_trials); 

        nTrials = size(singleTrial_betas, 1);

        % compute correlations for aligned × others
        T = table();


        for i = 1:numel(aligned_idx)
            t = aligned_idx(i); 
            beta_t = singleTrial_betas(t,:)'; %column vector for trial t

            for j = 1:nTrials
                if j == t, continue; end %skip self

                % angular difference modulo 60
                d_raw = abs(angle_data(t) - angle_data(j)); d_on_period = mod(d_raw, period_deg);
                diff_to_0 = min(d_on_period, period_deg - d_on_period);
                diff_to_half = abs(d_on_period - period_deg/2); diff_to_half = min(diff_to_half, period_deg - diff_to_half);
                % classify correlation
                if diff_to_0 < bin_half
                    alignment_label_curr_pair = 1;
                elseif diff_to_half < bin_half
                    alignment_label_curr_pair = 0;
                else
                    alignment_label_curr_pair = nan;
                end

                %6-bin similarity
                d180 = mod(d_raw, 180);
                for b = 1:length(bin_centers)
                    diff_bincenter = abs(d180 - bin_centers(b)); diff_bincenter = min(diff_bincenter, 180 - diff_bincenter);
                    if diff_bincenter < bin_half
                        bin_curr_pair = b;
                        break;
                    end
                end
                new_row = table(corr(beta_t, singleTrial_betas(j,:)'), alignment_label_curr_pair, bin_curr_pair, 'VariableNames', {'corr','alignment_label','bin'});
                clear alignment_label_curr_pair bin_curr_pair;

                distance_types = fieldnames(distance_data_alltypes_masked);
                for d = 1:length(distance_types)
                    dtype = distance_types{d};
                    new_row.(['phi_aligned_' dtype]) = distance_data_alltypes_masked.(dtype)(t);
                    new_row.(['other_' dtype])       = distance_data_alltypes_masked.(dtype)(j);
                end
                new_row.('phi_aligned_consistency') = consistency_data_all(t);
                new_row.('phi_aligned_comparison') = comparison_data_all(t);
                new_row.('other_consistency') = consistency_data_all(j);
                new_row.('other_comparison') = comparison_data_all(j);
                T = [T; new_row];
            end
        end
end