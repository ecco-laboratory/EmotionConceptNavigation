[project_dir, fmri_data_dir, ~, theta_dir, category_dir, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();


ICA_type = 'noICA'; 
spm1stlevel_dir = [spm1stlevel_dir, '_', ICA_type, '_wButton'];
theta_sources = {'Subavg', 'Subspec'};
periodicity = {6, 4, 5, 7, 8};
cross_validations = {'xModalitySingleRun'};%{'xModalityRun', 'xModality', 'xRun'};%'phi0''avgAllRuns',
phi_averaging_methods = 'voxelComponentAverage'; 
phi_calculation_methods = 'singleTrialBeta';
phi_dir_voxelcircularmean = fullfile(phi_dir, 'univariate', phi_calculation_methods, ICA_type);

base_output_dir = fullfile(project_dir, 'outputs', 'singleTrialBetaAnalysis', ICA_type, 'incl_all_subs_trials','onoffGridcontrast_multivariate');
phi_dir_voxelcomponentaverage = fullfile(project_dir, 'outputs', 'phi', 'univariate', [phi_calculation_methods, '_regionwise'], ICA_type);


if ~exist(base_output_dir, 'dir'), mkdir(base_output_dir); end
runs = [1, 2]; 
%
smooth = 'unsmoothed'; base_output_dir = fullfile(base_output_dir, smooth); phi_dir_voxelcomponentaverage = fullfile(phi_dir_voxelcomponentaverage, smooth); spm1stlevel_dir = fullfile(spm1stlevel_dir, smooth);
%phi_source_region = 'OFC2016ConstantinescuR5'; base_output_dir = fullfile(base_output_dir, ['phi_', phi_source_region]);%'current'
phi_source_region= 'current';
subsample = false;
n_samples = 1000;seed = 42;

category_formats = {'start', 'end'};
if_va_control = true;
if_category_control = true;
if_id_setting_control = true;
va_dir = fullfile(project_dir, 'data', 'beh', 'VA');
id_setting_dir = fullfile(project_dir, 'data', 'beh', 'id_setting');

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

if if_category_control, [category_data_concatSubs_allFormatsSources, category_info] = load_and_prepare_category_data(subjects, runs, spm1stlevel_dir, category_dir, theta_sources, category_formats); end
if if_va_control, [va_data_all, subject_ids, modality_run_ids] = load_and_prepare_va_data(subjects, runs, spm1stlevel_dir, va_dir, theta_sources); end
if if_va_control
    if ~isequal(brain_modality_run_idx, modality_run_ids), disp(brain_modality_run_idx); disp(modality_run_ids); error('brain_modality_run_idx and modality_run_ids are not the same'); end
    if ~isequal(brain_sub_ids, subject_ids), disp(brain_sub_ids); disp(subject_ids); error('brain_sub_ids and subject_ids are not the same'); end
    clear subject_ids; clear modality_run_ids;
end
if if_category_control
    if ~isequal(brain_modality_run_idx, category_info.(theta_sources{1}).(category_formats{1}).modality_run), disp(brain_modality_run_idx); disp(category_info.(theta_sources{1}).(category_formats{1}).modality_run);
        error('brain_modality_run_idx and category_info.(theta_sources{1}).(category_formats{1}).modality_run are not the same');
    end
    if ~isequal(brain_sub_ids, category_info.(theta_sources{1}).(category_formats{1}).sub_ids), disp(brain_sub_ids); disp(category_info.(theta_sources{1}).(category_formats{1}).sub_ids);
        error('brain_sub_ids and category_info.(theta_sources{1}).(category_formats{1}).sub_ids are not the same');
    end
    clear category_info;
end


for t = 1:length(theta_sources)
    current_angle_source = theta_sources{t};
    [angle_data_all, subject_ids, modality_run_ids] = load_and_prepare_angle_data(subjects, runs, spm1stlevel_dir, theta_dir, current_angle_source);
    if if_category_control, category_data_all = category_data_concatSubs_allFormatsSources.(current_angle_source); end
    if if_id_setting_control, [id_setting_data_all,~,~] = load_and_prepare_id_setting_data(subjects, runs, spm1stlevel_dir, id_setting_dir, current_angle_source); end


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

                        if if_va_control
                            data_mask = current_subject_idx & ismember(modality_run_ids, current_test_idx);
                            core_args = {brain_data_current_region(data_mask, :), angle_data_all(data_mask, :), phi_value, current_periodicity, ...
                                        va_data_all.(current_angle_source)(data_mask, :)};

                            optional_args = {};
                            if if_category_control, optional_args = [optional_args, {'StartCategory', category_data_all.start(data_mask, :), 'EndCategory', category_data_all.end(data_mask, :)}];end
                            if if_id_setting_control, optional_args = [optional_args, {'IDSetting', id_setting_data_all(data_mask, :)}]; end
                            optional_args = [optional_args, {'Subsample', subsample, 'NumSamples', n_samples, 'Seed', seed}];

                            [beta_alignment, contrast_value] = perform_multivariate_contrast_control(core_args{:}, optional_args{:});

                            if subsample == false
                                new_beta_alignment_row = table(contrast_value, beta_alignment, string(current_subject), string(current_region_name), string(strjoin(current_test_idx, '_')), ...
                                                        'VariableNames', {'contrast_value', 'beta_alignment', 'subject', 'region', 'test_runs'});
                            else
                                new_beta_alignment_row = table(contrast_value.mean, contrast_value.median, ...
                                                            contrast_value.std, contrast_value.ci95(1), contrast_value.ci95(2),...% {contrast_value.samples}, ...
                                                            beta_alignment.mean, beta_alignment.std, beta_alignment.ci95(1), beta_alignment.ci95(2),...% {beta_alignment.samples}, ...
                                                            string(current_subject), string(current_region_name), string(strjoin(current_test_idx, '_')), ...
                                                            'VariableNames', {'contrast_mean', 'contrast_median', ...
                                                                'contrast_std', 'ci_lower', 'ci_upper',...% 'contrast_samples', ...
                                                                'beta_alignment_mean', 'beta_alignment_std', 'beta_alignment_ci_lower', 'beta_alignment_ci_upper',...% 'beta_alignment_samples', ...
                                                                'subject', 'region', 'test_runs'});
                            end
                            path_to_save = fullfile(base_output_dir, current_angle_source, ['periodicity', num2str(current_periodicity)], current_cross_validation, ['phi_', phi_averaging_methods], phi_calculation_methods, 'csv');
                            if ~exist(path_to_save, 'dir'), mkdir(path_to_save); end

                            %get filename suffix based on control factors
                            suffix_parts = {};
                            if if_va_control, suffix_parts{end+1} = 'va'; end
                            if if_category_control, suffix_parts{end+1} = 'cat'; end
                            if if_id_setting_control, suffix_parts{end+1} = 'id'; end
                            if isequal(suffix_parts, {'va'}), suffix_str = ''; else suffix_str = ['_' strjoin(suffix_parts, '_')]; end

                            if subsample == false, csv_file_beta_alignment = fullfile(path_to_save, ['beta_alignment' suffix_str '.csv']); else csv_file_beta_alignment = fullfile(path_to_save, ['beta_alignment_subsample' suffix_str '.csv']); end
                            if exist(csv_file_beta_alignment, 'file'), writetable(new_beta_alignment_row, csv_file_beta_alignment, 'WriteMode', 'append'); else writetable(new_beta_alignment_row, csv_file_beta_alignment); end
                        else
                            [contrast_value, bin_means, bin_centers, bin_trial_counts, num_aligned_trials] = perform_multivariate_contrast(brain_data_current_region(current_subject_idx & ismember(modality_run_ids, current_test_idx), :), ...
                                                                        angle_data_all(current_subject_idx & ismember(modality_run_ids, current_test_idx), :), phi_value, current_periodicity, subsample, n_samples, seed);

                            
                            if subsample == false
                                %store bin means for this subject, region, and test runs
                                new_bin_means_row = table(bin_means', bin_centers', bin_trial_counts', repmat(string(current_subject), length(bin_means), 1), repmat(string(current_region_name), length(bin_means), 1), repmat(string(strjoin(current_test_idx, '_')), length(bin_means), 1), ...
                                                        'VariableNames', {'bin_means', 'bin_centers', 'bin_trial_counts', 'subject', 'region', 'test_runs'});
                                %store contrast value for this subject, region, and test runs
                                new_contrast_value_row = table(contrast_value, num_aligned_trials, string(current_subject), string(current_region_name), string(strjoin(current_test_idx, '_')), ...
                                                        'VariableNames', {'contrast_value', 'num_aligned_trials', 'subject', 'region', 'test_runs'});
                            else
                                new_contrast_value_row = table(contrast_value.mean, contrast_value.median, ...
                                                            contrast_value.std, contrast_value.ci95(1), contrast_value.ci95(2),...% {contrast_value.samples}, ...
                                                            num_aligned_trials, string(current_subject), string(current_region_name), string(strjoin(current_test_idx, '_')), ...
                                                    'VariableNames', {'contrast_mean', 'contrast_median', ...
                                                                    'contrast_std', 'ci_lower', 'ci_upper',...% 'contrast_samples', ...
                                                                    'num_aligned_trials', 'subject', 'region', 'test_runs'});
                            end
                            
                            %save csv files
                            path_to_save = fullfile(base_output_dir, current_angle_source, ['periodicity', num2str(current_periodicity)], current_cross_validation, ['phi_', phi_averaging_methods], phi_calculation_methods, 'csv');
                            if ~exist(path_to_save, 'dir'), mkdir(path_to_save); end
                            if subsample == false
                                csv_file_binmeans = fullfile(path_to_save, 'bin_means.csv'); 
                                if exist(csv_file_binmeans, 'file'), writetable(new_bin_means_row, csv_file_binmeans, 'WriteMode', 'append'); else writetable(new_bin_means_row, csv_file_binmeans); end
                            end
                            if subsample, csv_file_contrastvalues = fullfile(path_to_save, 'contrast_values_subsample.csv'); else csv_file_contrastvalues = fullfile(path_to_save, 'contrast_values.csv'); end
                            if exist(csv_file_contrastvalues, 'file'), writetable(new_contrast_value_row, csv_file_contrastvalues, 'WriteMode', 'append'); else writetable(new_contrast_value_row, csv_file_contrastvalues); end
                        end
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

function [angle_data_all, subject_ids, modality_run_ids] = load_and_prepare_angle_data(subjects, runs, spm1stlevel_dir, theta_dir, current_angle_source)

    angle_data_all = [];
    subject_ids = {};
    modality_run_ids = {};

    for s = 1:length(subjects)
        current_subject = subjects{s};
        modalities = get_modality_order(spm1stlevel_dir, current_subject);
        
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
    
                        Y = categories.(cat_format).(cat_source)';
                        subavg = categories.(cat_format).Subavg';
                        % fill NaNs in subspec with subavg
                        if strcmp(cat_source, 'Subspec')
                            nan_idx = cellfun(@(x) isnumeric(x) && isnan(x), Y);
                            Y(nan_idx) = subavg(nan_idx);
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

function [id_setting_all, subject_ids, modality_run_ids] = load_and_prepare_id_setting_data(subjects, runs, spm1stlevel_dir, id_setting_dir, id_setting_source)

    id_setting_all = {};
    subject_ids = {};
    modality_run_ids = {};

    for s = 1:length(subjects)
        current_subject = subjects{s};
        modalities = get_modality_order(spm1stlevel_dir, current_subject);
        
        % Load id setting data
        for m = 1:length(modalities)
            current_modality = modalities{m};
            for r = 1:length(runs)
                current_run = runs(r);
                id_setting_file = fullfile(id_setting_dir, ['sub', current_subject], ['sub', current_subject, '_', current_modality, '_idsetting_run', num2str(current_run), '.mat']);
                load(id_setting_file);

                idsetting_data = id_setting.(id_setting_source)';
                subavg = id_setting.Subavg';
                % fill NaNs in subspec with subavg
                if strcmp(id_setting_source, 'Subspec')
                    nan_idx = cellfun(@(x) isnumeric(x) && isnan(x), idsetting_data);
                    idsetting_data(nan_idx) = subavg(nan_idx);
                end
                id_setting_all = [id_setting_all; idsetting_data];
                subject_ids = [subject_ids; repmat({current_subject}, size(idsetting_data, 1), 1)];
                modality_run_ids = [modality_run_ids; repmat({[current_modality,num2str(current_run)]}, size(idsetting_data, 1), 1)];
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
    

function [contrast_value, bin_means, bin_centers, bin_trial_counts, num_aligned_trials] = perform_multivariate_contrast(singleTrial_betas, angle_data, phi_value,current_periodicity, subsample, n_samples, seed)

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
        aligned_idx = find(aligned_trials); num_aligned_trials = length(aligned_idx);

        nTrials = size(singleTrial_betas, 1);

        % compute correlations for aligned × others
        aligned_corrs = [];
        misaligned_corrs = [];
        bin_corrs = cell(1, length(bin_centers));
        for b = 1:length(bin_centers), bin_corrs{b} = []; end

        for i = 1:numel(aligned_idx)
            t = aligned_idx(i); 
            beta_t = singleTrial_betas(t,:)'; %column vector for trial t

            for j = 1:nTrials
                if j == t, continue; end %skip self

                r = corr(beta_t, singleTrial_betas(j,:)');

                % angular difference modulo 60
                d_raw = abs(angle_data(t) - angle_data(j)); d_on_period = mod(d_raw, period_deg);
                diff_to_0 = min(d_on_period, period_deg - d_on_period);
                diff_to_half = abs(d_on_period - period_deg/2); diff_to_half = min(diff_to_half, period_deg - diff_to_half);
                % classify correlation
                if diff_to_0 < bin_half
                    aligned_corrs(end+1) = r;
                elseif diff_to_half < bin_half
                    misaligned_corrs(end+1) = r;
                end

                %6-bin similarity
                d180 = mod(d_raw, 180);
                for b = 1:length(bin_centers)
                    diff_bincenter = abs(d180 - bin_centers(b)); diff_bincenter = min(diff_bincenter, 180 - diff_bincenter);
                    if diff_bincenter < bin_half
                        bin_corrs{b} = [bin_corrs{b}; r];
                        break;
                    end
                end


            end
        end

        %optional subsampling for balanced comparison 
        rng(seed); 
        if subsample
            n_min = min(numel(aligned_corrs), numel(misaligned_corrs));
            for s = 1:n_samples
                a_idx = randperm(numel(aligned_corrs), n_min);
                m_idx = randperm(numel(misaligned_corrs), n_min);
                contrast_samples(s) = mean(aligned_corrs(a_idx), 'omitnan') - mean(misaligned_corrs(m_idx), 'omitnan');
            end
            contrast_value.mean = mean(contrast_samples, 'omitnan');
            contrast_value.median = median(contrast_samples, 'omitnan');
            contrast_value.std = std(contrast_samples, 'omitnan');
            contrast_value.ci95 = prctile(contrast_samples, [2.5, 97.5]);
            contrast_value.samples = contrast_samples;
        else
            %contrast (aligned vs misaligned)
            contrast_value = mean(aligned_corrs, 'omitnan') - mean(misaligned_corrs, 'omitnan');
        end

        %mean corrs per bin
        bin_means = nan(1, length(bin_centers)); bin_trial_counts = nan(1, length(bin_centers));
        for b = 1:length(bin_centers)
            if ~isempty(bin_corrs{b})
                bin_means(b) = mean(bin_corrs{b} ,'omitnan');
                bin_trial_counts(b) = length(bin_corrs{b});
            end
        end
end



function [beta_alignment, contrast_value] = perform_multivariate_contrast_control(singleTrial_betas, angle_data, phi_value,current_periodicity, va_data, varargin)

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
    aligned_corrs = [];
    misaligned_corrs = [];
    bin_corrs = cell(1, length(bin_centers));
    for b = 1:length(bin_centers), bin_corrs{b} = []; end

    % get control factors
    v_start = va_data.startValence; a_start = va_data.startArousal; v_end = va_data.endValence; a_end = va_data.endArousal;
    %optional control factors
    p = inputParser;
    addParameter(p, 'StartCategory', []);addParameter(p, 'EndCategory', []);addParameter(p, 'IDSetting', []);
    addParameter(p, 'Subsample', false); addParameter(p, 'NumSamples', 1000); addParameter(p, 'Seed', 42);
    parse(p, varargin{:});
    startCat = p.Results.StartCategory; endCat = p.Results.EndCategory; idSet = p.Results.IDSetting;
    subsample = p.Results.Subsample; n_samples = p.Results.NumSamples; seed = p.Results.Seed;
    
    R_all = []; alignment_all = []; Dstart_all = []; Dend_all = [];
    same_start = []; same_end = []; same_id = [];
    
    %Compare aligned x aligned with aligned x misaligned
    for i = 1:numel(aligned_idx)
        t = aligned_idx(i); 
        beta_t = singleTrial_betas(t,:)'; %column vector for trial t

        for j = 1:nTrials
            if j == t, continue; end %skip self

            r = corr(beta_t, singleTrial_betas(j,:)');

            % angular difference modulo 60
            d_raw = abs(angle_data(t) - angle_data(j)); d_on_period = mod(d_raw, period_deg);
            diff_to_0 = min(d_on_period, period_deg - d_on_period);
            diff_to_half = abs(d_on_period - period_deg/2); diff_to_half = min(diff_to_half, period_deg - diff_to_half);
            % classify correlation
            if diff_to_0 < bin_half
                aligned_corrs(end+1) = r;
                alignment_label = 1;
            elseif diff_to_half < bin_half
                misaligned_corrs(end+1) = r;
                alignment_label = 0;
            end

            % calculate control factors
            D_start = sqrt((v_start(t)-v_start(j))^2 + (a_start(t)-a_start(j))^2);
            D_end   = sqrt((v_end(t)-v_end(j))^2     + (a_end(t)-a_end(j))^2);

            R_all(end+1) = r; alignment_all(end+1) = alignment_label;
            Dstart_all(end+1) = D_start; Dend_all(end+1) = D_end;

            % optional control factors
            if ~isempty(startCat), same_start(end+1,1) = double(strcmp(startCat{t}, startCat{j})); end
            if ~isempty(endCat), same_end(end+1,1) = double(strcmp(endCat{t}, endCat{j})); end
            if ~isempty(idSet), same_id(end+1,1) = double(strcmp(idSet{t}, idSet{j})); end

        end
    end
    R_all = R_all(:); alignment_all = alignment_all(:); Dstart_all = Dstart_all(:); Dend_all = Dend_all(:);
    if ~isempty(startCat), same_start = same_start(:); end 
    if ~isempty(endCat), same_end = same_end(:); end
    if ~isempty(idSet), same_id = same_id(:); end

   %regression with control: R ~ alignment + D_start + D_end + optional control factors (startCat, endCat, idSet)
   if ~isempty(R_all)
        if subsample
            rng(seed);
            aligned_idx_all = find(alignment_all == 1);
            misaligned_idx_all = find(alignment_all == 0);
            n_min = min(numel(aligned_idx_all), numel(misaligned_idx_all));

            beta_samples = nan(n_samples,1);
            contrast_samples = nan(n_samples,1);
            for s = 1:n_samples
                a_idx = randsample(aligned_idx_all, n_min);
                m_idx = randsample(misaligned_idx_all, n_min);
                subs_idx = [a_idx; m_idx];

                %subsampled regression
                tbl_sub = create_table_subset(R_all, alignment_all, Dstart_all, Dend_all, ...
                                            same_start, same_end, same_id, startCat, endCat, idSet, subs_idx);
                beta_samples(s) = run_lm_and_get_beta(tbl_sub);

                %subsampled contrast
                a_corrs = aligned_corrs(randperm(numel(aligned_corrs), n_min));
                m_corrs = misaligned_corrs(randperm(numel(misaligned_corrs), n_min));
                contrast_samples(s) = mean(a_corrs, 'omitnan') - mean(m_corrs, 'omitnan');
            end

            % beta_alignment summary
            beta_alignment.mean   = mean(beta_samples, 'omitnan');
            beta_alignment.median = median(beta_samples, 'omitnan');
            beta_alignment.std    = std(beta_samples, 'omitnan');
            beta_alignment.ci95   = prctile(beta_samples, [2.5, 97.5]);
            beta_alignment.samples = beta_samples;

            % contrast summary
            contrast_value.mean   = mean(contrast_samples, 'omitnan');
            contrast_value.median = median(contrast_samples, 'omitnan');
            contrast_value.std    = std(contrast_samples, 'omitnan');
            contrast_value.ci95   = prctile(contrast_samples, [2.5, 97.5]);
            contrast_value.samples = contrast_samples;
        else
            tbl = create_table_subset(R_all, alignment_all, Dstart_all, Dend_all, ...
                                    same_start, same_end, same_id, startCat, endCat, idSet, 1:length(R_all));
            beta_alignment = run_lm_and_get_beta(tbl);
            contrast_value = mean(aligned_corrs, 'omitnan') - mean(misaligned_corrs, 'omitnan');
        end
    else
        beta_alignment = NaN;
        contrast_value = NaN;
    end
end

function tbl = create_table_subset(R_all, alignment_all, Dstart_all, Dend_all, ...
                                   same_start, same_end, same_id, startCat, endCat, idSet, idx)
    tbl = table(R_all(idx), alignment_all(idx), Dstart_all(idx), Dend_all(idx), ...
                'VariableNames', {'R_all','alignment_all','Dstart_all','Dend_all'});
    tbl.Dstart_all = zscore(tbl.Dstart_all);
    tbl.Dend_all   = zscore(tbl.Dend_all);
    if ~isempty(startCat), tbl.sameStart = same_start(idx); end
    if ~isempty(endCat), tbl.sameEnd = same_end(idx); end
    if ~isempty(idSet), tbl.sameID = same_id(idx); end
end

function beta_alignment = run_lm_and_get_beta(tbl)
    predictors = setdiff(tbl.Properties.VariableNames, {'R_all'});
    formula = sprintf('R_all ~ %s', strjoin(predictors, ' + '));
    mdl = fitlm(tbl, formula);
    beta_alignment = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Properties.RowNames, 'alignment_all'));
end
