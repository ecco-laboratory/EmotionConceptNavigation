[project_dir, fmri_data_dir, ~, theta_dir, ~, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();

ICA_type = 'noICA'; 
spm1stlevel_dir = [spm1stlevel_dir, '_', ICA_type, '_wButton'];
theta_sources = {'Subavg', 'Subspec'};
periodicity = {6, 4, 5, 7, 8};
cross_validations = {'xModalityRun', 'xModality', 'xRun'};%{'trainHalfTestSingleRun'};%{'xModalityRun', 'xModality', 'xRun','xModalitySingleRun','phi0','avgAllRuns','own'}; %'phi0''avgAllRuns',{'own'};%{'xModalitySingleRun'};%
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
subsample = false;
n_samples = 1000;seed = 42;

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

for t = 1:length(theta_sources)
    current_angle_source = theta_sources{t};
    [angle_data_all, subject_ids, modality_run_ids] = load_and_prepare_angle_data(subjects, runs, spm1stlevel_dir, theta_dir, current_angle_source);

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
                    elseif strcmp(current_cross_validation, 'own')
                        training_idx = {{'face1'}, {'face2'}, {'word1'}, {'word2'}};
                        test_idx = {{'face1'}, {'face2'}, {'word1'}, {'word2'}};
                    elseif strcmp(current_cross_validation, 'LORO')
                        %leave one run out
                        training_idx = {{'face2', 'word1','word2'}, {'face1', 'word1','word2'}, {'face1', 'face2','word2'}, {'face1', 'face2','word1'}};
                        test_idx = {{'face1'}, {'face2'}, {'word1'}, {'word2'}};
                    elseif strcmp(current_cross_validation, 'trainHalfTestSingleRun')
                        %train on half of the runs and test on the other half
                        training_idx = {{'face2', 'word1'}, {'face2', 'word2'}, {'word1', 'word2'}, ...
                                      {'face1', 'word1'}, {'face1', 'word2'}, {'word1', 'word2'}, ...
                                      {'face1', 'face2'}, {'face1', 'word2'}, {'face2', 'word2'}...
                                      {'face1', 'face2'}, {'face1', 'word1'}, {'face2', 'word1'}};
                        test_idx = {{'face1'}, {'face1'}, {'face1'}, ...
                                    {'face2'}, {'face2'}, {'face2'}, ...
                                    {'word1'}, {'word1'}, {'word1'}, ...
                                    {'word2'}, {'word2'}, {'word2'}};
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
                            if strcmp(current_cross_validation, 'xModalitySingleRun') || strcmp(current_cross_validation, 'own')
                                phi_value = phiRadDivByPeriod_struct.crossValidations.('xRun').(strjoin(current_training_idx, '_'));
                            elseif strcmp(current_cross_validation, 'trainHalfTestSingleRun')
                                %find in xModalityRun or xModality
                                if isfield(phiRadDivByPeriod_struct.crossValidations.xModalityRun, strjoin(current_training_idx, '_'))
                                    phi_value = phiRadDivByPeriod_struct.crossValidations.xModalityRun.(strjoin(current_training_idx, '_'));
                                elseif isfield(phiRadDivByPeriod_struct.crossValidations.xModality, strjoin(current_training_idx, '_'))
                                    phi_value = phiRadDivByPeriod_struct.crossValidations.xModality.(strjoin(current_training_idx, '_'));
                                else
                                    error('Cross validation not found: %s', strjoin(current_training_idx, '_'));
                                end
                            else
                                phi_value = phiRadDivByPeriod_struct.crossValidations.(current_cross_validation).(strjoin(current_training_idx, '_'));
                            end
                        end

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