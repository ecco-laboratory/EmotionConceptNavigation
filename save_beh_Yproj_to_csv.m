[project_dir, fmri_data_dir, ~, theta_dir, category_dir, ~, ~, ~, spm_timing_dir, motion_regressors_dir, nuisance_regressors_dir, spm1stlevel_dir, ~, phi_dir, subjects] = set_up_dirs_constants();

ICA_type = 'noICA'; % 'noICA' or 'wICA'
spm1stlevel_dir = [spm1stlevel_dir, '_', ICA_type, '_wButton'];
beh_sources = {'Subavg','Subspec'};
smooth = 'unsmoothed'; spm1stlevel_dir = fullfile(spm1stlevel_dir, smooth);

modality_to_use = 'both_modalities';
runs = [1, 2];
distance_types = {'catMaxEuclidean','catDominantAvg','catDominantEuclidean','vSigned','aSigned','judgedSigned','va', 'judged', 'catCosine', 'catMax', 'v', 'a'};
inscan_judgment_dir = fullfile(project_dir, 'data', 'beh', 'inscan_judgment');


region_names = {'aHC','pHC','OFC2016ConstantinescuR5','vmPFCcurrentStudyR5','HC', 'ERC'};
yhat_output_dir = fullfile(project_dir, 'outputs', 'singleTrialBetaAnalysis', ICA_type, 'incl_all_subs_trials', 'distanceCentered_multivariate_encoding_yhat_ytest', smooth, modality_to_use); 

va_dir = fullfile(project_dir, 'data', 'beh', 'VA');
category_formats = {'startContinuous', 'endContinuous'};
[category_data_concatSubs_allFormatsSources, category_info] = load_and_prepare_category_data(subjects, runs, spm1stlevel_dir, category_dir, beh_sources, category_formats, modality_to_use); 
[va_data_allSources, subject_ids, modality_run_ids] = load_and_prepare_va_data(subjects, runs, spm1stlevel_dir, va_dir, beh_sources, modality_to_use); 

for t = 1:length(beh_sources)
    current_beh_source = beh_sources{t};
    [consistency_data_all, comparison_data_all,subject_ids, modality_run_ids] = load_and_prepare_scanbeh_data(subjects, runs, spm1stlevel_dir, inscan_judgment_dir, current_beh_source, modality_to_use); 
    category_data_all = category_data_concatSubs_allFormatsSources.(current_beh_source); 
    va_data_all = va_data_allSources.(current_beh_source); 
    for c = 1:length(distance_types)
        current_distance_type = distance_types{c};
        distance_data = prepare_distance_data(category_data_all, va_data_all, comparison_data_all, current_distance_type);
        
        current_yhat_output_dir = fullfile(yhat_output_dir, current_beh_source, current_distance_type);
        output_dir = fullfile(current_yhat_output_dir, 'csv');if ~exist(output_dir, 'dir'), mkdir(output_dir); end
        load(fullfile(current_yhat_output_dir, 'Ypredproj_allregions.mat'));
        for r = 1:length(region_names)
            current_region_name = region_names{r};
            current_region_Yproj = Ypredproj_allregions.(current_region_name);
            %make sure current_region_Yproj .test_subject is the same as the subject_ids
            if ~isequal(current_region_Yproj.test_subject, subject_ids)
                error('Test subject mismatch between current_region_Yproj and subject_ids.');
            end
            tbl = current_region_Yproj;
            %tbl.Yproj = tbl.Yproj - mean(tbl.Yproj);
            tbl.Ypredproj = tbl.Ypredproj;
            tbl.consistency = consistency_data_all;
            %save tbl to csv
            %writetable(tbl, fullfile(output_dir, ['beh_Yproj_', current_region_name, '.csv']));

            tbl.distance = distance_data;
            writetable(tbl, fullfile(output_dir, ['beh_Ypredproj_raw_wDistance_', current_region_name, '.csv']));

        end
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

function [distance_data] = prepare_distance_data(category_data_all, va_data_all, comparison_data_all, current_distance_type)
    vSignedDist = va_data_all.endValence - va_data_all.startValence;
    aSignedDist = va_data_all.endArousal - va_data_all.startArousal;
    vDist = abs(vSignedDist);
    aDist = abs(aSignedDist);
    vaDist = sqrt(vDist.^2 + aDist.^2);
    
    if strcmp(current_distance_type, 'va')
        distance_data = vaDist - mean(vaDist);
    elseif strcmp(current_distance_type, 'v')
        distance_data = vDist - mean(vDist);
    elseif strcmp(current_distance_type, 'a')
        distance_data = aDist - mean(aDist);
    elseif strcmp(current_distance_type, 'vSigned')
        distance_data = vSignedDist - mean(vSignedDist);
    elseif strcmp(current_distance_type, 'aSigned')
        distance_data = aSignedDist - mean(aSignedDist);
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
        if strcmp(current_distance_type, 'judged'), distance_data = judgedDist - mean(judgedDist); elseif strcmp(current_distance_type, 'judgedSigned'), distance_data = signed_judgedDist - mean(signed_judgedDist); end
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
            distance_data = cosineDist - mean(cosineDist);
        elseif strcmp(current_distance_type, 'catMax')
            distance_data = catMaxDist - mean(catMaxDist);
        elseif strcmp(current_distance_type, 'catMaxEuclidean')
            distance_data = catMaxEuclideanDist - mean(catMaxEuclideanDist);
        elseif strcmp(current_distance_type, 'catDominantEuclidean')
            distance_data = catDominantEuclideanDist - mean(catDominantEuclideanDist);
        elseif strcmp(current_distance_type, 'catDominantAvg')
            distance_data = catDominantAvgDist - mean(catDominantAvgDist);
        end
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