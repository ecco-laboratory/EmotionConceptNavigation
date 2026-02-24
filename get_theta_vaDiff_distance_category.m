clear; close all; clc;
[project_dir, ~, psychopy_csv_dir, theta_dir, category_dir, vaDiff_dir, avgVA_dir, distance_dir, ~, ~, ~,~, ~, ~, subjects] = set_up_dirs_constants();

overwrite = true;
runs = [1, 2];
modalities = {'face', 'word'};
category_names = {'hap','fear','sad','surprise','disgust','neutral','anger'};
cate_ori_rating_names = struct( ...
    'HAP', 'hap', ...
    'AFR', 'fear', ...
    'SAD', 'sad', ...
    'SUR', 'surprise', ...
    'DIS', 'disgust', ...
    'NEU', 'neutral', ...
    'ANG', 'anger');

num_trials_face = 84;
num_trials_word = 124;
if ~exist(theta_dir, 'dir')
    mkdir(theta_dir);
end
if ~exist(vaDiff_dir, 'dir')
    mkdir(vaDiff_dir);
end
if ~exist(avgVA_dir, 'dir')
    mkdir(avgVA_dir);
end
if ~exist(category_dir, 'dir')
    mkdir(category_dir);
end
if ~exist(distance_dir, 'dir')
    mkdir(distance_dir);
end


for i = 1:length(subjects)
    subject = subjects{i};
    %check if the file already exists
    if ~overwrite && exist(fullfile(theta_dir, ['sub', subject], ['sub', subject, '_face_thetas_run1.mat']), 'file') && ...
       exist(fullfile(theta_dir, ['sub', subject], ['sub', subject, '_face_thetas_run2.mat']), 'file') && ...
       exist(fullfile(theta_dir, ['sub', subject], ['sub', subject, '_word_thetas_run1.mat']), 'file') && ...
       exist(fullfile(theta_dir, ['sub', subject], ['sub', subject, '_word_thetas_run2.mat']), 'file') 
       
        disp(['Skipping subject ', subject, ' because files already exist']);
        continue;
    end
    for m = 1:length(modalities)
        modality = modalities{m};
        for r = 1:length(runs)
            run = runs(r);
            disp(['Processing subject: ', subject, ', run: ', num2str(run), ', modality: ', modality]);
            
            scanner_beh_files = dir(fullfile(psychopy_csv_dir, ['sub-', subject], ['sub-', subject, '_task-', modality, '_run-0', num2str(run), '_*.csv']));
            if isempty(scanner_beh_files) 
                warning(['Could not find scanner behavioral files for subject ', subject]);
                continue;
            end
            scanner_beh_file = readtable(fullfile(scanner_beh_files(1).folder, scanner_beh_files(1).name));
            % Load norm data
            norm_data_file = fullfile(psychopy_csv_dir, 'behTables',['norm_table_', modality, '.csv']);
            norm_data = readtable(norm_data_file);
            % Load subject-specific ratings
            subject_specific_ratings_file = fullfile(psychopy_csv_dir, 'behTables', ['sub', subject], ['sub', subject, '_table_', modality, '.csv']);
            if exist(subject_specific_ratings_file, 'file')
                subject_specific_ratings = readtable(subject_specific_ratings_file);
            else
                warning(['Could not find subject-specific ratings file for subject ', subject, ', so using sub0001 as template and filling with nans']);
                sub0001_data = readtable(fullfile(psychopy_csv_dir, 'behTables', 'sub0001', ['sub0001_table_', modality, '.csv']));
                subject_specific_ratings = array2table(nan(size(sub0001_data, 1), length(sub0001_data.Properties.VariableNames)), 'VariableNames', sub0001_data.Properties.VariableNames);
            end
            % Load subject-averaged ratings
            subject_averaged_ratings_file = fullfile(psychopy_csv_dir, 'behTables', ['subAvg_table_', modality, '.csv']);
            subject_averaged_ratings = readtable(subject_averaged_ratings_file);
            
            thetas_norm = [];
            thetas_subspec = [];
            thetas_subavg = [];
            valenceDiff_norm = [];
            valenceDiff_subspec = [];
            valenceDiff_subavg = [];
            arousalDiff_norm = [];
            arousalDiff_subspec = [];
            arousalDiff_subavg = [];
            avgValence_norm = [];
            avgValence_subspec = [];
            avgValence_subavg = [];
            avgArousal_norm = [];
            avgArousal_subspec = [];
            avgArousal_subavg = [];
            distances_norm = [];
            distances_subspec = [];
            distances_subavg = [];
            startCategory_norm = {};
            endCategory_norm = {};
            startCategory_subspec = {};
            endCategory_subspec = {};
            startCategory_subavg = {};
            endCategory_subavg = {};
            
            startCategoryOnehot_norm = table();
            endCategoryOnehot_norm = table();
            startCategoryOnehot_subspec = table();
            endCategoryOnehot_subspec = table();
            startCategoryContinuous_subspec = table();
            endCategoryContinuous_subspec = table();
            startCategoryOnehot_subavg = table();
            endCategoryOnehot_subavg = table();
            startCategoryContinuous_subavg = table();
            endCategoryContinuous_subavg = table();
            
            if strcmp(modality, 'face')
                % For face, find match using face_id, start_category, and end_category
                for trial = 1:size(scanner_beh_file, 1)
                    if isempty(scanner_beh_file.group{trial}) 
                        continue;  % Skip this trial
                    end

                    
                    face_id = scanner_beh_file.group{trial};
                    start_category = scanner_beh_file.start_category{trial};
                    end_category = scanner_beh_file.end_category{trial};
                    %Match in norm data
                    norm_row = find(strcmp(norm_data.face_id, face_id) & ...
                                    strcmp(norm_data.start_label, start_category) & ...
                                    strcmp(norm_data.end_label, end_category));
                    if ~isempty(norm_row)
                        thetas_norm(end+1) = norm_data.theta(norm_row);
                        valenceDiff_norm(end+1) = norm_data.end_valence(norm_row) - norm_data.start_valence(norm_row);
                        arousalDiff_norm(end+1) = norm_data.end_arousal(norm_row) - norm_data.start_arousal(norm_row);
                        avgValence_norm(end+1) = (norm_data.start_valence(norm_row) + norm_data.end_valence(norm_row)) / 2;
                        avgArousal_norm(end+1) = (norm_data.start_arousal(norm_row) + norm_data.end_arousal(norm_row)) / 2;
                        distances_norm(end+1) = norm_data.distance(norm_row);
                        startCategory_norm{end+1} = start_category;
                        endCategory_norm{end+1} = end_category;
                        %make one-hot tables with variable names category_names(start_category being 1 and 0 for all other categories)
                        %note start_category and category_names have this correspondence: cate_ori_rating_names = {'hap': 'HAP', 'fear': 'AFR', 'sad': 'SAD', 'surprise': 'SUR', 'anger': 'ANG', 'neutral': 'NEU', 'disgust': 'DIS'}
                        start_onehot = double(ismember(category_names, cate_ori_rating_names.(start_category)));
                        end_onehot = double(ismember(category_names, cate_ori_rating_names.(end_category)));    
                        startCategoryOnehot_norm = [startCategoryOnehot_norm; array2table(start_onehot, 'VariableNames', category_names)];
                        endCategoryOnehot_norm = [endCategoryOnehot_norm; array2table(end_onehot, 'VariableNames', category_names)];
                    else
                        thetas_norm(end+1) = NaN;
                        valenceDiff_norm(end+1) = NaN;
                        arousalDiff_norm(end+1) = NaN;
                        avgValence_norm(end+1) = NaN;
                        avgArousal_norm(end+1) = NaN;
                        distances_norm(end+1) = NaN;
                        startCategory_norm{end+1} = NaN;
                        endCategory_norm{end+1} = NaN;
                        start_onehot = NaN(1, length(category_names));
                        end_onehot = NaN(1, length(category_names));
                        startCategoryOnehot_norm = [startCategoryOnehot_norm; array2table(start_onehot, 'VariableNames', category_names)];
                        endCategoryOnehot_norm = [endCategoryOnehot_norm; array2table(end_onehot, 'VariableNames', category_names)];
                        warning(['No matching norm data for face_id: ', face_id, ', start_category: ', start_category, ', end_category: ', end_category]);
                    end
                    
                    %Match in subject-specific ratings
                    subspec_row = find(strcmp(subject_specific_ratings.group, face_id) & ...
                                      strcmp(subject_specific_ratings.ori_start_category, start_category) & ...
                                      strcmp(subject_specific_ratings.ori_end_category, end_category));
                    if ~isempty(subspec_row)
                        thetas_subspec(end+1) = subject_specific_ratings.theta(subspec_row);
                        valenceDiff_subspec(end+1) = subject_specific_ratings.end_valence(subspec_row) - subject_specific_ratings.start_valence(subspec_row);
                        arousalDiff_subspec(end+1) = subject_specific_ratings.end_arousal(subspec_row) - subject_specific_ratings.start_arousal(subspec_row);
                        avgValence_subspec(end+1) = (subject_specific_ratings.start_valence(subspec_row) + subject_specific_ratings.end_valence(subspec_row)) / 2;
                        avgArousal_subspec(end+1) = (subject_specific_ratings.start_arousal(subspec_row) + subject_specific_ratings.end_arousal(subspec_row)) / 2;
                        distances_subspec(end+1) = subject_specific_ratings.distance(subspec_row);
                        startCategory_subspec{end+1} = subject_specific_ratings.max_start_category{subspec_row};
                        endCategory_subspec{end+1} = subject_specific_ratings.max_end_category{subspec_row};

                        start_vars = strcat('start_', category_names);
                        end_vars = strcat('end_', category_names);
                        start_data = subject_specific_ratings{subspec_row, start_vars};
                        end_data = subject_specific_ratings{subspec_row, end_vars};
                        startCategoryContinuous_subspec = [startCategoryContinuous_subspec; array2table(start_data, 'VariableNames', category_names)];
                        endCategoryContinuous_subspec = [endCategoryContinuous_subspec; array2table(end_data, 'VariableNames', category_names)];
                    else
                        thetas_subspec(end+1) = NaN;
                        valenceDiff_subspec(end+1) = NaN;
                        arousalDiff_subspec(end+1) = NaN;
                        avgValence_subspec(end+1) = NaN;
                        avgArousal_subspec(end+1) = NaN;
                        distances_subspec(end+1) = NaN;
                        startCategory_subspec{end+1} = NaN;
                        endCategory_subspec{end+1} = NaN;
                        start_data = NaN(1, length(category_names));
                        end_data = NaN(1, length(category_names));
                        startCategoryContinuous_subspec = [startCategoryContinuous_subspec; array2table(start_data, 'VariableNames', category_names)];
                        endCategoryContinuous_subspec = [endCategoryContinuous_subspec; array2table(end_data, 'VariableNames', category_names)];
                        warning(['No matching subject-specific data for face_id: ', face_id, ', start_category: ', start_category, ', end_category: ', end_category]);
                    end

                    %Match in subject-averaged ratings
                    subavg_row = find(strcmp(subject_averaged_ratings.group, face_id) & ...
                                     strcmp(subject_averaged_ratings.ori_start_category, start_category) & ...
                                     strcmp(subject_averaged_ratings.ori_end_category, end_category));
                    if ~isempty(subavg_row)
                        thetas_subavg(end+1) = subject_averaged_ratings.theta(subavg_row);
                        valenceDiff_subavg(end+1) = subject_averaged_ratings.end_valence(subavg_row) - subject_averaged_ratings.start_valence(subavg_row);
                        arousalDiff_subavg(end+1) = subject_averaged_ratings.end_arousal(subavg_row) - subject_averaged_ratings.start_arousal(subavg_row);
                        avgValence_subavg(end+1) = (subject_averaged_ratings.start_valence(subavg_row) + subject_averaged_ratings.end_valence(subavg_row)) / 2;
                        avgArousal_subavg(end+1) = (subject_averaged_ratings.start_arousal(subavg_row) + subject_averaged_ratings.end_arousal(subavg_row)) / 2;
                        distances_subavg(end+1) = subject_averaged_ratings.distance(subavg_row);
                        startCategory_subavg{end+1} = subject_averaged_ratings.max_start_category{subavg_row};
                        endCategory_subavg{end+1} = subject_averaged_ratings.max_end_category{subavg_row};

                        start_vars = strcat('start_', category_names);
                        end_vars = strcat('end_', category_names);
                        start_data = subject_averaged_ratings{subavg_row, start_vars};
                        end_data = subject_averaged_ratings{subavg_row, end_vars};
                        startCategoryContinuous_subavg = [startCategoryContinuous_subavg; array2table(start_data, 'VariableNames', category_names)];
                        endCategoryContinuous_subavg = [endCategoryContinuous_subavg; array2table(end_data, 'VariableNames', category_names)];
                    else
                        thetas_subavg(end+1) = NaN;
                        valenceDiff_subavg(end+1) = NaN;
                        arousalDiff_subavg(end+1) = NaN;
                        avgValence_subavg(end+1) = NaN;
                        avgArousal_subavg(end+1) = NaN;
                        distances_subavg(end+1) = NaN;
                        startCategory_subavg{end+1} = NaN;
                        endCategory_subavg{end+1} = NaN;
                        start_data = NaN(1, length(category_names));
                        end_data = NaN(1, length(category_names));
                        startCategoryContinuous_subavg = [startCategoryContinuous_subavg; array2table(start_data, 'VariableNames', category_names)];
                        endCategoryContinuous_subavg = [endCategoryContinuous_subavg; array2table(end_data, 'VariableNames', category_names)];
                        warning(['No matching subject-averaged data for face_id: ', face_id, ', start_category: ', start_category, ', end_category: ', end_category]);
                    end
                end

                
                if ~isempty(thetas_norm)
                   %check length of thetas_norm, thetas_subspec, thetas_subavg
                   if length(scanner_beh_file.group(2:end-1)) ~= length(thetas_norm) || length(scanner_beh_file.group(2:end-1)) ~= length(thetas_subspec) || length(scanner_beh_file.group(2:end-1)) ~= length(thetas_subavg)
                       warning(['Length of scanner_beh_file(2:end-1) is not the same as the length of thetas_norm, thetas_subspec, thetas_subavg for subject ', subject]);
                       continue;
                   end
                end
                
            elseif strcmp(modality, 'word')
                % For word, find match using group, start_word, and end_word
                for trial = 1:size(scanner_beh_file, 1)
                    if isempty(scanner_beh_file.group{trial})
                        continue;  % Skip this trial
                    end
                    
                    group = scanner_beh_file.group{trial};
                    start_word = scanner_beh_file.start_word{trial};
                    end_word = scanner_beh_file.end_word{trial};
                    start_category = scanner_beh_file.start_category{trial};
                    end_category = scanner_beh_file.end_category{trial};
                    if isempty(start_word) || isempty(end_word)
                        continue;
                    end
                    
                    % Match in norm data
                    norm_row = find(strcmp(norm_data.group, group) & ...
                                    strcmp(norm_data.start_stimulus, start_word) & ...
                                    strcmp(norm_data.end_stimulus, end_word));
                    if ~isempty(norm_row)
                        thetas_norm(end+1) = norm_data.theta(norm_row);
                        valenceDiff_norm(end+1) = norm_data.end_valence(norm_row) - norm_data.start_valence(norm_row);
                        arousalDiff_norm(end+1) = norm_data.end_arousal(norm_row) - norm_data.start_arousal(norm_row);
                        avgValence_norm(end+1) = (norm_data.start_valence(norm_row) + norm_data.end_valence(norm_row)) / 2;
                        avgArousal_norm(end+1) = (norm_data.start_arousal(norm_row) + norm_data.end_arousal(norm_row)) / 2;
                        distances_norm(end+1) = norm_data.distance(norm_row);
                        startCategory_norm{end+1} = start_category;
                        endCategory_norm{end+1} = end_category;
                        start_onehot = double(ismember(category_names, cate_ori_rating_names.(start_category)));
                        end_onehot = double(ismember(category_names, cate_ori_rating_names.(end_category)));
                        startCategoryOnehot_norm = [startCategoryOnehot_norm; array2table(start_onehot, 'VariableNames', category_names)];
                        endCategoryOnehot_norm = [endCategoryOnehot_norm; array2table(end_onehot, 'VariableNames', category_names)];
                    else
                        thetas_norm(end+1) = NaN;
                        valenceDiff_norm(end+1) = NaN;
                        arousalDiff_norm(end+1) = NaN;
                        avgValence_norm(end+1) = NaN;
                        avgArousal_norm(end+1) = NaN;
                        distances_norm(end+1) = NaN;
                        startCategory_norm{end+1} = NaN;
                        endCategory_norm{end+1} = NaN;
                        start_onehot = NaN(1, length(category_names));
                        end_onehot = NaN(1, length(category_names));
                        startCategoryOnehot_norm = [startCategoryOnehot_norm; array2table(start_onehot, 'VariableNames', category_names)];
                        endCategoryOnehot_norm = [endCategoryOnehot_norm; array2table(end_onehot, 'VariableNames', category_names)];
                        warning(['No matching norm data for group: ', group, ', start_word: ', start_word, ', end_word: ', end_word]);
                    end
                    
                    % Match in subject-specific ratings
                    subspec_row = find(strcmp(subject_specific_ratings.group, group) & ...
                                      strcmp(subject_specific_ratings.start_stimulus, start_word) & ...
                                      strcmp(subject_specific_ratings.end_stimulus, end_word));
                    if ~isempty(subspec_row)
                        thetas_subspec(end+1) = subject_specific_ratings.theta(subspec_row);
                        valenceDiff_subspec(end+1) = subject_specific_ratings.end_valence(subspec_row) - subject_specific_ratings.start_valence(subspec_row);
                        arousalDiff_subspec(end+1) = subject_specific_ratings.end_arousal(subspec_row) - subject_specific_ratings.start_arousal(subspec_row);
                        avgValence_subspec(end+1) = (subject_specific_ratings.start_valence(subspec_row) + subject_specific_ratings.end_valence(subspec_row)) / 2;
                        avgArousal_subspec(end+1) = (subject_specific_ratings.start_arousal(subspec_row) + subject_specific_ratings.end_arousal(subspec_row)) / 2;
                        distances_subspec(end+1) = subject_specific_ratings.distance(subspec_row);
                        startCategory_subspec{end+1} = subject_specific_ratings.max_start_category{subspec_row};
                        endCategory_subspec{end+1} = subject_specific_ratings.max_end_category{subspec_row};

                        start_vars = strcat('start_', category_names);
                        end_vars = strcat('end_', category_names);
                        start_data = subject_specific_ratings{subspec_row, start_vars};
                        end_data = subject_specific_ratings{subspec_row, end_vars};
                        startCategoryContinuous_subspec = [startCategoryContinuous_subspec; array2table(start_data, 'VariableNames', category_names)];
                        endCategoryContinuous_subspec = [endCategoryContinuous_subspec; array2table(end_data, 'VariableNames', category_names)];

                    else
                        thetas_subspec(end+1) = NaN;
                        valenceDiff_subspec(end+1) = NaN;
                        arousalDiff_subspec(end+1) = NaN;
                        avgValence_subspec(end+1) = NaN;
                        avgArousal_subspec(end+1) = NaN;
                        distances_subspec(end+1) = NaN;
                        startCategory_subspec{end+1} = NaN;
                        endCategory_subspec{end+1} = NaN;
                        startCategoryContinuous_subspec = [startCategoryContinuous_subspec; array2table(start_data, 'VariableNames', category_names)];
                        endCategoryContinuous_subspec = [endCategoryContinuous_subspec; array2table(end_data, 'VariableNames', category_names)];
                        warning(['No matching subject-specific data for group: ', group, ', start_word: ', start_word, ', end_word: ', end_word]);
                    end
                    
                    % Match in subject-averaged ratings
                    subavg_row = find(strcmp(subject_averaged_ratings.group, group) & ...
                                     strcmp(subject_averaged_ratings.start_stimulus, start_word) & ...
                                     strcmp(subject_averaged_ratings.end_stimulus, end_word));
                    if ~isempty(subavg_row)
                        thetas_subavg(end+1) = subject_averaged_ratings.theta(subavg_row);
                        valenceDiff_subavg(end+1) = subject_averaged_ratings.end_valence(subavg_row) - subject_averaged_ratings.start_valence(subavg_row);
                        arousalDiff_subavg(end+1) = subject_averaged_ratings.end_arousal(subavg_row) - subject_averaged_ratings.start_arousal(subavg_row);
                        avgValence_subavg(end+1) = (subject_averaged_ratings.start_valence(subavg_row) + subject_averaged_ratings.end_valence(subavg_row)) / 2;
                        avgArousal_subavg(end+1) = (subject_averaged_ratings.start_arousal(subavg_row) + subject_averaged_ratings.end_arousal(subavg_row)) / 2;
                        distances_subavg(end+1) = subject_averaged_ratings.distance(subavg_row);
                        startCategory_subavg{end+1} = subject_averaged_ratings.max_start_category{subavg_row};
                        endCategory_subavg{end+1} = subject_averaged_ratings.max_end_category{subavg_row};

                        start_vars = strcat('start_', category_names);
                        end_vars = strcat('end_', category_names);
                        start_data = subject_averaged_ratings{subavg_row, start_vars};
                        end_data = subject_averaged_ratings{subavg_row, end_vars};
                        startCategoryContinuous_subavg = [startCategoryContinuous_subavg; array2table(start_data, 'VariableNames', category_names)];
                        endCategoryContinuous_subavg = [endCategoryContinuous_subavg; array2table(end_data, 'VariableNames', category_names)];
                    else
                        thetas_subavg(end+1) = NaN;
                        valenceDiff_subavg(end+1) = NaN;
                        arousalDiff_subavg(end+1) = NaN;
                        avgValence_subavg(end+1) = NaN;
                        avgArousal_subavg(end+1) = NaN;
                        distances_subavg(end+1) = NaN;
                        startCategory_subavg{end+1} = NaN;
                        endCategory_subavg{end+1} = NaN;
                        start_data = NaN(1, length(category_names));
                        end_data = NaN(1, length(category_names));
                        startCategoryContinuous_subavg = [startCategoryContinuous_subavg; array2table(start_data, 'VariableNames', category_names)];
                        endCategoryContinuous_subavg = [endCategoryContinuous_subavg; array2table(end_data, 'VariableNames', category_names)];
                        warning(['No matching subject-averaged data for group: ', group, ', start_word: ', start_word, ', end_word: ', end_word]);
                    end
                end
                

                if ~isempty(thetas_norm)
                   %check length of thetas_norm, thetas_subspec, thetas_subavg
                   if length(scanner_beh_file.group(2:end-1)) ~= length(thetas_norm) || length(scanner_beh_file.group(2:end-1)) ~= length(thetas_subspec) || length(scanner_beh_file.group(2:end-1)) ~= length(thetas_subavg)
                       warning(['Length of scanner_beh_file(2:end-1) is not the same as the length of thetas_norm, thetas_subspec, thetas_subavg for subject ', subject]);
                       continue;
                   end
                end
            end
            
            thetas = struct();
            thetas.Norm = thetas_norm;
            thetas.Subspec = thetas_subspec;
            thetas.Subavg = thetas_subavg;
            theta_file = fullfile(theta_dir, ['sub', subject], ['sub', subject, '_', modality, '_thetas_run', num2str(run), '.mat']);
            if ~exist(fullfile(theta_dir, ['sub', subject]), 'dir')
                mkdir(fullfile(theta_dir, ['sub', subject]));
            end
            save(theta_file, 'thetas');
            disp(['Saved thetas to ', theta_file]);
            
            vaDiffs = struct();
            vaDiffs.valenceDiff.Norm = valenceDiff_norm;
            vaDiffs.valenceDiff.Subspec = valenceDiff_subspec;
            vaDiffs.valenceDiff.Subavg = valenceDiff_subavg;
            vaDiffs.arousalDiff.Norm = arousalDiff_norm;
            vaDiffs.arousalDiff.Subspec = arousalDiff_subspec;
            vaDiffs.arousalDiff.Subavg = arousalDiff_subavg;
            vaDiff_file = fullfile(vaDiff_dir, ['sub', subject], ['sub', subject, '_', modality, '_vaDiffs_run', num2str(run), '.mat']);
            if ~exist(fullfile(vaDiff_dir, ['sub', subject]), 'dir')
                mkdir(fullfile(vaDiff_dir, ['sub', subject]));
            end
            save(vaDiff_file, 'vaDiffs');
            disp(['Saved vaDiffs to ', vaDiff_file]);

            avgVAs = struct();
            avgVAs.valence.Norm = avgValence_norm;
            avgVAs.valence.Subspec = avgValence_subspec;
            avgVAs.valence.Subavg = avgValence_subavg;
            avgVAs.arousal.Norm = avgArousal_norm;
            avgVAs.arousal.Subspec = avgArousal_subspec;
            avgVAs.arousal.Subavg = avgArousal_subavg;
            avgVA_file = fullfile(avgVA_dir, ['sub', subject], ['sub', subject, '_', modality, '_avgVAs_run', num2str(run), '.mat']);
            if ~exist(fullfile(avgVA_dir, ['sub', subject]), 'dir')
                mkdir(fullfile(avgVA_dir, ['sub', subject]));
            end
            save(avgVA_file, 'avgVAs');
            disp(['Saved avgVAs to ', avgVA_file]);

            distances = struct();
            distances.euclideanDistance.Norm = distances_norm;
            distances.euclideanDistance.Subspec = distances_subspec;
            distances.euclideanDistance.Subavg = distances_subavg;
            distances.valenceAbsDistance.Norm = abs(valenceDiff_norm);
            distances.valenceAbsDistance.Subspec = abs(valenceDiff_subspec);
            distances.valenceAbsDistance.Subavg = abs(valenceDiff_subavg);
            distances.arousalAbsDistance.Norm = abs(arousalDiff_norm);
            distances.arousalAbsDistance.Subspec = abs(arousalDiff_subspec);
            distances.arousalAbsDistance.Subavg = abs(arousalDiff_subavg);
            distances.valenceSignedDistance.Norm = valenceDiff_norm;
            distances.valenceSignedDistance.Subspec = valenceDiff_subspec;
            distances.valenceSignedDistance.Subavg = valenceDiff_subavg;
            distances.arousalSignedDistance.Norm = arousalDiff_norm;
            distances.arousalSignedDistance.Subspec = arousalDiff_subspec;
            distances.arousalSignedDistance.Subavg = arousalDiff_subavg;
            distance_file = fullfile(distance_dir, ['sub', subject], ['sub', subject, '_', modality, '_distances_run', num2str(run), '.mat']);
            if ~exist(fullfile(distance_dir, ['sub', subject]), 'dir')
                mkdir(fullfile(distance_dir, ['sub', subject]));
            end
            save(distance_file, 'distances');
            disp(['Saved distances to ', distance_file]);

            categories = struct();
            categories.start.Norm = startCategory_norm;
            categories.end.Norm = endCategory_norm;
            categories.start.Subspec = startCategory_subspec;
            categories.end.Subspec = endCategory_subspec;
            categories.start.Subavg = startCategory_subavg;
            categories.end.Subavg = endCategory_subavg;
            categories.startContinuous.Subspec = startCategoryContinuous_subspec;
            categories.endContinuous.Subspec = endCategoryContinuous_subspec;
            categories.startContinuous.Subavg = startCategoryContinuous_subavg;
            categories.endContinuous.Subavg = endCategoryContinuous_subavg;
            categories.startOnehot.Subspec = getMaxOnehotTable(startCategoryContinuous_subspec);
            categories.endOnehot.Subspec = getMaxOnehotTable(endCategoryContinuous_subspec);

            categories.startOnehot.Subavg = getMaxOnehotTable(startCategoryContinuous_subavg);
            categories.endOnehot.Subavg = getMaxOnehotTable(endCategoryContinuous_subavg);
            categories.startOnehot.Norm = startCategoryOnehot_norm;
            categories.endOnehot.Norm = endCategoryOnehot_norm;

            category_file = fullfile(category_dir, ['sub', subject], ['sub', subject, '_', modality, '_categories_run', num2str(run), '.mat']);
            if ~exist(fullfile(category_dir, ['sub', subject]), 'dir')
                mkdir(fullfile(category_dir, ['sub', subject]));
            end
            save(category_file, 'categories');
            disp(['Saved categories to ', category_file]);
        end
    end
end

disp('Processing completed for all subjects, modalities, and runs.');

function Max = getMaxOnehotTable(tableContinuous)
    
    data = table2array(tableContinuous);  
    [~, maxIdx] = max(data, [], 2);  %get index of max value per row
    
    onehotData = zeros(size(data));
    for i = 1:size(data, 1)
        onehotData(i, maxIdx(i)) = 1;
    end
    
    Max = array2table(onehotData, 'VariableNames', tableContinuous.Properties.VariableNames);
end

