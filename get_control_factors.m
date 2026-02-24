clear; close all; clc;
[project_dir, ~, psychopy_csv_dir, theta_dir, category_dir, vaDiff_dir, avgVA_dir, distance_dir, ~, ~, ~,~, ~, ~, subjects] = set_up_dirs_constants();

id_setting_dir = fullfile(project_dir, 'data', 'beh', 'id_setting');


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
if ~exist(id_setting_dir, 'dir')
    mkdir(id_setting_dir);
end

for i = 1:length(subjects)
    subject = subjects{i};

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
            
            id_setting_norm = {}; id_setting_subspec = {}; id_setting_subavg = {};
            
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
                        id_setting_norm{end+1} = face_id;
                    else
                        id_setting_norm{end+1} = NaN;
                        warning(['No matching norm data for face_id: ', face_id, ', start_category: ', start_category, ', end_category: ', end_category]);
                    end
                    
                    %Match in subject-specific ratings
                    subspec_row = find(strcmp(subject_specific_ratings.group, face_id) & ...
                                      strcmp(subject_specific_ratings.ori_start_category, start_category) & ...
                                      strcmp(subject_specific_ratings.ori_end_category, end_category));
                    if ~isempty(subspec_row)
                        id_setting_subspec{end+1} = face_id;
                    else
                        id_setting_subspec{end+1} = NaN;
                        warning(['No matching subject-specific data for face_id: ', face_id, ', start_category: ', start_category, ', end_category: ', end_category]);
                    end

                    %Match in subject-averaged ratings
                    subavg_row = find(strcmp(subject_averaged_ratings.group, face_id) & ...
                                     strcmp(subject_averaged_ratings.ori_start_category, start_category) & ...
                                     strcmp(subject_averaged_ratings.ori_end_category, end_category));
                    if ~isempty(subavg_row)
                        id_setting_subavg{end+1} = face_id;
                    else
                        id_setting_subavg{end+1} = NaN;
                        warning(['No matching subject-averaged data for face_id: ', face_id, ', start_category: ', start_category, ', end_category: ', end_category]);
                    end
                end

                
                if ~isempty(id_setting_norm)
                   %check length of id_setting_norm, id_setting_subspec, id_setting_subavg
                   if length(scanner_beh_file.group(2:end-1)) ~= length(id_setting_norm) || length(scanner_beh_file.group(2:end-1)) ~= length(id_setting_subspec) || length(scanner_beh_file.group(2:end-1)) ~= length(id_setting_subavg)
                       warning(['Length of scanner_beh_file(2:end-1) is not the same as the length of id_setting_norm, id_setting_subspec, id_setting_subavg for subject ', subject]);
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
                        id_setting_norm{end+1} = group;
                    else
                        id_setting_norm{end+1} = NaN;
                        warning(['No matching norm data for group: ', group, ', start_word: ', start_word, ', end_word: ', end_word]);
                    end
                    
                    % Match in subject-specific ratings
                    subspec_row = find(strcmp(subject_specific_ratings.group, group) & ...
                                      strcmp(subject_specific_ratings.start_stimulus, start_word) & ...
                                      strcmp(subject_specific_ratings.end_stimulus, end_word));
                    if ~isempty(subspec_row)
                        id_setting_subspec{end+1} = group;
                    else
                        id_setting_subspec{end+1} = NaN;
                        warning(['No matching subject-specific data for group: ', group, ', start_word: ', start_word, ', end_word: ', end_word]);
                    end
                    
                    % Match in subject-averaged ratings
                    subavg_row = find(strcmp(subject_averaged_ratings.group, group) & ...
                                     strcmp(subject_averaged_ratings.start_stimulus, start_word) & ...
                                     strcmp(subject_averaged_ratings.end_stimulus, end_word));
                    if ~isempty(subavg_row)
                        id_setting_subavg{end+1} = group;
                    else
                        id_setting_subavg{end+1} = NaN;
                        warning(['No matching subject-averaged data for group: ', group, ', start_word: ', start_word, ', end_word: ', end_word]);
                    end
                end
                

                if ~isempty(id_setting_norm)
                   %check length of id_setting_norm, id_setting_subspec, id_setting_subavg
                   if length(scanner_beh_file.group(2:end-1)) ~= length(id_setting_norm) || length(scanner_beh_file.group(2:end-1)) ~= length(id_setting_subspec) || length(scanner_beh_file.group(2:end-1)) ~= length(id_setting_subavg)
                       warning(['Length of scanner_beh_file(2:end-1) is not the same as the length of id_setting_norm, id_setting_subspec, id_setting_subavg for subject ', subject]);
                       continue;
                   end
                end
            end
            

            id_setting = struct();
            id_setting.Norm = id_setting_norm;
            id_setting.Subspec = id_setting_subspec;
            id_setting.Subavg = id_setting_subavg;
            id_setting_file = fullfile(id_setting_dir, ['sub', subject], ['sub', subject, '_', modality, '_idsetting_run', num2str(run), '.mat']);
            if ~exist(fullfile(id_setting_dir, ['sub', subject]), 'dir')
                mkdir(fullfile(id_setting_dir, ['sub', subject]));
            end
            save(id_setting_file, 'id_setting');
            disp(['Saved id_setting to ', id_setting_file]);
        end
    end
end

disp('Processing completed for all subjects, modalities, and runs.');


