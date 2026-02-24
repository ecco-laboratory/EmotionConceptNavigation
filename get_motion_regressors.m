clear; close all; clc;
[project_dir, fmri_data_dir, ~, ~, ~, ~, ~, ~, ~, motion_regressors_dir, ~, ~,~, ~, subjects]   = set_up_dirs_constants();

confounds_dir = fmri_data_dir;

modalities = {'word', 'face'};
runs = [1, 2]; % Two runs for each modality
TR = 1.25; % seconds
discard_time = 7.5; % Time to discard in seconds
discard_volumes = discard_time / TR; % Number of volumes to discard 

for s = 1:length(subjects)
    subject = subjects{s};

    for m = 1:length(modalities)
        modality = modalities{m};

        for r = 1:length(runs)
            run = runs(r);

            % Load motion regressors from fMRIprep
            confounds_file = fullfile(confounds_dir, ['sub-', subject], 'func', ...
            ['sub-', subject, '_task-', modality, '_run-0', num2str(run), '_desc-confounds_timeseries.tsv']);


            if exist(confounds_file, 'file')
                % extract the 6 motion parameters and their derivatives
                try
                    %skip if already exists
                    output_file = fullfile(motion_regressors_dir, ['sub', subject], ['sub', subject, '_motionRegressors_', modality, '_run', num2str(run), '.txt']);
                    if exist(output_file, 'file')
                        disp(['Skipping subject ', subject, ', modality ', modality, ', run ', num2str(run), ' because it is already processed']);
                        continue;
                    end
                    
                    if ~exist(fullfile(motion_regressors_dir, ['sub', subject]), 'dir')
                        mkdir(fullfile(motion_regressors_dir, ['sub', subject]));
                    end

                    motion_regressors = extract_motion_regressors(confounds_file, discard_volumes);
                    motion_table = array2table(motion_regressors);
                    writetable(motion_table, output_file, 'Delimiter', '\t', 'WriteVariableNames', false);
                    
                catch ME
                    error(['Error processing motion parameters for subject ', subject, ...
                        ', modality ', modality, ', run ', num2str(run), ': ', ME.message]);
                end
            else
                warning(['Motion regressor file not found: ', confounds_file]);
            end
        end
    end
end



% Function to extract motion parameters from fMRIprep confounds file
function motion_regressors = extract_motion_regressors(filename, discard_volumes)
    T = readtable(filename, 'FileType', 'text', 'Delimiter', '\t');
    
    motion_param_names = T.Properties.VariableNames(contains(T.Properties.VariableNames, 'trans') | ...
                                                    contains(T.Properties.VariableNames, 'rot'));

    motion_regressors = T{(discard_volumes+1):end, motion_param_names};
end