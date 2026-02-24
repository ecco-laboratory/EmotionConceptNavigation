clear; close all; clc;
[project_dir, ~, psychopy_csv_dir, theta_dir, category_dir, vaDiff_dir, avgVA_dir, distance_dir, ~, ~, ~,~, ~, ~, subjects] = set_up_dirs_constants();

runs = {'01', '02'};
modalities = {'face', 'word'};
sources = {'subavg', 'subspec'};
beh_file = fullfile(project_dir, 'outputs', 'inoutscan_beh_allsubs_wMaxCategoryDistance.csv');
beh_data = readtable(beh_file);
if isnumeric(beh_data.sub_id), beh_data.sub_id = cellstr(num2str(beh_data.sub_id, '%04d'));end % 0001, 0002, etc.
if isnumeric(beh_data.run), beh_data.run = cellstr(num2str(beh_data.run, '%02d')); end% 01, 02, etc. end
output_dir = fullfile(project_dir, 'data', 'beh', 'inscan_judgment');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
for i = 1:length(subjects)
    subject = subjects{i};
    for m = 1:length(modalities)
        modality = modalities{m};
        for r = 1:length(runs)
            run = runs{r};
            disp(['Processing subject: ', subject, ', run: ', num2str(run), ', modality: ', modality]);

            %save outscan_consistency
            inscan_consistency_struct = struct();
            inscan_consistency_struct.Subspec = beh_data.outscan_consistency(strcmp(beh_data.sub_id, subject) & strcmp(beh_data.modality, modality) & strcmp(beh_data.run, run) & strcmp(beh_data.source, 'subspec'))';
            inscan_consistency_struct.Subavg = beh_data.outscan_consistency(strcmp(beh_data.sub_id, subject) & strcmp(beh_data.modality, modality) & strcmp(beh_data.run, run) & strcmp(beh_data.source, 'subavg'))';
            %make sure face run 01 or run0 02 have length 41 or 43,  word run 01 has length 62, word run 02 has length 62
            if strcmp(modality, 'face') && (strcmp(run, '01') || strcmp(run, '02'))
                if (length(inscan_consistency_struct.Subspec) ~= 41 && length(inscan_consistency_struct.Subspec) ~= 43) || (length(inscan_consistency_struct.Subavg) ~= 41 && length(inscan_consistency_struct.Subavg) ~= 43)
                    error('Face run 01 should have length 41 or 43 for subject %s', subject);
                end
            elseif strcmp(modality, 'word') && (strcmp(run, '01') || strcmp(run, '02'))
                if length(inscan_consistency_struct.Subspec) ~= 62 || length(inscan_consistency_struct.Subavg) ~= 62
                    error('Word run 01 or 02 should have length 62 for subject %s', subject);
                end
            end
            if ~exist(fullfile(output_dir, ['sub', subject]), 'dir')
                mkdir(fullfile(output_dir, ['sub', subject]));
            end
            save(fullfile(output_dir, ['sub', subject], ['sub', subject, '_', modality, '_consistency_run', num2str(str2num(run)), '.mat']), 'inscan_consistency_struct');
            disp(['Saved inscan judgment to ', fullfile(output_dir, ['sub', subject], ['sub', subject, '_', modality, '_consistency_run', num2str(str2num(run)), '.mat'])]);
            %save comparison_type as well
            inscan_comparison_struct = struct();
            inscan_comparison_struct.Subspec = beh_data.comparison_type(strcmp(beh_data.sub_id, subject) & strcmp(beh_data.modality, modality) & strcmp(beh_data.run, run) & strcmp(beh_data.source, 'subspec'))';
            inscan_comparison_struct.Subavg = beh_data.comparison_type(strcmp(beh_data.sub_id, subject) & strcmp(beh_data.modality, modality) & strcmp(beh_data.run, run) & strcmp(beh_data.source, 'subavg'))';
            %make sure Subspec and Subavg contain the same thing
            if ~isequal(inscan_comparison_struct.Subspec, inscan_comparison_struct.Subavg)
                error('Subspec and Subavg contain different things for subject %s', subject);
            end
            %make sure face run 01 has length 41, face run 02 has length 43, word run 01 has length 62, word run 02 has length 62   
            if strcmp(modality, 'face') && (strcmp(run, '01') || strcmp(run, '02'))
                if (length(inscan_comparison_struct.Subspec) ~= 41 && length(inscan_comparison_struct.Subspec) ~= 43) || (length(inscan_comparison_struct.Subavg) ~= 41 && length(inscan_comparison_struct.Subavg) ~= 43)
                   error('Face run 01 or 02 should have length 41 or 43 for subject %s', subject);
                end
            end
            if strcmp(modality, 'word') && (strcmp(run, '01') || strcmp(run, '02'))
                if length(inscan_comparison_struct.Subspec) ~= 62 || length(inscan_comparison_struct.Subavg) ~= 62
                    error('Word run 01 or 02 should have length 62 or 64 for subject %s', subject);
                end
            end
            save(fullfile(output_dir, ['sub', subject], ['sub', subject, '_', modality, '_comparison_run', num2str(str2num(run)), '.mat']), 'inscan_comparison_struct');
            disp(['Saved inscan comparison to ', fullfile(output_dir, ['sub', subject], ['sub', subject, '_', modality, '_comparison_run', num2str(str2num(run)), '.mat'])]);
        end
    end
end