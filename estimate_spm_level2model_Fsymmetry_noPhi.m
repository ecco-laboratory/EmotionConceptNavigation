[project_dir, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, spm1stlevel_dir, spm2ndlevel_dir, ~, subjects] = set_up_dirs_constants();

spm1stlevel_dir = [spm1stlevel_dir, '_noICA_wButton'];
spm2ndlevel_dir = [spm2ndlevel_dir, '_noICA_wButton'];

%periodicity = {6, 4, 5, 7, 8};
contrasts_to_analyze = {
    'ftest_CosSinTheta', ...
    'face_ftest_SinCosTheta', ...
    'word_ftest_SinCosTheta'
};
theta_sources = {'Subavg', 'Subspec'};


base_output_dir = fullfile(spm2ndlevel_dir, 'symmetry_noPhi_Z_spm', 'includeNonNorm');
if ~exist(base_output_dir, 'dir'), mkdir(base_output_dir); end
template_nifti_object = fmri_data(fullfile(spm1stlevel_dir, ['sub', subjects{1}], 'singleTrial', 'beta_0001.nii'));
save_zmaps = 0;

for c = 1:length(contrasts_to_analyze)
    current_contrast_name = contrasts_to_analyze{c};
    
    disp(['Setting up 2nd level for: Contrast ', current_contrast_name]);

    for p = 1:length(periodicity)
        current_periodicity = periodicity{p};

        for t = 1:length(theta_sources)
            current_theta_source = theta_sources{t};

            if save_zmaps

                Fimg_files_allsubs = {}; ddf_allsubs = zeros(1, length(subjects));
            
                for s_idx = 1:length(subjects)
                    subject = subjects{s_idx};
                    
                    subj_1st_level_output_dir = fullfile(spm1stlevel_dir, ['sub', subject], 'symmetry_noPhi', 'includeNonNorm', ['periodicity', num2str(current_periodicity)], ['theta', current_theta_source]);
                    
                    if ~exist(subj_1st_level_output_dir, 'dir')
                        warning(['1st-level output directory not found for subject ', subject, ' at: ', subj_1st_level_output_dir, '. Skipping subject.']);
                        continue; 
                    end
                    %load SPM.mat to find the correct contrast index
                    spm_mat_path = fullfile(subj_1st_level_output_dir, 'SPM.mat');
                    if exist(spm_mat_path, 'file')
                        SPM_load = load(spm_mat_path, 'SPM');
                        SPM = SPM_load.SPM;

                        contrast_idx = -1; 
                        
                        for k = 1:length(SPM.xCon)
                            if strcmp(SPM.xCon(k).name, current_contrast_name) 
                                contrast_idx = k;
                                break;
                            end
                        end
                        
                        if contrast_idx ~= -1
                            Fmap_image_name = sprintf('spmF_%04d.nii', contrast_idx);
                            Fmap_image_path = fullfile(subj_1st_level_output_dir, Fmap_image_name);

                            if exist(Fmap_image_path, 'file')
                                Fimg_files_allsubs{end+1} = Fmap_image_path;
                                ddf_allsubs(s_idx) = SPM.xX.erdf;
                            else
                                warning(['Contrast image (', contrast_image_name, ') not found for subject ', subject, ' in ', contrast_image_path]);
                            end
                        else
                            warning(['Contrast "', current_contrast_name, '" not found in SPM.xCon for subject ', subject, ' in ', spm_mat_path]);
                        end


                    else
                        warning(['SPM.mat not found for subject ', subject, ' in ', subj_1st_level_output_dir]);
                    end
                end
                Fstats_data = fmri_data(Fimg_files_allsubs).dat; %voxel*subjects
                disp(size(Fstats_data))
                %convert F map to z-statistic map
                ndf = 2; %cos and sin are two independent effects tested
                %Zstats_data = f2z(Fstats_data, ndf, ddf_allsubs);
                for s = 1:size(Fstats_data, 2)
                    Fvals = Fstats_data(:, s);

                    Z = nan(size(Fvals));
                    ok = isfinite(Fvals) & Fvals > 0;
                    Z(ok) = f2z(Fvals(ok), ndf, ddf_allsubs(s));
                    %save Z map as nifti file for this subject
                    template_nifti_object.dat = Z;
                    nifti_dir = fullfile(base_output_dir, 'Zmaps', ['periodicity', num2str(current_periodicity)], ['theta', current_theta_source]);
                    if ~exist(nifti_dir, 'dir'), mkdir(nifti_dir); end
                    nifti_file_path = fullfile(nifti_dir, [current_contrast_name, '_Z_sub', subjects{s}, '.nii']);
                    template_nifti_object.fullpath = nifti_file_path;
                    template_nifti_object.write;
                end
            end


            output_dir_2nd_level = fullfile(base_output_dir, 'SPM', ['periodicity', num2str(current_periodicity)], ['theta', current_theta_source], current_contrast_name);
            if ~exist(output_dir_2nd_level, 'dir'), mkdir(output_dir_2nd_level); end

            matlabbatch = {};

            matlabbatch{1}.spm.stats.factorial_design.dir = {output_dir_2nd_level};
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = {}; 
            for s_idx = 1:length(subjects)
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{end+1, 1} = fullfile(base_output_dir, 'Zmaps', ['periodicity', num2str(current_periodicity)], ['theta', current_theta_source], [current_contrast_name, '_Z_sub', subjects{s_idx}, '.nii']);
            end

            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

            % Model Estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat = {[output_dir_2nd_level, filesep, 'SPM.mat']};
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

            % Contrast specification 
            matlabbatch{3}.spm.stats.con.spmmat = {[output_dir_2nd_level, filesep, 'SPM.mat']};
            matlabbatch{3}.spm.stats.con.delete = 0;

            matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = [current_contrast_name, '_GroupTtest'];
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1; % one-sample t-test on a con image
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

            try
                spm_jobman('run', matlabbatch);
                disp(['Successfully completed 2nd-level analysis for: Contrast ', current_contrast_name]);
            catch ME
                warning(['Error during 2nd-level analysis for contrast ', current_contrast_name, ': ', ME.message]);
            end
            
        end 
    end
end