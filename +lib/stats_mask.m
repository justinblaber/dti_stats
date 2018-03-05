function stats_mask_file = stats_mask(job_dir,dwmri_file,bval_file,mask_file,fsl_exec)
    % "stats mask" is a mask of the brain with CSF removed and is used
    % when computing the chi-squared table. Diffusion tensor fits in the
    % CSF are generally very poor.
 
    disp('---');
    disp('Forming stats mask...');
    
    % Create "STATS_MASK" directory
    stats_mask_dir = system_utils.directory(job_dir,'STATS_MASK');
    stats_mask_dir.mkdir();
                
    % Create averages of every dwi for each bval, then use FSL's FAST to 
    % segment CSF
    mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_file.get_path(),'logical');
    csf_vol = zeros(size(mask_vol)); % Initialize
    bvals = bval_file.dlmread();
    bvals_dwi_unique = unique(bvals(bvals ~= 0));
    for i = 1:length(bvals_dwi_unique)
        % Grab dwi
        dwi_file = system_utils.file(stats_mask_dir,['dwi' num2str(bvals_dwi_unique(i)) '.nii.gz']);
        nifti_utils.idx_untouch_nii4D(dwmri_file.get_path(),find(bval_file.dlmread() == bvals_dwi_unique(i)),dwi_file.get_path()); %#ok<FNDSB>
        
        % Get mean        
        dwi_mean_file = system_utils.file(stats_mask_dir,['dwi' num2str(bvals_dwi_unique(i)) '_mean.nii.gz']);
        nifti_utils.mean_untouch_nii4D(dwi_file.get_path(),dwi_mean_file.get_path());
        
        % Mask out region - use dwi_mean as a template
        dwi_mask_file = system_utils.file(stats_mask_dir,['dwi' num2str(bvals_dwi_unique(i)) '_mask.nii.gz']);
        dwi_mask_nii = load_untouch_nii(dwi_mean_file.get_path());
        dwi_mask_nii.img = nifti_utils.load_untouch_nii_vol_scaled(dwi_mean_file.get_path(),'double');
        dwi_mask_nii.img(~mask_vol) = 0;
        nifti_utils.save_untouch_nii_using_scaled_img_info(dwi_mask_file.get_path(),dwi_mask_nii,'double');
    
        % Run FSL's FAST
        fast_basename = fullfile(stats_mask_dir.get_path(),['fast_' num2str(bvals_dwi_unique(i))]);
        csf_file = system_utils.file(stats_mask_dir,['fast_' num2str(bvals_dwi_unique(i)) '_pve_0.nii.gz']); % Fast technically thinks the input is T1, so this is what corresponds to CSF
        system_utils.system_with_errorcheck([fsl_exec.get_path('fast') ' -o ' fast_basename ' -v ' dwi_mask_file.get_path()],'Failed to run FAST on averaged DWI');
        
        % Add to csf_vol
        csf_vol = csf_vol + nifti_utils.load_untouch_nii_vol_scaled(csf_file.get_path(),'double');
        
        % Remove temporary files                
        dwi_file.rm();
        dwi_mean_file.rm();
        dwi_mask_file.rm();
        system_utils.file(stats_mask_dir,['fast_' num2str(bvals_dwi_unique(i)) '_mixeltype.nii.gz']).rm();
        system_utils.file(stats_mask_dir,['fast_' num2str(bvals_dwi_unique(i)) '_pve_0.nii.gz']).rm();
        system_utils.file(stats_mask_dir,['fast_' num2str(bvals_dwi_unique(i)) '_pve_1.nii.gz']).rm();
        system_utils.file(stats_mask_dir,['fast_' num2str(bvals_dwi_unique(i)) '_pve_2.nii.gz']).rm();
        system_utils.file(stats_mask_dir,['fast_' num2str(bvals_dwi_unique(i)) '_pveseg.nii.gz']).rm();
        system_utils.file(stats_mask_dir,['fast_' num2str(bvals_dwi_unique(i)) '_seg.nii.gz']).rm();        
    end   
    
    % Get csf_vol - take average and threshold    
    csf_vol = csf_vol./length(bvals_dwi_unique) > 0.15; 
    
    % Save stats mask
    stats_mask_file = system_utils.file(stats_mask_dir,'stats_mask.nii.gz');
    stats_mask_nii = load_untouch_nii(mask_file.get_path());
    stats_mask_nii.img = imerode(mask_vol,ones(3)) & ~csf_vol; % Remove csf from from slightly eroded mask
    nifti_utils.save_untouch_nii_using_scaled_img_info(stats_mask_file.get_path(),stats_mask_nii,'logical');
end