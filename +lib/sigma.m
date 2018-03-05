function sigma_est = sigma(job_dir,dwmri_file,bval_file,stats_mask_file,csf_info,bet_params,fsl_exec)
    % Computes standard deviation in CSF as an input to RESTORE
    
    disp('---');
    disp('Calculating sigma (standard deviation)...');
    
    % Create "SIGMA" directory
    sigma_dir = system_utils.directory(job_dir,'SIGMA');
    sigma_dir.mkdir();
    
    % Get average B0 as target
    b0_all_file = system_utils.file(sigma_dir,'b0_all.nii.gz');
    b0_mean_file = system_utils.file(sigma_dir,'b0_mean.nii.gz');
    b0_mean_mask_file = system_utils.file(sigma_dir,'b0_mean_mask.nii.gz');
    nifti_utils.idx_untouch_nii4D(dwmri_file.get_path(),find(bval_file.dlmread() == 0),b0_all_file.get_path()); %#ok<FNDSB>
    nifti_utils.mean_untouch_nii4D(b0_all_file.get_path(),b0_mean_file.get_path());
    system_utils.system_with_errorcheck([fsl_exec.get_path('bet') ' ' b0_mean_file.get_path() ' ' b0_mean_mask_file.get_path() ' ' bet_params],'Failed to generate mask with BET');
        
    % First do 12 DOF registration of masked T2 template to masked averaged B0.
    % Note that both inputs are skull stripped.
    T2_b0_xform_flirt_file = system_utils.file(sigma_dir,'T2_b0_xform_flirt.mat');
    T2_b0_flirt_file = system_utils.file(sigma_dir,'T2_b0_flirt.nii.gz');
    system_utils.system_with_errorcheck([fsl_exec.get_path('flirt') ' -cost normmi -searchcost normmi -dof 12 -omat ' T2_b0_xform_flirt_file.get_path() ' -in  ' csf_info.template_masked_path ' -ref ' b0_mean_mask_file.get_path() ' -out ' T2_b0_flirt_file.get_path() ' -interp sinc'],'Failed to register template to averaged B0.');
        
    % Use 12 DOF as initial guess to fnirt. Note that neither inputs are
    % skull stripped.
    T2_b0_xform_fnirt_file = system_utils.file(sigma_dir,'T2_b0_xform_fnirt.nii.gz');
    T2_b0_fnirt_file = system_utils.file(sigma_dir,'T2_b0_fnirt.nii.gz');
    system_utils.system_with_errorcheck([fsl_exec.get_path('fnirt') ' --aff=' T2_b0_xform_flirt_file.get_path() ' --cout=' T2_b0_xform_fnirt_file.get_path() ' --in=' csf_info.template_path ' --ref=' b0_mean_file.get_path() ' --iout=' T2_b0_fnirt_file.get_path()],'Failed to apply fnirt.');

    % Apply transformation to CSF label
    CSF_label_B0_file = system_utils.file(sigma_dir,'CSF_label_b0.nii.gz');
    system_utils.system_with_errorcheck([fsl_exec.get_path('applywarp') ' --ref=' T2_b0_fnirt_file.get_path() ' --in=' csf_info.label_path ' --warp=' T2_b0_xform_fnirt_file.get_path() ' --out=' CSF_label_B0_file.get_path() ' --interp=trilinear'],'Failed to apply transformation to labels.');
    
    % Get standard deviation in DWI within CSF template - this is the sigma
    % of the BACKGROUND (or close to it) noise
    bvals = bval_file.dlmread();
    dwmri_vols = nifti_utils.load_untouch_nii4D_vol_scaled(dwmri_file.get_path(),'double');    
    % Threshold csf_vol since it's a probability atlas
    csf_vol = nifti_utils.load_untouch_nii_vol_scaled(CSF_label_B0_file.get_path(),'double') > 0.1;
    
    % Figure out what to do to csf volume
    if length(find(csf_vol)) > 20
        % Decent number of csf voxels - attempt to union csf results from 
        % atlas and stats mask (FAST) to refine results
        csf_union_vol = csf_vol & ~nifti_utils.load_untouch_nii_vol_scaled(stats_mask_file.get_path(),'logical');
        if length(find(csf_union_vol)) > 20
            % Set csf_vol as the union
            csf_vol = csf_union_vol;
        else
            warning(['csf atlas and FAST disagree on CSF location... ' ...
                     'sigma computation might be poor. Using atlas ' ...
                     'for csf location.']);
        end
    elseif ~any(csf_vol(:))
        % No csf voxels
        warning(['No CSF voxels were found after registration of csf ' ...
                 'atlas. This probably means registration failed. ' ...
                 'Using the whole brain for csf mask which will ' ...
                 'overestimate sigma but will allow pipeline to complete.']);
        csf_vol = nifti_utils.load_untouch_nii_vol_scaled(b0_mean_mask_file.get_path(),'logical');
    else
        % Very few csf voxels
        warning(['Very few voxels in csf atlas after registration. ' ...
                 'Sigma computation might be poor. Using atlas for csf ' ...
                 'location.']);
    end
    
    sigma_dwi = [];
    for i = 1:length(bvals)
        if bvals(i) ~= 0
           % Use MAD to estimate sigma
           dwi_vol = dwmri_vols(:,:,:,i);
           MAD = nanmedian(abs(dwi_vol(csf_vol) - nanmedian(dwi_vol(csf_vol))));
           sigma_dwi = [sigma_dwi 1.4826*MAD]; %#ok<AGROW>
        end
    end
       
    % Get estimated sigma
    sigma_est = nanmedian(sigma_dwi);    
    disp(['Background sigma estimated to be: ' num2str(sigma_est)]);

    % Write sigma file sigma value
    sigma_file = system_utils.file(sigma_dir,'sigma_est.txt');
    sigma_file.dlmwrite(sigma_est,' ');
    
    % Clear out other files
    b0_all_file.rm();
    b0_mean_file.rm();
    b0_mean_mask_file.rm();
    T2_b0_xform_flirt_file.rm();
    T2_b0_flirt_file.rm();
    T2_b0_xform_fnirt_file.rm();
    T2_b0_fnirt_file.rm();
    CSF_label_B0_file.rm();
end