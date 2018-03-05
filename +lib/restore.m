function [dt_file, fa_file, md_file, v1_file, outliers_file] = restore(job_dir,dwmri_file,bvec_file,bval_file,mask_file,sigma,camino_exec)
    % Computes diffusion tensor and associated statistics with camino's
    % RESTORE algorithm.
    
    disp('---');
    disp('Running RESTORE analysis...');
    
    % Create "RESTORE" directory
    restore_dir = system_utils.directory(job_dir,'RESTORE');
    restore_dir.mkdir();
            
    % Create scheme file
    scheme_file = system_utils.file(restore_dir,'dwmri.scheme');
    system_utils.system_with_errorcheck([camino_exec.get_path('fsl2scheme') ' -bvalfile ' bval_file.get_path() ' -bvecfile ' bvec_file.get_path() ' -bscale 1 > ' scheme_file.get_path()],'failed to create scheme file');
    
    % Convert dwmri to bfloat
    dwmri_bfloat_file = system_utils.file(restore_dir,'dwmri.Bfloat');
    system_utils.system_with_errorcheck([camino_exec.get_path('image2voxel') ' -4dimage ' dwmri_file.get_path() ' -outputfile ' dwmri_bfloat_file.get_path()],'failed to convert dwmri to bfloat.');    
    
    % Call restore - use entire mask to compute diffusion tensor
    outliers_bin_file = system_utils.file(restore_dir,'outliers.bin'); % stored directly as bytes
    dt_bdouble_file = system_utils.file(restore_dir,'dt.Bdouble');
    system_utils.system_with_errorcheck([camino_exec.get_path('restore') ' ' dwmri_bfloat_file.get_path() ' ' scheme_file.get_path() ' ' num2str(sigma) ' ' outliers_bin_file.get_path() ' -bgmask ' mask_file.get_path() ' > ' dt_bdouble_file.get_path()],'Failed to run RESTORE.');
                
    % Calculate FA    
    fa_file = system_utils.file(restore_dir,'fa.nii.gz');
    system_utils.system_with_errorcheck(['cat ' dt_bdouble_file.get_path() ' | ' camino_exec.get_path('fa') ' | ' camino_exec.get_path('voxel2image') ' -outputroot ' fullfile(restore_dir.get_path(),'fa') ' -header ' dwmri_file.get_path()],'Failed to calculate FA.');
    
    % Calculate MD
    md_file = system_utils.file(restore_dir,'md.nii.gz');
    system_utils.system_with_errorcheck(['cat ' dt_bdouble_file.get_path() ' | ' camino_exec.get_path('md') ' | ' camino_exec.get_path('voxel2image') ' -outputroot ' fullfile(restore_dir.get_path(),'md') ' -header ' dwmri_file.get_path()],'Failed to calculate MD.');
            
    % Calculate eigenvectors/eigenvalues 
    dteig_bdouble_file = system_utils.file(restore_dir,'dteig.Bdouble');
    system_utils.system_with_errorcheck(['cat ' dt_bdouble_file.get_path() ' | ' camino_exec.get_path('dteig') ' > ' dteig_bdouble_file.get_path()],'Failed to calculate eigenvalues and eigenvectors.');
    
    % Load data - use mask_vol to get dimensions of data
    mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_file.get_path(),'logical');
    % dt stored as: [exitcode, ln(S(0)), D_xx, D_xy, D_xz, D_yy, D_yz, D_zz]    
    dt_bdouble_file.open('r','b');    
    dt_vol = double(permute(reshape(dt_bdouble_file.read('double'),[8 size(mask_vol,1) size(mask_vol,2) size(mask_vol,3)]),[2 3 4 1])); % Do reshape and permute to get it into suitable nifti format
    dt_bdouble_file.close();
    exitcode_vol = dt_vol(:,:,:,1);
    dt_vol = dt_vol(:,:,:,3:end);
    
    % Load fa and md
    fa_vol = nifti_utils.load_untouch_nii_vol_scaled(fa_file.get_path(),'double');
    md_vol = nifti_utils.load_untouch_nii_vol_scaled(md_file.get_path(),'double');
    
    % dt_eig stored as: [l1, e11, e12, e13, l2, e21, e22, e23, l3, e31, e32, e33]
    dteig_bdouble_file.open('r','b');    
    dteig_vol = double(permute(reshape(dteig_bdouble_file.read('double'),[12 size(mask_vol,1) size(mask_vol,2) size(mask_vol,3)]),[2 3 4 1])); % Do reshape and permute to get it into suitable nifti format
    dteig_bdouble_file.close();
        
    % Get v1_vol from dteig_vol
    v1_vol = dteig_vol(:,:,:,2:4);
    
    % Get outliers_vol    
    outliers_bin_file.open('r','b');    
    num_dwi = length(find(bval_file.dlmread() ~= 0)); 
    outliers_vol = double(permute(reshape(outliers_bin_file.read(),[num_dwi size(mask_vol,1) size(mask_vol,2) size(mask_vol,3)]),[2 3 4 1])); % Do reshape and permute to get it into suitable nifti format
    outliers_bin_file.close();
    
    % Remove camino-related files
    scheme_file.rm();
    dwmri_bfloat_file.rm();
    outliers_bin_file.rm();
    dt_bdouble_file.rm();
    dteig_bdouble_file.rm();
        
    % Save outputs as niftis
    template_nii = load_untouch_nii(mask_file.get_path());
    
    % exitcode
    exitcode_file = system_utils.file(restore_dir,'exitcode.nii.gz');
    template_nii.img = exitcode_vol;    
    nifti_utils.save_untouch_nii_using_scaled_img_info(exitcode_file.get_path(),template_nii,'double');
    
    % dt
    dt_file = system_utils.file(restore_dir,'dt.nii.gz');
    template_nii.img = dt_vol;    
    nifti_utils.save_untouch_nii_using_scaled_img_info(dt_file.get_path(),template_nii,'double');
    
    % fa
    fa_file = system_utils.file(restore_dir,'fa.nii.gz');
    template_nii.img = fa_vol;    
    nifti_utils.save_untouch_nii_using_scaled_img_info(fa_file.get_path(),template_nii,'double');
    
    % md
    md_file = system_utils.file(restore_dir,'md.nii.gz');
    template_nii.img = md_vol;    
    nifti_utils.save_untouch_nii_using_scaled_img_info(md_file.get_path(),template_nii,'double');
    
    % v1
    v1_file = system_utils.file(restore_dir,'v1.nii.gz');
    template_nii.img = v1_vol;    
    nifti_utils.save_untouch_nii_using_scaled_img_info(v1_file.get_path(),template_nii,'double');
  
    % outliers
    outliers_file = system_utils.file(restore_dir,'outliers.nii.gz');
    template_nii.img = outliers_vol;    
    nifti_utils.save_untouch_nii_using_scaled_img_info(outliers_file.get_path(),template_nii,'logical');
end