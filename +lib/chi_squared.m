function [chi_squared_plot, dwi_model_file] = chi_squared(job_dir,dwmri_file,bvec_file,bval_file,mask_file,dt_file,stats_mask_file)
    % Compute chi squared plot to indicate a goodness of fit of the
    % diffusion tensor. chi squared plot is done slice-wise where slices
    % are 3rd dimension of volume WRT data storage.

    disp('---');
    disp('Running chi squared analysis...');
    
    % Create "CHISQUARED" directory
    chi_squared_dir = system_utils.directory(job_dir,'CHISQUARED');
    chi_squared_dir.mkdir();
       
    % Load data
    dwmri_vol = nifti_utils.load_untouch_nii4D_vol_scaled(dwmri_file.get_path(),'double');
    bvecs = bvec_file.dlmread();
    bvals = bval_file.dlmread();
    mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_file.get_path(),'logical');
    dt_vol = nifti_utils.load_untouch_nii4D_vol_scaled(dt_file.get_path(),'double');
    stats_mask_vol = nifti_utils.load_untouch_nii_vol_scaled(stats_mask_file.get_path(),'logical');
       
    % Diffusion tensor stored as: [Dxx, Dxy, Dxz, Dyy, Dyz, Dzz]
    %
    % Signal reconstructed as:
    %   D = [Dxx Dxy Dxz;
    %        Dxy Dyy Dyz;
    %        Dzx Dyz Dzz]
    %   S = S0*exp(-b*g'*D*g) = S0*exp(-b*(Dxx*gxx^2 + 2*Dxy*gx*gy + 2*Dxz*gx*gz + Dyy*gy^2 + 2*Dyz*gy*gz + Dzz*gz^2))
    
    % Get relevent volumes to reconstruct estimated signal from diffusion
    % tensor
    b0_vol = nanmean(dwmri_vol(:,:,:,bvals == 0),4); % TODO: could possibly use geometric mean
    dwi_vol = dwmri_vol(:,:,:,bvals ~= 0);
    bvecs_dwi = bvecs(:,bvals ~= 0);
    bvals_dwi = bvals(bvals ~= 0);   
        
    % For simplicity, do one gradient direction at a time
    dwi_model_vol = zeros(size(dwi_vol));    
    for i = 1:size(dwi_vol,4)
        % Reshape to make getting S a vector equation. The exp() here can
        % cause Infs to appear if tensor fit is bad (i.e. in the CSF for
        % multi-shell dwi).
        dwi_model_vec = b0_vol(:).*exp(-bvals_dwi(i)*(...
            reshape(dt_vol(:,:,:,1),[],1)*bvecs_dwi(1,i)^2 + ...
            2*reshape(dt_vol(:,:,:,2),[],1)*bvecs_dwi(1,i)*bvecs_dwi(2,i) + ...
            2*reshape(dt_vol(:,:,:,3),[],1)*bvecs_dwi(1,i)*bvecs_dwi(3,i) + ...
            reshape(dt_vol(:,:,:,4),[],1)*bvecs_dwi(2,i)^2 + ...
            2*reshape(dt_vol(:,:,:,5),[],1)*bvecs_dwi(2,i)*bvecs_dwi(3,i) + ...
            reshape(dt_vol(:,:,:,6),[],1)*bvecs_dwi(3,i)^2));

        % Shape it back into a volume and store it
        dwi_model_vec = reshape(dwi_model_vec,size(dwi_vol,1),size(dwi_vol,2),size(dwi_vol,3));
        dwi_model_vec(~mask_vol) = 0;
        dwi_model_vol(:,:,:,i) = dwi_model_vec;
    end
    
    % Save dwi_model
    dwi_model_file = system_utils.file(chi_squared_dir,'dwi_model.nii.gz');
    dwi_model_nii = load_untouch_nii(mask_file.get_path());   
    dwi_model_nii.img = dwi_model_vol;
    nifti_utils.save_untouch_nii_using_scaled_img_info(dwi_model_file.get_path(),dwi_model_nii,'double');
        
    % chi squared plot - do plot with integration across slices, since
    % corruption usually affects entire slice (assuming dwmri are stored 
    % as collected from the scanner). Use stats_mask for this.
    chi_squared_plot = zeros(size(dwi_vol,3),size(dwi_vol,4));
    for i = 1:size(dwi_vol,3)
        for j = 1:size(dwi_vol,4)            
            % Get squared difference
            chi_squared_slice = squeeze(dwi_model_vol(:,:,i,j) - dwi_vol(:,:,i,j)).^2;
            chi_squared_slice(~stats_mask_vol(:,:,i)) = 0;
            chi_squared_slice(~isfinite(chi_squared_slice)) = 0;
                        
            % Normalize based on measured dwi slice "energy"
            dwi_slice = squeeze(dwi_vol(:,:,i,j)).^2;
            dwi_slice(~stats_mask_vol(:,:,i)) = 0;
            dwi_slice(~isfinite(dwi_slice)) = 0;
            
            % Sum and store
            chi_squared_plot(i,j) = sum(chi_squared_slice(:))./sum(dwi_slice(:));            
        end
    end    
    % Set anything non-finite to NaN.
    chi_squared_plot(~isfinite(chi_squared_plot)) = NaN;
    
    % Save chi squared plot
    chi_squared_file = system_utils.file(chi_squared_dir,'chi_squared_plot.txt');
    chi_squared_file.dlmwrite(chi_squared_plot,' ');
end