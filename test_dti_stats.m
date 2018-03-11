%% Set environment
addpath('system_utils');
addpath(genpath('nifti_utils'));
addpath(genpath('dwmri_visualizer'));
addpath('dti_stats');

%% Run dti stats pipeline

% Set job directory path
job_dir_path = 'test_dti_stats';

% Use outputs from PREPROCESSED folder of topup_eddy_preprocess pipeline
dwmri_path = 'PREPROCESSED/dwmri.nii.gz';
bvec_path = 'PREPROCESSED/dwmri.bvec';
bval_path = 'PREPROCESSED/dwmri.bval';
mask_path = 'PREPROCESSED/mask.nii.gz';
movement_params_path = 'PREPROCESSED/eddy_params.txt';

% Set FSL path
fsl_path = '~/fsl_5_0_10_eddy_5_0_11';

% Set camino path
camino_path = '~/camino';

% BET params
bet_params = '-f 0.3 -R';

% CSF info
csf_info.label_path = 'dti_stats/csf_label/mni_icbm152_csf_tal_nlin_sym_09a_trunc.nii.gz';
csf_info.template_path = 'dti_stats/csf_label/mni_icbm152_t2_tal_nlin_sym_09a.nii.gz';
csf_info.template_masked_path = 'dti_stats/csf_label/mni_icbm152_t2_tal_nlin_sym_09a_mask.nii.gz';

% Perform dti stats pipeline
dti_stats_pdf_path = dti_stats(job_dir_path, ...
                               dwmri_path, ...
                               bvec_path, ...
                               bval_path, ...
                               mask_path, ...
                               fsl_path, ...
                               camino_path, ...
                               csf_info, ...
                               bet_params, ...
                               movement_params_path);
