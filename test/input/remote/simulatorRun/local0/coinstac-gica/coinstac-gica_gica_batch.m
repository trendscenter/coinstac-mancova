%% Batch script for running gica

%% Performance type
perfType = 1;
%% Reliability analysis
which_analysis = 1;
%% Output directory
outputDir = '/transfer/local0/simulatorRun/coinstac-gica';

%% Output files prefix
prefix = 'coinstac-gica';

dataSelectionMethod = 4;

%% Input file patterns
input_data_file_patterns = {'/input/local0/simulatorRun/rtswa_1TR_M511.nii';
'/input/local0/simulatorRun/rtswa_1YM_2X11.nii';
'/input/local0/simulatorRun/rtswa_34F_YG11.nii';
'/input/local0/simulatorRun/rtswa_1DH_NL11.nii';
'/input/local0/simulatorRun/rtswa_1MD_F111.nii';
'/input/local0/simulatorRun/rtswa_2CF_YT11.nii';
'/input/local0/simulatorRun/rtswa_14R_3F11.nii';
'/input/local0/simulatorRun/rtswa_1SF_MB11.nii';
'/input/local0/simulatorRun/rtswa_2RG_F111.nii';
'/input/local0/simulatorRun/rtswa_1NB_NL11.nii';
'/input/local0/simulatorRun/rtswa_2LV_YG11.nii';
'/input/local0/simulatorRun/rtswa_2EX_YG11.nii';
'/input/local0/simulatorRun/rtswa_3FW_6J11.nii';
'/input/local0/simulatorRun/rtswa_33N_7P11.nii';
'/input/local0/simulatorRun/rtswa_1JR_S911.nii';
'/input/local0/simulatorRun/rtswa_38R_F111.nii';
'/input/local0/simulatorRun/rtswa_1VX_F111.nii';
'/input/local0/simulatorRun/rtswa_2ZG_F111.nii';
'/input/local0/simulatorRun/rtswa_3JG_NL11.nii';
'/input/local0/simulatorRun/rtswa_383_F111.nii';
};

%% Dummy scans
dummy_scans = 0;
%% Input mask
maskFile = '/computation/local_data/mask.nii';

%% Group PCA type 
group_pca_type = 'subject specific';
%% PCA Algorithm
pcaType = 'Standard';
%% ICA Algorithm
algoType = 16;
%% Back-reconstruction type
backReconType = 1;
%% Pre-processing type
preproc_type = 1;
%% Number of data reduction steps
numReductionSteps = 1;
%% MDL Estimation 
doEstimation = 0;
%% Number of PC in the first PCA step
numOfPC1 = 53;
%% Scaling type 
scaleType = 1;
%% Spatial references 
refFiles = {'/transfer/local0/simulatorRun/coinstac-gica/NeuroMark_INTERP_17_20_17.nii';
};

%% Report generator 
display_results = 1;
