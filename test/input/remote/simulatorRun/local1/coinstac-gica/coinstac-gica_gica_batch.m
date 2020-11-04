%% Batch script for running gica

%% Performance type
perfType = 1;
%% Reliability analysis
which_analysis = 1;
%% Output directory
outputDir = '/transfer/local1/simulatorRun/coinstac-gica';

%% Output files prefix
prefix = 'coinstac-gica';

dataSelectionMethod = 4;

%% Input file patterns
input_data_file_patterns = {'/input/local1/simulatorRun/rtswa_1FA_NL11.nii';
'/input/local1/simulatorRun/rtswa_2WJ_E511.nii';
'/input/local1/simulatorRun/rtswa_255_NL11.nii';
'/input/local1/simulatorRun/rtswa_3G7_F111.nii';
'/input/local1/simulatorRun/rtswa_11G_1G11.nii';
'/input/local1/simulatorRun/rtswa_26F_BA11.nii';
'/input/local1/simulatorRun/rtswa_1TN_F111.nii';
'/input/local1/simulatorRun/rtswa_2TD_MB11.nii';
'/input/local1/simulatorRun/rtswa_18S_NL11.nii';
'/input/local1/simulatorRun/rtswa_1CJ_NL11.nii';
'/input/local1/simulatorRun/rtswa_156_NL11.nii';
'/input/local1/simulatorRun/rtswa_2P7_2X31.nii';
'/input/local1/simulatorRun/rtswa_2A5_NL11.nii';
'/input/local1/simulatorRun/rtswa_1BW_1G11.nii';
'/input/local1/simulatorRun/rtswa_1GG_NX11.nii';
'/input/local1/simulatorRun/rtswa_1NF_BN11.nii';
'/input/local1/simulatorRun/rtswa_2HA_F111.nii';
'/input/local1/simulatorRun/rtswa_1D6_4311.nii';
'/input/local1/simulatorRun/rtswa_1RB_B211.nii';
'/input/local1/simulatorRun/rtswa_1PH_5811.nii';
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
refFiles = {'/transfer/local1/simulatorRun/coinstac-gica/NeuroMark_INTERP_17_20_17.nii';
};

%% Report generator 
display_results = 1;
