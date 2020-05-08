import os
from nipype.interfaces import gift

# GICA DEFAULTS
DEFAULT_DIM = 100
DEFAULT_ALG = 16
DEFAULT_ICA_PARAM_FILE = ''
DEFAULT_OUT_DIR = '.'
DEFAULT_DISPLAY_RESULTS = 1
DEFAULT_REFS = []
DEFAULT_RUN_NAME = 'COINSTAC_SCICA'
DEFAULT_GROUP_PCA_TYPE = 0
DEFAULT_BACK_RECON_TYPE = 5
DEFAULT_PREPROC_TYPE = 1
DEFAULT_NUM_REDUCTION_STEPS = 1
DEFAULT_SCALE_TYPE = 0
DEFAULT_GROUP_ICA_TYPE = 'spatial'
DEFAULT_WHICH_ANALYSIS = 1
DEFAULT_MASK = ''

matlab_cmd = '/computation/groupicatv4.0b/GroupICATv4.0b_standalone_sep_10_2019/run_groupica.sh /usr/local/MATLAB/MATLAB_Runtime/v91/'


def gift_gica(
    in_files=[], dim=DEFAULT_DIM, algoType=DEFAULT_ALG, refFiles=DEFAULT_REFS,
    run_name=DEFAULT_RUN_NAME, out_dir=DEFAULT_OUT_DIR, group_pca_type=DEFAULT_GROUP_PCA_TYPE,
    backReconType=DEFAULT_BACK_RECON_TYPE, preproc_type=DEFAULT_PREPROC_TYPE,
    numReductionSteps=DEFAULT_NUM_REDUCTION_STEPS, scaleType=DEFAULT_SCALE_TYPE,
    group_ica_type=DEFAULT_GROUP_ICA_TYPE, display_results=DEFAULT_DISPLAY_RESULTS,
    which_analysis=DEFAULT_WHICH_ANALYSIS, mask=DEFAULT_MASK
):
    """
    Wrapper for initializing GIFT nipype interface to run Group ICA.

    Args:
        in_files            (List [Str])    :   Input file names (either single file name or a list)
        dim                 (Int)           :   Dimensionality reduction into #num dimensions
        algoType            (Int)           :   options are 1 - Infomax, 2 - Fast ica , ...
        refFiles            (List [Str])    :   file names for reference templates (either single file name or a list)
        run_name            (Str)           :   Name of the analysis run
        out_dir             (Str)           :   Full file path of the results directory
        group_pca_type      (Str)           :   options are 'subject specific' and 'grand mean'
        backReconType       (Int)           :   options are 1 - regular, 2 - spatial-temporal regression, 3 - gica3, 4 - gica, 5 - gig-ica
        preproc_type        (Int)           :   options are 1 - remove mean per timepoint, 2 - remove mean per voxel, 3 - intensity norm, 4 - variance norm
        numReductionSteps   (Int)           :   Number of reduction steps used in the first pca step
        scaleType           (Int)           :   options are 0 - No scaling, 1 - percent signal change, 2 - Z-scores
        group_ica_type      (Str)           :   1 - Spatial ica, 2 - Temporal ica.
        display_results     (Int)           :   0 - No display, 1 - HTML report, 2 - PDF
        which_analysis      (Int)           :   Options are 1, 2, and 3. 1 - standard group ica, 2 - ICASSO and 3 - MST.
        mask                (Str)           :   Enter file names using full path of the mask. If you wish to use default mask leave it empty

        algoType full options:
        1           2           3       4           5       6
        'Infomax'   'Fast ICA'  'Erica' 'Simbec'    'Evd'   'Jade Opac',
        7           8           9                   10 
        'Amuse'     'SDD ICA'   'Semi-blind'        'Constrained ICA (Spatial)' 
        11              12      13          14      15          16          17
        'Radical ICA'   'Combi' 'ICA-EBM'   'ERBM'  'IVA-GL'    'GIG-ICA'   'IVA-L'

    Args (not supported here, but available for nipype):
        perfType            (Int)           :   Options are 1, 2, and 3. 1 - maximize performance, 2 - less memory usage  and 3 - user specified settings.
        prefix              (Str)           :   Enter prefix to be appended with the output files
        dummy_scans         (Int)           :   enter dummy scans
        numWorkers          (Int)           :   Number of parallel workers    
        doEstimation        (Int)           :   options are 0 and 1 


    """
    gift.GICACommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

    gc = gift.GICACommand()
    gc.inputs.in_files = in_files
    gc.inputs.algoType = algoType
    #gc.inputs.group_pca_type = group_pca_type
    gc.inputs.backReconType = backReconType
    gc.inputs.preproc_type = preproc_type
    gc.inputs.numReductionSteps = numReductionSteps
    gc.inputs.scaleType = scaleType
    gc.inputs.group_ica_type = group_ica_type
    gc.inputs.which_analysis = which_analysis
    gc.inputs.refFiles = refFiles
    gc.inputs.display_results = display_results
    gc.inputs.mask = mask

    if dim > 0:
        gc.inputs.dim = dim

    gc.inputs.out_dir = out_dir
    output = {}
    
    output = gc.run()
    
    
    return output
