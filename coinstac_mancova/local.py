import sys
import ujson as json
from .run_gift import gift_mancova, gift_gica
from utils import listRecursive
import utils as ut
import os
import glob
import pandas as pd
import scipy.io as sio

TC_SEARCH_STRING = 'gica_cmd_sub*%d_timecourses_ica_s1_.nii'


def convert_covariates(covariate_filename, state, covariate_types=None,  N=None):
    df = pd.read_csv(covariate_filename)
    out_dir = state["outputDirectory"]
    cov_types = {"name": [], "type": []}
    if covariate_types is not None:
        cov_types = pd.read_csv(covariate_types)
    covariates = {}
    for covariate_name in df.columns:
        covariate_series = df[covariate_name]
        if covariate_name == 'filename':
            continue
        if N is not None:
            covariate_series = covariate_series[:N]
        if covariate_name in cov_types['name']:
            cov_type = cov_types['types'][cov_types['name'].index(
                covariate_name)]
        else:
            cov_type = "continuous"
        fname = os.path.join(out_dir, "COINSTAC_COVAR_%s.txt" % covariate_name)
        with open(fname, 'w') as file:
            file.write('\n'.join([str(s) for s in list(covariate_series)]))
        ut.log("Wrote covariates %s to file %s" %
               (covariate_name, fname), state)
        covariates[covariate_name] = [cov_type, fname]
    return covariates


def parse_stats(filename):
    stats_loaded = sio.loadmat(filename, struct_as_record=False)
    mult_output = {}
    mult = stats_loaded['MULT'][0][0].__dict__

    stats = mult['stats'][0][0].__dict__
    mult_output['stats'] = stats
    mult_output['X'] = mult['X']
    mult_output['p'] = mult['p']
    mult_output['t'] = mult['t']

    uni_output = {}
    uni = stats_loaded['UNI'][0][0].__dict__
    stats = uni['stats'][0][0].__dict__
    uni_output['stats'] = stats
    uni_output['X'] = uni['X']
    uni_output['p'] = uni['p']
    uni_output['t'] = uni['t']

    return {'MULT': mult, 'UNI': uni}


def local_run_mancova(args):
    state = args["state"]
    ut.log("Got input %s" % (args["input"]), state)
    in_files = [os.path.join(state['baseDirectory'], f)
                for f in args["input"]["data"]]
    ut.log("Loaded files %s" % ', '.join(in_files), state)
    ext = '.csv'
    csv_filename = [i for i in args["input"]["data"] if ext in i]
    covariate_file = os.path.join(state["baseDirectory"], csv_filename[0])
    covariates = convert_covariates(covariate_file, state, covariate_types=None, N=len(in_files))

    ica_parameters = os.path.join(
        state['outputDirectory'], 'gica_cmd_ica_parameter_info.mat')
    maskfile = args["input"]["mask"]
    interactions = args["input"]["interactions"]
    ut.log("Interpolating", state)
    template = ut.get_interpolated_nifti(
        in_files[0], args["input"]["scica_template"], destination_dir=state['outputDirectory'])
    ut.log("Interpolated template at file %s" % template, state)
    pyscript = os.path.join(state["outputDirectory"], "pyscript_gicacommand.m")
    if os.path.exists(pyscript):
        os.remove(pyscript)
    ut.log("Running group ICA", state)
    output = gift_gica(
        in_files=in_files,
        refFiles=[template],
        mask=maskfile,
        out_dir=state["outputDirectory"],
        group_pca_type="subject specific",
        algoType=16,
    )
    subject_sms = list(glob.glob(os.path.join(
        state["outputDirectory"], 'gica_cmd_sub*_component_ica_s1_*.nii')))

    subject_tcs = []
    for i in range(1, (len(subject_sms)+1)):
        fn = os.path.join(state["outputDirectory"], TC_SEARCH_STRING) % i
        found = glob.glob(fn)
        if len(found) > 0:
            subject_tcs.append(found[0])
        else:
            break
    other_matfiles = [f for f in list(
        glob.glob(os.path.join(state['outputDirectory'], '*.mat')))]
    ut.log("Running Mancova", state)
    gift_mancova(
        ica_param_file=ica_parameters,
        out_dir=state["outputDirectory"],
        # run_name=DEFAULT_RUN_NAME,
        # comp_network_names=DEFAULT_COMP_NETWORK_NAMES,
        TR=2,
        features=["spatial maps", "timecourses spectra", "fnc correlations"],
        covariates=covariates,
        interactions=interactions,
        numOfPCs=[53],
        # feature_params=DEFAULT_FEATURE_PARAMS,
    )
    ut.log("Collecting Mancova results", state)
    FNC_FILE_1 = os.path.join(state['outputDirectory'],
                              'dfnc_stats',
                              'gica_cmd_mancovan_results_fnc.mat')
    fnc_output = {}
    if os.path.exists(FNC_FILE_1):
        fnc_output = parse_stats(FNC_FILE_1)
    FNC_FILE_2 = os.path.join(state['outputDirectory'],
                              'dfnc_stats',
                              'gica_cmd_mancovan_results_fnc_domain_avg.mat')
    fnc_output_avg = {}
    if os.path.exists(FNC_FILE_2):
        fnc_output_avg = parse_stats(FNC_FILE_2)
    SPECTRA_FILES = glob.glob(os.path.join(state['outputDirectory'],
                                           'spectra_stats',
                                           '*.mat'))
    spectra_results = {}
    for i, spectra_file in enumerate(SPECTRA_FILES):
        spectra_results[i] = parse_stats(spectra_file)

    sm_results = {}
    SM_FILES = glob.glob(os.path.join(state['outputDirectory'],
                                      'sm_stats',
                                      '*.mat'))
    for sm_file in SM_FILES:
        sm_results[i] = parse_stats(sm_file)
    output_dict = {
        'subject_sms': subject_sms,
        'subject_tcs': subject_tcs,
        'other_matfiles': other_matfiles,
        'fnc': fnc_output,
        'fnc_avg': fnc_output_avg,
        'spectra': spectra_results,
        'sm': sm_results,
        'computation_phase': 'scica_mancova_1'
    }
    cache_dict = {}
    computation_output = {
        "output": output_dict,
        "cache": cache_dict,
        "state": state
    }

    return computation_output


def scica_check_out(args):
    state = args["state"]
    gigica_dir = os.path.join(
        state["outputDirectory"], "gica_gmc_gica_results")
    output_dict = {
        "gigica_output": None
    }
    if os.path.exists(gigica_dir):
        output_dict["gigica_output"] = gigica_dir
    cache_dict = {}
    computation_output = {
        "output": output_dict,
        "cache": cache_dict,
        "state": state
    }
    return computation_output
