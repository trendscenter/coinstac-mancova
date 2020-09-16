import sys
import ujson as json

# import json
from .run_gift import gift_mancova, gift_gica
from utils import listRecursive
import utils as ut
import os
import glob
import pandas as pd
import scipy.io as sio
import shutil
import h5py
import numpy as np

TC_SEARCH_STRING = "gica_cmd_sub*%d_timecourses_ica_s1_.nii"

NEUROMARK_NETWORKS = {
    "SC": [1, 2, 3, 4, 5],
    "AUD": [6, 7],
    "SM": [8, 9, 10, 11, 12, 13, 14, 15, 16],
    "VIS": [17, 18, 19, 20, 21, 22, 23, 24, 25],
    "CC": [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
    "DMN": [43, 44, 45, 46, 47, 48, 49],
    "CR": [50, 51, 52, 53],
}


def convert_covariates(covariate_filename, state, covariate_types=None, N=None):
    df = pd.read_csv(covariate_filename)
    out_dir = state["outputDirectory"]
    cov_types = {"name": [], "type": []}
    if covariate_types is not None:
        cov_types = pd.read_csv(covariate_types)
    covariates = {}
    ut.log("Covariate Types %s" % (cov_types), state)
    ut.log("Covariate Types name %s" % (cov_types["name"].values), state)
    for covariate_name in df.columns:
        if covariate_types is not None:
            if covariate_name not in cov_types["name"].values:
                continue
        covariate_series = df[covariate_name]
        if covariate_name == "filename":
            continue
        if N is not None:
            covariate_series = covariate_series[:N]
        if covariate_name in cov_types["name"]:
            cov_type = cov_types["types"][cov_types["name"].index(covariate_name)]
        else:
            cov_type = "continuous"
        fname = os.path.join(out_dir, "COINSTAC_COVAR_%s.txt" % covariate_name)
        with open(fname, "w") as file:
            file.write("\n".join([str(s) for s in list(covariate_series)]))
        ut.log("Wrote covariates %s to file %s" % (covariate_name, fname), state)
        covariates[covariate_name] = [cov_type, fname]
    ut.log("Covariates are %s" % (covariates), state)
    return covariates, df


def loadmat(filename):
    """
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """
    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def _check_keys(dict):
    """
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    """
    for key in dict:
        if isinstance(dict[key], sio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict


def _todict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries
    """
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, sio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict


def parse_stats(filename):
    using_h5py = False
    try:
        stats_loaded = loadmat(filename)
    except Exception:
        # f = h5py.
        using_h5py = True
        stats_loaded = h5py.File(filename, "r")
        # stats_loaded = _check_keys(stats_loaded)
    mult_output = {}
    if not using_h5py:
        mult = stats_loaded["MULT"][0][0].__dict__
    else:
        mult = stats_loaded["MULT"]
    if not using_h5py:
        stats = mult["stats"][0][0].__dict__
    else:
        stats = {k: v[()] for k, v in mult["stats"].items()}
    mult_output["stats"] = stats
    mult_output["X"] = mult["X"][()]
    mult_output["p"] = mult["p"][()]
    mult_output["t"] = mult["t"][()]

    uni_output = {}
    if not using_h5py:
        uni = stats_loaded["UNI"][0][0].__dict__
        stats = uni["stats"][0][0].__dict__
    else:
        uni = stats_loaded["UNI"]
        stats = uni["stats"][()]  # {k:v[()] for k,v in uni['stats'].items()}
    # uni_output['stats'] = stats
    # uni_output['X'] = uni['X']
    uni_output["p"] = uni["p"][()]
    uni_output["t"] = uni["t"][()]

    return {"MULT": mult_output, "UNI": uni_output}


def local_run_mancova(args):
    state = args["state"]
    ut.log("Got input %s" % (args["input"]), state)
    csv_filename = [i for i in args["input"]["data"] if ".csv" in i]
    univariate_test_list = args["input"]["univariate_test_list"]
    covariate_file = os.path.join(state["baseDirectory"], csv_filename[0])
    covariate_type_file = os.path.join(state["baseDirectory"], csv_filename[1])
    ut.log("Covariate File Name:" + covariate_file, state)
    file_list = args["input"]["data"]
    file_list.remove(csv_filename[0])
    file_list.remove(csv_filename[1])

    in_files = [os.path.join(state["baseDirectory"], f) for f in file_list]
    ut.log("Loaded files %s" % ", ".join(in_files), state)
    covariates, covariates_df = convert_covariates(
        covariate_file, state, covariate_types=covariate_type_file, N=len(in_files)
    )

    maskfile = args["input"]["mask"]

    pyscript = os.path.join(state["outputDirectory"], "pyscript_gicacommand.m")
    if os.path.exists(pyscript):
        os.remove(pyscript)
    if not os.path.exists(
        os.path.join(state["baseDirectory"], args["input"]["gica_input_dir"])
    ):
        ut.log(
            "The input directory %s, does not exist"
            % (
                str(
                    os.path.join(
                        state["baseDirectory"], args["input"]["gica_input_dir"]
                    )
                )
            ),
            state,
        )
        args["input"]["gica_input_dir"] = None
    if args["input"]["gica_input_dir"] is None:
        ut.log("Interpolating", state)
        template = ut.get_interpolated_nifti(
            in_files[0],
            args["input"]["scica_template"],
            destination_dir=state["outputDirectory"],
        )
        ut.log("Interpolated template at file %s" % template, state)
        ut.log("Running group ICA", state)
        gica_out_dir = os.path.join(state["outputDirectory"], "coinstac-gica")
        os.makedirs(gica_out_dir, exist_ok=True)
        output = gift_gica(
            in_files=in_files,
            refFiles=template,
            mask=maskfile,
            out_dir=gica_out_dir,
            group_pca_type="subject specific",
            algoType=16,
            run_name="coinstac-gica",
        )
    else:
        ut.log(
            "Copying preexisting output from %s to %s"
            % (
                os.path.join(state["baseDirectory"], args["input"]["gica_input_dir"]),
                os.path.join(state["outputDirectory"], args["input"]["gica_input_dir"]),
            ),
            state,
        )
        for filename in glob.glob(
            os.path.join(state["baseDirectory"], args["input"]["gica_input_dir"])
        ):
            if os.path.isdir(filename):
                shutil.copytree(
                    filename,
                    os.path.join(state["outputDirectory"], os.path.basename(filename),),
                )
            else:
                shutil.copy(
                    filename,
                    os.path.join(
                        state["outputDirectory"],
                        args["input"]["gica_input_dir"],
                        os.path.basename(filename),
                    ),
                )
    ica_parameter_list = list(
        glob.iglob(
            os.path.join(state["outputDirectory"], "**", "*parameter_info.mat"),
            recursive=True,
        )
    )
    if len(ica_parameter_list) == 0:
        raise (ValueError("Ica parameter file could not be found in the outupt"))
    ica_parameters = ica_parameter_list[0]
    ut.log("ICA Parameters file is %s" % (str(ica_parameters)), state)
    ut.log("Checking covariates %s" % (str(covariates)), state)
    maskfile = args["input"]["mask"]
    interactions = args["input"]["interactions"]
    if args["input"]["run_mancova"]:
        ut.log("Running Mancova", state)
        try:
            multivariate_out_dir = os.path.join(
                state["outputDirectory"], "coinstac-multivariate"
            )
            os.makedirs(multivariate_out_dir, exist_ok=True)
            gift_mancova(
                ica_param_file=ica_parameters,
                out_dir=multivariate_out_dir,
                TR=args["input"].get("TR", 2),
                features=args["input"]["features"],
                covariates=covariates,
                interactions=interactions,
                numOfPCs=args["input"].get("numOfPCs", [4, 4, 4]),
                run_name="coinstac-mancovan-multivariate",
                comp_network_names=args["input"].get(
                    "comp_network_names", NEUROMARK_NETWORKS
                ),
                freq_limits=args["input"].get("freq_limits", [0.1, 0.15]),
                t_threshold=args["input"].get("t_threshold", 0.05),
                image_values=args["input"].get("image_values", "positive"),
                threshdesc=args["input"].get("threshdesc", "fdr"),
                p_threshold=args["input"].get("p_threshold", 0.05),
                display_p_threshold=args["input"].get("display_p_threshold", 0.05),
            )
        except Exception as e:
            ut.log(
                "Multivariate analysis raised an exception, likely because of bad conditioning. More subjects are required. Full error string {err}".format(
                    err=str(e)
                ),
                state,
            )
    else:
        ut.log("Skipping Mancova", state)
    if args["input"]["run_univariate_tests"]:
        ut.log("Running univariate tests", state)
        for univariate_test in univariate_test_list:
            key = list(univariate_test.keys())[0]

            ut.log(
                "Test with name %s has parameters %s" % (key, univariate_test[key]),
                state,
            )
            if key != "regression":
                variable = univariate_test[key].pop("variable")
                if type(univariate_test[key]) is dict:
                    univariate_test[key]["datasets"] = [
                        list(np.argwhere(covariates_df[variable] == name).flatten() + 1)
                        for name in univariate_test[key]["name"]
                    ]
            else:
                variable = ""
            ut.log(
                "After parsing, test with name %s has parameters %s"
                % (key, univariate_test[key]),
                state,
            )
            univariate_out_dir = os.path.join(
                state["outputDirectory"], "coinstac-univariate-%s-%s" % (key, variable)
            )
            os.makedirs(univariate_out_dir, exist_ok=True)
            if key == "regression":
                univariate_test = univariate_test[key]
            try:
                gift_mancova(
                    ica_param_file=ica_parameters,
                    out_dir=univariate_out_dir,
                    TR=args["input"].get("TR", 2),
                    features=args["input"]["features"],
                    comp_network_names=args["input"].get(
                        "comp_network_names", NEUROMARK_NETWORKS
                    ),
                    univariate_tests=univariate_test,
                    run_name="coinstac-mancovan-univariate",
                    numOfPCs=args["input"].get("numOfPCs", [4, 4, 4]),
                    freq_limits=args["input"].get("freq_limits", [0.1, 0.15]),
                    t_threshold=args["input"].get("t_threshold", 0.05),
                    image_values=args["input"].get("image_values", "positive"),
                    threshdesc=args["input"].get("threshdesc", "none"),
                    p_threshold=args["input"].get("p_threshold", 0.05),
                    display_p_threshold=args["input"].get("display_p_threshold", 0.05),
                )
            except Exception as e:
                ut.log(
                    "Univariate analysis ({key}, {variable}) raised an exception. Full error string {err}".format(
                        key=key, variable=variable, err=str(e)
                    ),
                    state,
                )
    else:
        ut.log("Skipping mancova", state)
    ut.log("Collecting Mancova results", state)
    output_dict = {"computation_phase": "scica_mancova_1"}
    ut.log("Output %s" % (output_dict), state)
    cache_dict = {}
    computation_output = {"output": output_dict, "cache": cache_dict, "state": state}

    return computation_output


def scica_check_out(args):
    state = args["state"]
    gigica_dir = os.path.join(state["outputDirectory"], "gica_gmc_gica_results")
    output_dict = {"gigica_output": None}
    if os.path.exists(gigica_dir):
        output_dict["gigica_output"] = gigica_dir
    cache_dict = {}
    computation_output = {"output": output_dict, "cache": cache_dict, "state": state}
    return computation_output
