import sys
import ujson as json
import hdfdict

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
import base64

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
    out_dir = state["transferDirectory"]
    shutil.copy(
        covariate_filename, os.path.join(out_dir, os.path.basename(covariate_filename))
    )
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
    ut.log("Covariates are %s -df %s" % (covariates, df), state)
    if N is not None:
        df = df.head(N)
    return covariates, df, cov_types


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
    covariates, covariates_df, covariate_types = convert_covariates(
        covariate_file, state, covariate_types=covariate_type_file, N=len(in_files)
    )

    maskfile = args["input"]["mask"]

    pyscript = os.path.join(state["transferDirectory"], "pyscript_gicacommand.m")
    if os.path.exists(pyscript):
        os.remove(pyscript)
    if (
        not os.path.exists(
            os.path.join(state["baseDirectory"], args["input"]["gica_input_dir"])
        )
        or len(args["input"]["gica_input_dir"]) == 0
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
        """
        ut.log("Interpolating", state)
        template = ut.get_interpolated_nifti(
            in_files[0],
            args["input"]["scica_template"],
            destination_dir=state["outputDirectory"],
        )
        """
        template=args["input"]["scica_template"]
        ut.log("Interpolated template at file %s" % template, state)
        ut.log("Running group ICA", state)
        gica_out_dir = os.path.join(state["transferDirectory"], "coinstac-gica")
        os.makedirs(gica_out_dir, exist_ok=True)
        os.makedirs(
            os.path.join(state["outputDirectory"], "coinstac-gica"), exist_ok=True
        )
        output = gift_gica(
            in_files=in_files,
            refFiles=template,
            mask=maskfile,
            out_dir=gica_out_dir,
            group_pca_type="subject specific",
            algoType=16,
            run_name="coinstac-gica",
        )
        ut.log(
            "Copying  output from %s to %s"
            % (gica_out_dir, os.path.join(state["outputDirectory"]),),
            state,
        )
        for filename in glob.glob(gica_out_dir):
            if os.path.isdir(filename):
                shutil.copytree(
                    filename,
                    os.path.join(
                        state["outputDirectory"],
                        "coinstac-gica",
                        os.path.basename(filename),
                    ),
                )
            else:
                shutil.copy(
                    filename,
                    os.path.join(
                        state["outputDirectory"],
                        "coinstac-gica",
                        os.path.basename(filename),
                    ),
                )
    else:
        ut.log(
            "Copying preexisting output from %s to %s"
            % (
                os.path.join(state["baseDirectory"], args["input"]["gica_input_dir"]),
                os.path.join(
                    state["transferDirectory"], args["input"]["gica_input_dir"]
                ),
            ),
            state,
        )
        shutil.copytree(
            os.path.join(state["baseDirectory"], args["input"]["gica_input_dir"]),
            os.path.join(state["transferDirectory"], args["input"]["gica_input_dir"])
        )

    output_dict = dict(
        computation_phase="scica_mancova_1",
        univariate_test_list=univariate_test_list,
        covariates=covariates,
        covariates_df=covariates_df.to_dict(),
        covariate_types=covariate_types.to_dict(),
        TR=args["input"].get("TR", 2),
        features=args["input"]["features"],
        interactions=args["input"].get("interactions", []),
        numOfPCs=args["input"].get("numOfPCs", [4, 4, 4]),
        run_name="coinstac-mancovan-multivariate",
        comp_network_names=args["input"].get("comp_network_names", NEUROMARK_NETWORKS),
        freq_limits=args["input"].get("freq_limits", [0.1, 0.15]),
        t_threshold=args["input"].get("t_threshold", 0.05),
        image_values=args["input"].get("image_values", "positive"),
        threshdesc=args["input"].get("threshdesc", "fdr"),
        p_threshold=args["input"].get("p_threshold", 0.05),
        display_p_threshold=args["input"].get("display_p_threshold", 0.05),
        run_univariate_tests=args["input"].get("run_univariate_tests", True),
        run_mancova=args["input"].get("run_mancova", True),
    )
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
