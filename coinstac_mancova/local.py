import sys
import ujson as json
import hdfdict

# import json
from .run_gift import gift_mancova, gift_gica, gift_run_matlab_script
from utils import listRecursive
import utils as ut
import os
import glob
import copy
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
    #out_dir = state["transferDirectory"]
    out_dir = state["outputDirectory"]
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
            ut.log("Covariate series for covariate_name %s" % (covariate_name), state)
            if N>0 :
                covariate_series = covariate_series[:N]
            else:
                covariate_series = covariate_series[:]

        if covariate_name in cov_types["name"].values:
            cov_type = cov_types["type"][list(cov_types["name"].values).index(covariate_name)]
        else:
            cov_type = "continuous"

        ut.log("Covariate name %s , type %s " % (covariate_name, cov_type), state)
        ut.log("Covariate values %s " % (",".join([str(s) for s in list(covariate_series)])), state)
        fname = os.path.join(out_dir, "COINSTAC_COVAR_%s.txt" % covariate_name)
        with open(fname, "w") as file:
            file.write("\n".join([str(s) for s in list(covariate_series)]))
        ut.log("Wrote covariates %s to file %s" % (covariate_name, fname), state)
        covariates[covariate_name] = [cov_type, fname]
    ut.log("Covariates are %s -df %s" % (covariates, df), state)
    if N is not None and N>0:
        df = df.head(N)
    ut.log("AFTER Covariates are %s -df %s" % (covariates, df), state)
    return covariates, df, cov_types

def chmod_dir_recursive(dir_name):
    try:
        for dir_root, dirs, files in os.walk(dir_name):
            for d in dirs:
                os.chmod(os.path.join(dir_root, d), 0o777)
            for f in files:
                os.chmod(os.path.join(dir_root, f), 0o777)
    except:
        #ut.log("Error while changing permissions for files in dir: %s"%dir_name)
        a=""

def local_run_mancova(args):
    state = args["state"]
    ut.log("Got input %s" % (args["input"]), state)

    args["cache"]["skip_gica"] = args["input"]["skip_gica"]
    args["cache"]["gica_input_dir"] = args["input"]["gica_input_dir"]
    args["cache"]["univariate_test_list"] = args["input"]["univariate_test_list"]

    cov_filename = [i for i in args["input"]["data"] if "covariates.csv" in i]
    ctype_filename = [i for i in args["input"]["data"] if "covariate_keys.csv" in i]
    covariate_file = os.path.join(state["baseDirectory"], cov_filename[0])
    covariate_type_file = os.path.join(state["baseDirectory"], ctype_filename[0])
    ut.log("Covariate File Name:" + covariate_file, state)
    file_list = args["input"]["data"]
    file_list.remove(cov_filename[0])
    file_list.remove(ctype_filename[0])

    in_files = [os.path.join(state["baseDirectory"], f) for f in file_list]
    ut.log("Loaded files %s" % ", ".join(in_files), state)
    covariates, covariates_df, covariate_types = convert_covariates(
        covariate_file, state, covariate_types=covariate_type_file, N=len(in_files)
    )

    maskfile = args["input"]["mask"]

    TR = args["input"]["TR"]

    #pyscript = os.path.join(state["transferDirectory"], "pyscript_gicacommand.m")
    pyscript = os.path.join(state["outputDirectory"], "pyscript_gicacommand.m")
    if args["cache"]["skip_gica"] is False:
        if os.path.exists(pyscript):
            os.remove(pyscript)
        if (
            not os.path.exists(
                os.path.join(state["baseDirectory"], args["cache"]["gica_input_dir"])
            )
            or len(args["cache"]["gica_input_dir"]) == 0
        ):
            ut.log(
                "The input directory %s, does not exist"
                % (
                    str(
                        os.path.join(
                            state["baseDirectory"], args["cache"]["gica_input_dir"]
                        )
                    )
                ),
                state,
            )
            args["cache"]["gica_input_dir"] = None
    if args["cache"]["gica_input_dir"] is None:
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
            scaleType=2,
            TR=args["input"].get("TR", 2),
        )
        baseDir=gica_out_dir
        """
        # Copying output of coinstac-gica to transfer directory.
        os.makedirs(
            os.path.join(state["transferDirectory"], "coinstac-gica"))
        ut.log(
            "Copying  output from %s to %s"
            % (gica_out_dir, os.path.join(state["transferDirectory"]),),
            state,
        )
        for filename in glob.glob(gica_out_dir):
            if os.path.isdir(filename):
                shutil.copytree(
                    filename,
                    os.path.join(
                        state["transferDirectory"],
                        "coinstac-gica",
                        os.path.basename(filename),
                    ),
                )
            else:
                shutil.copy(
                    filename,
                    os.path.join(
                        state["transferDirectory"],
                        "coinstac-gica",
                        os.path.basename(filename),
                    ),
                )
        """
    else:
        if args["cache"]["skip_gica"] is False:
            baseDir = os.path.join(state["baseDirectory"], args["cache"]["gica_input_dir"])
        else:
            baseDir = os.path.join(state["baseDirectory"], args["input"]["gica_input_dir"])
        ut.log(
                "Using preexisting GICA output from %s "
                % ( baseDir ),
                state,
            )
        """
        ut.log(
            "Copying preexisting output from %s to %s"
            % (
                baseDir,
                os.path.join(
                    state["transferDirectory"], args["cache"]["gica_input_dir"]
                ),
            ),
            state,
        )
        if os.path.isdir(os.path.join(state["transferDirectory"], args["cache"]["gica_input_dir"])) is False:
            shutil.copytree(
                baseDir,
                os.path.join(state["transferDirectory"], args["cache"]["gica_input_dir"])
            )
        """

    ut.log( "Done with GICA ", state)

    if args["input"]["run_univariate_tests"]:
        stat_results = dict()
        ica_parameters = [f for f in list(
            glob.iglob(
                os.path.join(baseDir, "**", "*parameter_info.mat"),
                recursive=True,
            )
        ) if os.path.exists(f)]
        ica_parameters.sort()

        ut.log("Running univariate tests", state)
        univariate_test_list = copy.deepcopy(args["input"]["univariate_test_list"])
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
            if univariate_out_dir[-1] == "-":
                univariate_out_dir = univariate_out_dir[: (len(univariate_out_dir) - 1)]
            stat_results[univariate_out_dir] = dict()
            ut.log(
                "threshdesc setting:"+json.dumps(args["input"]['threshdesc']),
                state,
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
                    covariates=covariates,
                    univariate_tests=univariate_test,
                    run_name="coinstac-mancovan-univariate",
                    numOfPCs=args["input"].get("numOfPCs", [4, 4, 4]),
                    freq_limits=args["input"].get("freq_limits", [0.1, 0.15]),
                    t_threshold=args["input"].get("t_threshold", 0.05),
                    image_values=args["input"].get(
                        "image_values", "positive"
                    ),
                    threshdesc=args["input"].get("threshdesc", "fdr"),
                    p_threshold=args["input"].get("p_threshold", 0.05),
                    display_p_threshold=args["input"].get(
                        "display_p_threshold", 0.05
                    ),
                )
                stat_results[univariate_out_dir][key] = list(
                    glob.glob(os.path.join(univariate_out_dir, "**", "*.html"))
                )
                #shutil.copytree(univariate_out_dir, os.path.join(state["transferDirectory"], os.path.basename(univariate_out_dir)))
                extract_univariate_stats(state, args["input"]["features"], univariate_out_dir)

            except Exception as e:
                ut.log(
                    "Univariate analysis ({key}, {variable}) raised an exception. Full error string {err}".format(
                        key=key, variable=variable, err=str(e)
                    ),
                    state,
                )
    else:
        ut.log("Skipping univariate tests", state)



    output_dict = dict(
        computation_phase="scica_mancova_1",
        univariate_test_list=args["input"]["univariate_test_list"],
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

    # Gives permission errors while copying if file already exists
    chmod_dir_recursive(state["transferDirectory"])
    chmod_dir_recursive(state["outputDirectory"])

    return computation_output


def extract_univariate_stats(state, features, univariate_out_dir):
    import hdf5storage
    matfiledata = {}
    matfiledata[u'clientID'] = state["clientId"]
    matfiledata[u'outputDirs'] = [state["transferDirectory"], state["outputDirectory"]]
    matfiledata[u'featureStats'] = features
    matfiledata[u'univariateResultsDir'] = univariate_out_dir
    matfiledata[u'mancovanInputFilesPrefix'] = 'coinstac-gica_'#'rest_hcp_' #'coinstac-gica_' OR 'coinstac-gica_merge_'
    output_file = '/computation/'+state["clientId"]+'_mancova_structure_info.mat';

    ut.log("Writing mat data %s" % (matfiledata), state)
    hdf5storage.write(matfiledata, '', output_file, matlab_compatible=True)
    ut.log("Done writing mat data", state)

    #orig_mat_script='/computation/local_mat_script.m'
    orig_mat_script='/computation/coinstac_mancova/matcode/local_mat_script.m'
    temp_mat_script='/computation/temp_%s_mat_script.m'%(state["clientId"])
    temp_mat_script=os.path.join(state["cacheDirectory"], "temp_%s_mat_script.m"%(state["clientId"]))
    # open both files
    with open(orig_mat_script, 'r') as orig_file, open(temp_mat_script, 'a') as tmp_file:
        tmp_file.write("mat_info=load('"+output_file+"');")
        # read content from first file
        for line in orig_file:
            # append content to second file
            tmp_file.write(line)

    ut.log("Running local matlab script: %s"%temp_mat_script, state)
    gift_run_matlab_script(temp_mat_script)



def scica_check_out(args):
    state = args["state"]
    gigica_dir = os.path.join(state["outputDirectory"], "gica_gmc_gica_results")
    output_dict = {"gigica_output": None}
    if os.path.exists(gigica_dir):
        output_dict["gigica_output"] = gigica_dir
    cache_dict = {}
    computation_output = {"output": output_dict, "cache": cache_dict, "state": state}
    return computation_output
