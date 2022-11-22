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

def chmod_dir_recursive(dir_name):
    try:
        for dir_root, dirs, files in os.walk(dir_name):
            for d in dirs:
                os.chmod(os.path.join(dir_root, d), 0o777)
            for f in files:
                os.chmod(os.path.join(dir_root, f), 0o777)
    except:
        pass;


def convert_covariates(covariate_df, state, covariate_types=None, N=None):

    out_dir = state["outputDirectory"]
    #shutil.copy(
    #    covariate_filename, os.path.join(out_dir, os.path.basename(covariate_filename))
    #)

    df = covariate_df
    cov_types = covariate_types
    covariates = {}

    ut.log("Covariate Types %s" % (cov_types), state)
    ut.log("Covariate Types name %s" % (cov_types["name"].values), state)

    for covariate_name in df.columns:
        if covariate_types is not None:
            if covariate_name not in cov_types["name"].values:
                continue
        covariate_series = df[covariate_name]
        if covariate_name == "filename" or covariate_name == "niftifilename":
            continue
        if N is not None:
            ut.log("Covariate series for covariate_name %s" % (covariate_name), state)
            covariate_series = covariate_series[:N] if N>0 else covariate_series[:]

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

    if N is not None and N > 0:
        df = df.head(N)

    ut.log("AFTER Covariates are %s -df %s" % (covariates, df), state)

    #print(covariates)
    return covariates, df, cov_types


def extract_covariates(covariate_dict):
    from pandas.api.types import is_numeric_dtype

    # filename key has GICA output filenames; the rest of the columns are
    # covariate columns with : separting the covariate name and its type
    df = pd.DataFrame.from_dict(covariate_dict).T

    # Remove index as it correspond to GICA filenames and not relevant with covariates
    df.reset_index(drop=True, inplace=True)

    # GICA files are more than number of subjects/covariates; when we remove index there will be nan rows.
    df = df.replace(r'^\s*$', np.nan, regex=True);
    df = df.dropna()

    # Correcting datatypes
    covar_df = df.apply(pd.to_numeric, errors='ignore')

    # extract covariate keys - name and its type
    cov_types = {"name": [], "type": []}
    rename_cols = {}
    for colname in df.columns:
        if 'filename' not in colname:
            vals = colname.strip().split(':')
            pd_col_type = 'continuous' if is_numeric_dtype(df.dtypes[colname]) else 'categorical'
            if len(vals) > 0:
                covar_name = vals[0].strip();covar_type = vals[1].strip()
                rename_cols[colname] = covar_name
                cov_types['name'].append(covar_name)
                cov_types['type'].append(covar_type)

                # check column type in pandas, correct column type if string column is converted to numeric
                if covar_type == 'categorical' and pd_col_type != covar_type:
                    covar_df[colname] = covar_df[colname].astype(str)

            else:
                cov_types['name'].append(vals)
                cov_types['type'].append(pd_col_type)
    covar_df.rename(rename_cols, axis='columns', inplace=True)

    return covar_df, pd.DataFrame(cov_types)



def local_run_mancova(args):
    state = args["state"]
    ut.log("Got input %s" % (args["input"]), state)

    args["cache"]["skip_gica"] = args["input"]["skip_gica"]
    args["cache"]["gica_input_dir"] = args["input"]["gica_input_dir"]
    args["cache"]["univariate_test_list"] = args["input"]["univariate_test_list"]

    covariates_df, covariate_types_df = extract_covariates(args["input"]["covariates"])
    file_list=covariates_df['niftifilename'].tolist()

    covariates, covariates_df, covariate_types_df = convert_covariates(
        covariates_df, state, covariate_types_df, N=len(file_list)
    )

    maskfile = args["input"]["mask"]
    TR = args["input"]["TR"]

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

        in_files = [os.path.join(state["baseDirectory"], f) for f in file_list]
        ut.log("Loaded files %s" % ", ".join(in_files), state)

        template=args["input"]["scica_template"]
        ut.log("Interpolated template at file %s" % template, state)
        ut.log("Running group ICA", state)
        gica_out_dir = os.path.join(state["outputDirectory"], "coinstac-gica")
        os.makedirs(gica_out_dir, exist_ok=True)
        curr_TR = args["input"].get("TR", [2])
        curr_TR = curr_TR if type(curr_TR) is list else [curr_TR]

        output = gift_gica(
            in_files=in_files,
            refFiles=template,
            mask=maskfile,
            out_dir=gica_out_dir,
            group_pca_type="subject specific",
            algoType=16,
            run_name="coinstac-gica",
            scaleType=2,
            TR=curr_TR,
        )
        baseDir=gica_out_dir

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

        ut.log("ICA params file" + ''.join(ica_parameters),  state)

        if len(ica_parameters) == 0:
            raise (ValueError("Ica parameter files could not be found in the outupt"))
        ut.log("Running univariate tests", state)
        univariate_test_list = copy.deepcopy(args["input"]["univariate_test_list"])
        for univariate_test in univariate_test_list:
            write_stats_info = 0
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
                write_stats_info = 1
            mancova_params_dict={'ica_param_file':ica_parameters,
                                        'out_dir':univariate_out_dir,
                                        'TR':args["input"].get("TR", 2),
                                        'features':args["input"]["features"],
                                        'comp_network_names':args["input"].get("comp_network_names", NEUROMARK_NETWORKS),
                                        'covariates':covariates,
                                        'univariate_tests':univariate_test,
                                        'run_name':"coinstac-mancovan-univariate",
                                        'numOfPCs':args["input"].get("numOfPCs", [4, 4, 4]),
                                        'freq_limits':args["input"].get("freq_limits", [0.1, 0.15]),
                                        't_threshold':args["input"].get("t_threshold", 0.05),
                                        'image_values':args["input"].get("image_values", "positive"),
                                        'threshdesc':args["input"].get("threshdesc", "fdr"),
                                        'p_threshold':args["input"].get("p_threshold", 0.05),
                                        'display_p_threshold':args["input"].get("display_p_threshold", 0.05),
                                        'display_local_result_summary':False
                                     }
            ut.log("Calling GIFT MANCOVA Univariate with params: %s"% (str(mancova_params_dict)), state )
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
                    display_local_result_summary=False,
                    write_stats_info=write_stats_info
                )
            ut.log("Done with  GIFT MANCOVA call", state)
            transfer_univariate_stats(state, univariate_out_dir)
    else:
        ut.log("Skipping univariate tests", state)

    output_dict = dict(
        computation_phase="scica_mancova_1",
        univariate_test_list=args["input"]["univariate_test_list"],
        covariates=covariates,
        covariates_df=covariates_df.to_dict(),
        covariate_types=covariate_types_df.to_dict(),
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


def transfer_univariate_stats(state, univariate_out_dir):
    ut.log("Extracting stats info files", state)
    stats_files = [f for f in list(
                        glob.iglob(
                            os.path.join(univariate_out_dir, "**", "*mancovan_stats_info.mat"),
                            recursive=True,
                        )) if os.path.exists(f)]
    test_name=os.path.basename(univariate_out_dir)
    ut.log("mancovan_stats_info files generated: %s" % (''.join(stats_files)), state)

    ut.log("Copying mancovan_stats_info to remote", state)
    for f in stats_files:
        shutil.copy(
            f, os.path.join(state["transferDirectory"],
              os.path.basename(f))
        )


def scica_check_out(args):
    state = args["state"]
    gigica_dir = os.path.join(state["outputDirectory"], "gica_gmc_gica_results")
    output_dict = {"gigica_output": None}
    if os.path.exists(gigica_dir):
        output_dict["gigica_output"] = gigica_dir
    cache_dict = {}
    computation_output = {"output": output_dict, "cache": cache_dict, "state": state}
    return computation_output
