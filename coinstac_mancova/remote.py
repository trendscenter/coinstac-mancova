import sys
import os
import glob
import pandas as pd
import shutil
import ujson as json
from utils import listRecursive
import utils as ut
from .run_gift import gift_mancova, gift_gica


NEUROMARK_NETWORKS = {
    "SC": [1, 2, 3, 4, 5],
    "AUD": [6, 7],
    "SM": [8, 9, 10, 11, 12, 13, 14, 15, 16],
    "VIS": [17, 18, 19, 20, 21, 22, 23, 24, 25],
    "CC": [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42],
    "DMN": [43, 44, 45, 46, 47, 48, 49],
    "CR": [50, 51, 52, 53],
}


import numpy as np


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
    ut.log("Covariate Types types %s" % (cov_types["type"]), state)
    ut.log("Covariate df columns %s" % (df.columns), state)
    for covariate_name in df.columns:
        if covariate_types is not None:
            if covariate_name not in cov_types["name"].values:
                continue
        covariate_series = df[covariate_name]
        if covariate_name == "filename":
            continue
        if N is not None:
            covariate_series = covariate_series[:N]
        ut.log(
            "Checking logic %s %s"
            % (covariate_name, covariate_name in cov_types["name"].values),
            state,
        )
        if covariate_name in cov_types["name"].values:
            cov_type = cov_types["type"][list(cov_types["name"].values).index(covariate_name)]
        else:
            cov_type = "continuous"
        fname = os.path.join(out_dir, "COINSTAC_COVAR_%s.txt" % covariate_name)
        with open(fname, "w") as file:
            file.write("\n".join([str(s) for s in list(covariate_series)]))
        ut.log("Wrote covariates %s to file %s" % (covariate_name, fname), state)
        covariates[covariate_name] = [cov_type, fname]
    ut.log("Covariates are %s" % (covariates), state)
    if N is not None:
        df = df.head(N)
    return covariates, df


def mancova_aggregate(args):
    inputs = args["input"]
    state = args["state"]

    ut.log("Checking the inputs on remote: " + str(inputs), state)
    ut.log("Checking the state: " + str(state), state)
    ica_parameter_list = [f for f in list(
        glob.iglob(
            os.path.join(state["baseDirectory"], "**", "*parameter_info.mat"),
            recursive=True,
        )
    ) if os.path.exists(f)]
    ica_parameter_list.sort()
    ut.log("Checking ICA parameters %s" % (str(ica_parameter_list)), state)
    if len(ica_parameter_list) == 0:
        raise (ValueError("Ica parameter files could not be found in the outupt"))
    covariate_file = os.path.join(state["outputDirectory"], "covariates.csv")
    covariate_type_file = os.path.join(state["outputDirectory"], "covariate_types.csv")
    covariate_dfs = [pd.DataFrame(c["covariates_df"]) for c in inputs.values()]
    covariate_types = [pd.DataFrame(c["covariate_types"]) for c in inputs.values()][0]
    univariate_test_list = [c["univariate_test_list"] for c in inputs.values()][0]
    covariate_types.to_csv(covariate_type_file)
    covariate_df = pd.concat(covariate_dfs)
    covariate_df.to_csv(covariate_file, index=False)
    covariates, covariates_df = convert_covariates(
        covariate_file, state, covariate_types=covariate_type_file, N=len(covariate_df)
    )

    ica_parameters = ica_parameter_list
    ut.log("ICA Parameters file is %s" % (str(ica_parameters)), state)
    ut.log("Checking covariates %s" % (str(covariates)), state)
    stat_results = dict()
    if args["input"]["local0"]["run_mancova"]:
        ut.log("Running Mancova", state)
        try:
            multivariate_out_dir = os.path.join(
                state["outputDirectory"], "coinstac-multivariate"
            )
            stat_results[multivariate_out_dir] = dict()
            os.makedirs(multivariate_out_dir, exist_ok=True)
            gift_mancova(
                ica_param_file=ica_parameters,
                out_dir=multivariate_out_dir,
                TR=args["input"]["local0"].get("TR", 2),
                features=args["input"]["local0"]["features"],
                covariates=covariates,
                interactions=args["input"]["local0"].get("interactions", []),
                numOfPCs=args["input"]["local0"].get("numOfPCs", [4, 4, 4]),
                run_name="coinstac-mancovan-multivariate",
                comp_network_names=args["input"]["local0"].get(
                    "comp_network_names", NEUROMARK_NETWORKS
                ),
                freq_limits=args["input"]["local0"].get("freq_limits", [0.1, 0.15]),
                t_threshold=args["input"]["local0"].get("t_threshold", 0.05),
                image_values=args["input"]["local0"].get("image_values", "positive"),
                threshdesc=args["input"]["local0"].get("threshdesc", "fdr"),
                p_threshold=args["input"]["local0"].get("p_threshold", 0.05),
                display_p_threshold=args["input"]["local0"].get(
                    "display_p_threshold", 0.05
                ),
            )
            stat_results[multivariate_out_dir] = list(
                glob.glob(os.path.join(multivariate_out_dir, "**", "*.html"))
            )
            shutil.copytree(multivariate_out_dir, os.path.join(state["transferDirectory"], os.path.basename(multivariate_out_dir))
        except Exception as e:
            ut.log(
                "Multivariate analysis raised an exception, likely because of bad conditioning. More subjects are required. Full error string {err}".format(
                    err=str(e)
                ),
                state,
            )
    else:
        ut.log("Skipping Mancova", state)
    if args["input"]["local0"]["run_univariate_tests"]:
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
            if univariate_out_dir[-1] == "-":
                univariate_out_dir = univariate_out_dir[: (len(univariate_out_dir) - 1)]
            stat_results[univariate_out_dir] = dict()

            os.makedirs(univariate_out_dir, exist_ok=True)
            if key == "regression":
                univariate_test = univariate_test[key]
            try:
                gift_mancova(
                    ica_param_file=ica_parameters,
                    out_dir=univariate_out_dir,
                    TR=args["input"]["local0"].get("TR", 2),
                    features=args["input"]["local0"]["features"],
                    comp_network_names=args["input"]["local0"].get(
                        "comp_network_names", NEUROMARK_NETWORKS
                    ),
                    covariates=covariates,
                    univariate_tests=univariate_test,
                    run_name="coinstac-mancovan-univariate",
                    numOfPCs=args["input"]["local0"].get("numOfPCs", [4, 4, 4]),
                    freq_limits=args["input"]["local0"].get("freq_limits", [0.1, 0.15]),
                    t_threshold=args["input"]["local0"].get("t_threshold", 0.05),
                    image_values=args["input"]["local0"].get(
                        "image_values", "positive"
                    ),
                    threshdesc=args["input"]["local0"].get("threshdesc", "none"),
                    p_threshold=args["input"]["local0"].get("p_threshold", 0.05),
                    display_p_threshold=args["input"]["local0"].get(
                        "display_p_threshold", 0.05
                    ),
                )
                stat_results[univariate_out_dir][key] = list(
                    glob.glob(os.path.join(univariate_out_dir, "**", "*.html"))
                )
                shutil.copytree(univariate_out_dir, os.path.join(state["transferDirectory"], os.path.basename(univariate_out_dir))
            except Exception as e:
                ut.log(
                    "Univariate analysis ({key}, {variable}) raised an exception. Full error string {err}".format(
                        key=key, variable=variable, err=str(e)
                    ),
                    state,
                )
    else:
        ut.log("Skipping mancova", state)
    output_dict = {
        "computation_phase": "scica_mancova_remote",
        "stat_results": stat_results,
    }
    ut.log("Output %s" % (output_dict), state)
    cache_dict = {}
    computation_output = {"output": output_dict, "cache": cache_dict, "state": state}
    return computation_output
