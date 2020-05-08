#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 16:07:00 2018 (MDT)

@author: Rogers F. Silva
"""

import json
import os
import sys
import numpy as np
import utils as ut


def drm_remote_1(args):

    state = args['state']
    inputs = args['input']
    cache = args['cache']

    # Start remote computation:

    # Compute global row sum over all local row_sums
    row_sum_global = np.hstack(
        tuple(np.array(su['row_sum']) for (localID, su) in inputs.items())
    ).sum(axis=1)

    # Compoute global number of columns over all local num_cols
    num_cols_global = np.array([su['num_cols']
                                for (localID, su) in inputs.items()]).sum()

    # Compute global row mean
    row_mean_global = row_sum_global/num_cols_global

    # Compile results to be transmitted to local sites and cached for reuse in next iteration
    computation_output = {
        "state": state,
        "output": {
            "row_mean_global": row_mean_global.tolist(),
            "computation_phase": "drm_remote_1"
        },
        "cache": cache,

    }

    #output_file = os.path.join(state['outputDirectory'], 'row_mean_global.data')
    #np.savetxt(output_file, row_mean_global, fmt='%.6f')

    #import pdb; pdb.set_trace()
    return computation_output


if __name__ == '__main__':

    parsed_args = json.loads(sys.stdin.read())
    phase_key = list(ut.listRecursive(parsed_args, 'computation_phase'))

    if not phase_key:
        raise ValueError(
            "Error occurred at Remote: missing phase key from local site(s).")
    elif "local_1" in phase_key:
        computation_output = drm_remote_1(parsed_args)
        # Transmit results to remote
        # as file (for large volumes of data; OS overhead):
        # as JSON string (for smaller volumes of data; JSON conversion overhead):
        sys.stdout.write(json.dumps(computation_output))
    else:
        raise ValueError("Error occurred at Remote: unknown phase key.")
