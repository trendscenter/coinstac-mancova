#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 16:08:00 2018 (MDT)

@author: Rogers F. Silva
"""

import json
import os
import sys
import numpy as np
import utils as ut


def mask_nifti(data, mask):
    return data[mask == 1, :]


def create_mask(data):
    mask = np.ones((data.shape[0], 1))
    for t in range(data.shape[1]):
        mask = mask * np.array(data[:, t] >= np.mean(data[:, t]), dtype=np.int)
    return mask


def apply_mask(data, mask, state):
    ut.log("Data Size: %s, Mask Shape: %s" %
           (str(data.shape), str(mask.shape)), state)
    D = data[mask.flatten() == 1, :]
    ut.log("Masked Data Size: %s, Mask Shape: %s" %
           (str(D.shape), str(mask.shape)), state)
    ut.log("Sum Masked Data: %s, Sum Mask: %s" %
           (str(sum(sum(D))), str(sum(sum(mask)))), state)

    return D


def masking_local_1(args):

    state = args['state']
    inputs = args['input']
    cache = args['cache']
    files_loaded = inputs["datasets"]

    mask_file = inputs['mask'][0]
    mask_file_type = "nii"
    mask = ut.read_data([os.path.join(mask_file)],
                        mask_file_type, state["clientId"])["0"]
    flat_mask = ut.flatten_data(mask, state)
    masked_data = {idx: apply_mask(ut.flatten_data(data, state), flat_mask, state)
                   for idx, data in files_loaded.items()}

    # Compile results to be transmitted to remote and cached for reuse in next iteration
    computation_output = {
        "state": state,
        "output": {
            "datasets": masked_data,
            "computation_phase": 'masking_local_1'
        },
        "cache": cache,
    }

    return computation_output


def local_2(args):
    state = args['state']
    inputs = args['input']
    cache = args['cache']
    mask = args['mask']

    file_list = inputs['data'][0]
    file_list = [
        os.path.join(state["baseDirectory"], file) for file in file_list
    ]
    file_type = inputs['data'][1][0]
    files_loaded = [ut.read_data(f, file_type) for f in file_list]


if __name__ == '__main__':
    parsed_args = json.loads(sys.stdin.read())
    phase_key = list(ut.listRecursive(parsed_args, 'computation_phase'))

    if not phase_key:
        computation_output = masking_local_1(parsed_args)
        # Transmit results to remote
        # as file (for large volumes of data, OS overhead):
        # as JSON string (for smaller volumes of data, JSON conversion overhead):
        sys.stdout.write(json.dumps(computation_output))
    else:
        raise ValueError("Error occurred at Local")
