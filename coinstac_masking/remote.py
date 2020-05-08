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


def masking_remote_1(args):

    state = args['state']
    inputs = args['input']
    cache = args['cache']

    # Start remote computation:
    computation_output = {
        "output": {
            "computation_phase": "masking_remote_1"
        },
        "cache": dict(),
        "success": True
    }

    #import pdb; pdb.set_trace()
    return computation_output


if __name__ == '__main__':

    parsed_args = json.loads(sys.stdin.read())
    phase_key = list(ut.listRecursive(parsed_args, 'computation_phase'))

    if not phase_key:
        raise ValueError(
            "Error occurred at Remote: missing phase key from local site(s).")
    elif "local_1" in phase_key:
        computation_output = masking_remote_1(parsed_args)
        # Transmit results to remote
        # as file (for large volumes of data; OS overhead):
        # as JSON string (for smaller volumes of data; JSON conversion overhead):
        sys.stdout.write(json.dumps(computation_output))
    else:
        raise ValueError("Error occurred at Remote: unknown phase key.")
