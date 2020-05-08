import sys
import ujson as json
from .BackRecon import gift_gica
from utils import listRecursive
import utils as ut
import os
import glob

TC_SEARCH_STRING = 'gica_cmd_sub*%d_timecourses_ica_s1_.nii'


def scica_local_1(args):
    state = args["state"]
    in_files = [os.path.join(state['baseDirectory'], f)
                for f in args["input"]["data"]]
    maskfile = os.path.join('/computation', 'local_data', 'mask.nii')
    template = os.path.join('/computation', 'local_data', 'NeuroMark.nii')
    template = ut.get_interpolated_nifti(
        in_files[0], template, destination_dir=state["outputDirectory"])
    pyscript = os.path.join(state["outputDirectory"], 'pyscript_gicacommand.m')
    if os.path.exists(pyscript):
        os.remove(pyscript)
    output = gift_gica(
            in_files=in_files,
            refFiles=[template],
            mask=maskfile,
            out_dir=state["outputDirectory"],
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

    output_dict = {
        'subject_sms': subject_sms,
        'subject_tcs': subject_tcs,
        'other_matfiles': other_matfiles,
        'computation_phase': 'scica_local_1'
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
    gigica_dir = os.path.join(state["outputDirectory"], "gica_gmc_gica_results")
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
       

if __name__ == '__main__':
    parsed_args = json.loads(sys.stdin.read())
    phase_key = list(listRecursive(parsed_args, 'computation_phase'))

    if not phase_key:
        computation_output = scica_local_1(parsed_args)
        sys.stdout.write(computation_output)
    else:
        raise ValueError("Error occurred at Local")
