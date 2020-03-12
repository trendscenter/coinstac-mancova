import sys
import ujson as json
from utils import listRecursive


def scica_remote_noop(args):
    output_dict = dict(computation_phase="scica_remote_noop")
    cache_dict = {}
    computation_output = {
        "output": output_dict,
        "cache": cache_dict,
        "success": True
    }
    return json.dumps(computation_output)


if __name__ == '__main__':
    parsed_args = json.loads(sys.stdin.read())
    phase_key = list(listRecursive(parsed_args, 'computation_phase'))

    if phase_key == 'scica_local_1':
        computation_output = spatially_constrained_ica_remote_noop(parsed_args)
        sys.stdout.write(computation_output)
    else:
        raise ValueError("Error occurred at Local")
