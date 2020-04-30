import sys
import ujson as json
from utils import listRecursive


def aggregate_sub_stats(state_type, output_so_far):
    variable_keys = ['mult', 'uni']
    for vk in variable_keys:
        if vk not in output_so_far.keys():
            output_so_far[vk] = {'X': state_type[vk]['X'], 'p': state_type[vk]
                                 ['p'], 't': state_type[vk]['t'], 'stats': state_type[vk]['stats']}
        else:
            output_so_far[vk] = {'X': (output_so_far[vk]['X']+state_type[vk]['X'])/2,
                                 'p': (output_so_far[vk]['p']+state_type[vk]['p'])/2,
                                 't': (output_so_far[vk]['t']+state_type[vk]['t'])/2,
                                 'stats': {k: (v+state_type[vk]['stats'][k])/2 for k, v in output_so_far[vk]['stats'].items()}}
    return output_so_far


def aggregate_stats(all_stats):
    output = {'fnc': {}, 'fnc_avg': {}, 'spectra': {}, 'sm': {}}
    for localID, local_results in all_stats.items():
        local_fnc = local_results['fnc']
        output['fnc'] = aggregate_sub_stats(local_fnc, output['fnc'])
        local_fnc_avg = local_results['fnc_avg']
        output['fnc_avg'] = aggregate_sub_stats(
            local_fnc_avg, output['fnc_avg'])
        spectra = local_results['spectra']
        for sk in spectra.keys():
            if sk not in output['spectra'].keys():
                output['spectra'][sk] = {}
            output['spectra'][sk] = aggregate_sub_stats(
                spectra, output['spectra'][sk])
        sms = local_results['sm']
        for sk in sms.keys():
            if sk not in output['sm'].keys():
                output['sm'][sk] = {}
            output['sm'][sk] = aggregate_sub_stats(sms, output['sm'][sk])
    return output


def mancova_aggregate(args):
    inputs = args['input']
    output_dict = dict(computation_phase="mancova_remote_aggregate")
    aggregated = aggregate_stats(inputs)
    for k, v in aggregated:
        output_dict[k] = v
    cache_dict = {}
    computation_output = {
        "output": output_dict,
        "cache": cache_dict,
        "success": True
    }
    return json.dumps(computation_output)
