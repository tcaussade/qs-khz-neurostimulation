import json


def load_params_json(jsonfile, user_param):
    # This function is intended to update user's param_dict with default
    # parameters read from a file. The idea is not to replace user's
    # parameters but fill with defaults those missing parameters

    with open(jsonfile, 'r') as fp:
        defaults = json.load(fp)

    context = defaults.copy()
    context.update(user_param)
    return context


def load_defaults(sim_param=None, waveform_param=None, field_param=None):

    if sim_param is not None:
        sim_param = load_params_json(
            "./default_sim_param.json", sim_param)
    if waveform_param is not None:
        waveform_param = load_params_json(
            "./default_waveform_param.json", waveform_param)
    if field_param is not None:
        field_param = load_params_json(
            "./default_field_param.json", field_param)

    return sim_param, waveform_param, field_param
