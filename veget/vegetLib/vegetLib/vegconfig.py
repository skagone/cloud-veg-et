import os
import sys
import yaml

def return_veget_params(config_directory):
        
    # this allows for the config to be created from a preexisting file
    config_run_file_path = os.path.join(config_directory,'run_param.yml')
    if os.path.exists(config_run_file_path):
        with open(config_run_file_path, 'r') as cfgpath:
            run_config_dict = yaml.safe_load(cfgpath)
            # print(run_config_dict)
    else:
        print('the path does not exist, check the path you gave VegET()')
        sys.exit(1)

    config_model_file_path = os.path.join(config_directory, 'model_param.yml')
    if os.path.exists(config_model_file_path):
        with open(config_model_file_path, 'r') as cfgpath:
            model_config_dict = yaml.safe_load(cfgpath)
            # print(model_config_dict)
    else:
        print('the path does not exist, check the path you gave VegET()')
        sys.exit(1)

    config_path_file_path = os.path.join(config_directory, 'path_param.yml')
    if os.path.exists(config_path_file_path):
        with open(config_path_file_path, 'r') as cfgpath:
            path_config_dict = yaml.safe_load(cfgpath)
            # print(path_config_dict)
    else:
        print('the path does not exist, check the path you gave VegET()')
        sys.exit(1)

    config_dict = {**run_config_dict, **model_config_dict, **path_config_dict}
    # print(config_dict)

    return config_dict
