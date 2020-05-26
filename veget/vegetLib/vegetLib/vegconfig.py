import os
import yaml

def return_veget_params(config_directory):
        
    # this allows for the config to be created from a preexisting file
    config_run_file_path = config_directory + '/run_param.yml'
    if os.path.exists(config_run_file_path):
        with open(config_run_file_path, 'r') as cfgpath:
            run_config_dict = yaml.safe_load(cfgpath)
            print(run_config_dict)
    else:
        print('the path does not exist, check the path you gave VegET()')
        sys.exit(1)

    config_model_file_path = config_directory + '/model_param.yml'
    if os.path.exists(config_model_file_path):
        with open(config_model_file_path, 'r') as cfgpath:
            model_config_dict = yaml.safe_load(cfgpath)
            print(model_config_dict)
    else:
        print('the path does not exist, check the path you gave VegET()')
        sys.exit(1)

    config_dict = {**run_config_dict, **model_config_dict}
    print(config_dict)

    return config_dict
