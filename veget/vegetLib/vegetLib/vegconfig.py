import os
import sys
import yaml
import boto3

def return_veget_params(config_directory):
        
    # this allows for the config_dict to be created from a preexisting file
    config_run_file_path = os.path.join(config_directory, 'run_param.yml')
    if os.path.exists(config_run_file_path):
        with open(config_run_file_path, 'r') as cfgpath:
            run_config_dict = yaml.safe_load(cfgpath)
            # print(run_config_dict)
    else:
        print('the path does not exist, check the path you gave VegET() {}'.format(config_directory))
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

def s3_save_config_files(config_directory, s3_output_path):
    s3 = boto3.client('s3')
    file_list = ['model_param.yml',  'path_param.yml',  'run_param.yml']
    for file in file_list:
        local_file = config_directory + '/' + file
        with open(local_file, "rb") as f:
            bucket = s3_output_path.split('/')[0]
            prefix = '/'.join(s3_output_path.split('/')[1:])
            bucket_filepath = prefix + '/aaalog/' + file
            print(bucket, bucket_filepath)
            s3.upload_fileobj(f, bucket, bucket_filepath)
        f.close()


