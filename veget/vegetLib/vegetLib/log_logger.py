import logging
import os
import boto3

def log_make_logger(nameV):
   
    LOGGER = logging.getLogger(nameV)

    LOGGER.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
   
    myStreamTypicallySTDOUT = logging.StreamHandler()
    myStreamTypicallySTDOUT.setFormatter(formatter)
    LOGGER.addHandler(myStreamTypicallySTDOUT)

    log_file_name = './log/' + nameV + '.log'

    fh = logging.FileHandler(log_file_name, mode='w')
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    LOGGER.addHandler(fh)
    
    LOGGER.info("Logging Begins")
    return LOGGER

from inspect import currentframe

def log_get_line_number():
    cf = currentframe()
    return cf.f_back.f_lineno


def s3_save_log_file(s3_output_path):

        s3 = boto3.client('s3')
        log_files = os.listdir("./log")
        for lf in log_files:
            if '.log' in lf:
                local_file = './log/' + lf
                with open(local_file, "rb") as f:
                    bucket = s3_output_path.split('/')[0]
                    prefix = '/'.join(s3_output_path.split('/')[1:])
                    bucket_filepath = prefix + '/aaalog/' + lf
                    print(bucket, bucket_filepath)
                    s3.upload_fileobj(f, bucket, bucket_filepath)
                f.close()
