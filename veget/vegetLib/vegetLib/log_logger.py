import logging
import boto3

def log_make_logger(nameV):
   
    LOGGER = logging.getLogger(nameV)

    LOGGER.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
   
    myStreamTypicallySTDOUT = logging.StreamHandler()
    myStreamTypicallySTDOUT.setFormatter(formatter)
    LOGGER.addHandler(myStreamTypicallySTDOUT)

    fh = logging.FileHandler('./log/run.log')
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
        local_file = './log/run.log'
        with open(local_file, "rb") as f:
            bucket = s3_output_path.split('/')[0]
            prefix = '/'.join(s3_output_path.split('/')[1:])
            bucket_filepath = prefix + '/aaalog/run.log'
            print(bucket, bucket_filepath)
            s3.upload_fileobj(f, bucket, bucket_filepath)
