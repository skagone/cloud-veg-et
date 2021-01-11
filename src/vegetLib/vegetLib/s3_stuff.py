import os
import boto3

def s3_hello(person_name):
    print('Hello There Person:', person_name)

def s3_push_delete_local(local_file, bucket, bucket_filepath):
        s3 = boto3.client('s3')
        with open(local_file, "rb") as f:
            if 'vsis3' in bucket:
                bucket = bucket.split('/')[-1]
                print(bucket, bucket_filepath)
            s3.upload_fileobj(f, bucket, bucket_filepath)
        os.remove(local_file)

def return_s3_list(working_bucket, prefix):
        aws_list = []
        s3 = boto3.resource('s3')
        bucket_name = working_bucket
        bucket = s3.Bucket(bucket_name)
        for obj in bucket.objects.filter(Prefix=prefix):
            obj_key = obj.key
            obj_key = working_bucket + '/' + obj_key
            aws_list.append((obj_key, obj.size))
        return aws_list


def s3_list_pseudo_subdirs(bucket, prefix_with_slash):

    subfolder_list = []
    #Make sure you provide / in the end
    prefix = prefix_with_slash 

    client = boto3.client('s3')
    result = client.list_objects(Bucket=bucket, Prefix=prefix, Delimiter='/')
    for o in result.get('CommonPrefixes'):
        #print ('sub folder : ', o.get('Prefix'))
        subfolder_list.append(o.get('Prefix'))
    return subfolder_list




