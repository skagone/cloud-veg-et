import os
from s3fs.core import S3FileSystem as s3
import boto3
import sys

class PathManager:
    """
    This class creates the input paths for the dynamic and static data needed in the model.
    """


    # configuration non_std_inputs = True, accept parameters such as ndvi_fmt = userndvi_{}_vers1.tif
    # where {} is the date and date_fmt YYYYmmdd or YYYYdd etc for a filepath. A dictionary can relate user-input
    # name formatting to model std formats. Right now, non_std_inputs must be false and only 1 type of naming
    # convention per parameter file is allowed.

    ndvif = None
    pptf = None
    petf = None
    tavgf = None
    tminf = None
    tmaxf = None

    def __init__(self, config):
            self.config = config


    def get_dynamic_data(self, today, settings): # DOY=None, year_doy=None
        """
        This gets dynamic data, such as NDVI, precipitation, temperature, etc..
        :param today: datetime object to retrieve year, month, day, etc. from today
        :param settings: data set characteristics from the configurations file
        :return:
        """

        name_key = 'name_fmt'
        loc_key = 'dir_loc'
        dt_key = 'dt_fmt'
        clim_key = 'climatology'
        doy = today.timetuple().tm_yday

        print('settings', settings)

        if settings[clim_key]:
            # for climatology then we expect a DOY format
            if settings[dt_key] == 'doy':
                dynamic_key = '{:03d}'.format(doy)
            else:
                print('{} is set to climatology but date format from config is {}'.format(settings[name_key],
                                                                                          settings[dt_key]))
                sys.exit(0)
        elif settings[dt_key] == 'YYYYdoy':
            dynamic_key = '{}{:03d}'.format(today.year, doy)
        else:
            print('Hey user, the format of the dt_fmt configuration you gave: {} is not supported at '
                  'this time'.format(settings[dt_key]))
            sys.exit(0)

        fpath = os.path.join(settings[loc_key], settings[name_key].format(dynamic_key))
        if self.config.path_mode == 'local':
            pass
        elif self.config.path_mode == 'aws':
            fpath = s3.open(fpath)
        elif self.config.path_mode == 'google':
            # TODO
            print('google cloud bucket is not implemented. Returning as if the path were for a local')
            pass
        return fpath


    def get_static_data(self, settings):
        """
        This gets static data, such as soil data sets.
        :param settings: data set characteristics from the configurations file
        :return:
        """
        fpath = settings['file_loc']

        if self.config.path_mode == 'local':
            pass
        elif self.config.path_mode == 'aws':
            fpath = s3.open(fpath)
        elif self.config.path_mode == 'google':
            # TODO
            print('google cloud bucket is not implemented. Returning as if the path were for a local')
            pass
        return fpath


    def make_folder(self, folder_path):
        """"""
        # TODO - if user wants to put outputs in dir/subdir, how to handle with boto3???
        if self.config.path_mode == 'local':
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
        elif self.config.path_mode == 'aws':
            pass
            # if not os.path.exists(folder_path):
                # bucket= os.path.split(folder_path)[0]
                # KeyFileName = os.path.split(folder_path)[1]
                # s3 = boto3.client("s3")
                # with open(myfilename) as f:
                #     s3.
                #     s3.client.put_object(Bucket=bucket, Key=KeyFileName)





