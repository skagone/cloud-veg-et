import os
from s3fs.core import S3FileSystem as s3
import boto3
import yaml
import sys
import ast
from datetime import date, datetime, timedelta

from .log_logger import log_make_logger

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

    def __init__(self, config_dict):
            self.config_dict = config_dict
            self.log = log_make_logger('PATHMANAGER')

    def get_dynamic_settings(self, today, settings):
        """"""

        # ===
        file = 'dynamic_file'
        interval_key = 'interval_list'
        dynamic_key = 'dynamic_keys'
        # ===

        dynamic_keys = settings[dynamic_key]
        interval_lst = settings[interval_key]

        # check to make sure the intervals and dynamic keys are the same len
        if not len(dynamic_keys) == len(interval_lst):
            print(f'there should me the same number of'
                  f' {interval_key} as {dynamic_key} in the config')
            sys.exit(0)

        # TODO - enable date ranges that are more precise than annual...
        for i in range(len(interval_lst)):
            year_tuple = ast.literal_eval(interval_lst[i])
            dk = dynamic_keys[i]

            dt_early = datetime(year=int(year_tuple[0]), month=1, day=1)
            dt_late = datetime(year=int(year_tuple[1]), month=1, day=1)

            # if today is between the start and end date of the period, use it to get the settings
            # from the yaml file.
            if dt_early <= today < dt_late:
                
                with open(settings[file], 'r') as rfile:
                    dynamic_settings_dict = yaml.safe_load(rfile)
#                 print('dynamic settings dictionary \n', dynamic_settings_dict)
#                 print('of type: ', type(dynamic_settings_dict))
                new_settings = dynamic_settings_dict[dk]
                return new_settings
            # go on down the list to find a day that is within the interval
            else:
                pass




    def get_dynamic_data(self, today, settings): # DOY=None, year_doy=None
        """
        This gets dynamic data, such as NDVI, precipitation, temperature, etc..
        :param today: datetime object to retrieve year, month, day, etc. from today
        :param settings: data set characteristics from the configurations file
        :return:
        """

        # Dynamic settings are where the settings and data-sets change based on the day or year.
        # Useful for including fore-casting or back-casting data with historical
        # satellite (observational) data-sets within a single veg-et run.
        try:
            dynamic_settings = settings['dynamic_settings']
            if dynamic_settings == None:
                dynamic_settings = False
        except KeyError:
            dynamic_settings = False

        name_key = 'name_fmt'
        loc_key = 'dir_loc'
        dt_key = 'dt_fmt'
        clim_key = 'climatology'
        doy = today.timetuple().tm_yday

        # print('settings', settings)

        if dynamic_settings:
            # the settings for the dynamic data are modified based on the date.
            settings = self.get_dynamic_settings(today=today, settings=settings)

        # the date finding tool is flexible enough so the non-climatological
        # data doesnt have to be organized in year folders. It CAN be though.
        # final_path = None
        if self.config_dict['path_mode'] == 'local':
            # print('local is happening')
            if settings[clim_key]:
                'clim is happening'
                # for climatology then we expect a DOY format
                if settings[dt_key] == 'doy':
                    dynamic_key = '{:03d}'.format(doy)
                    # we assume that there are no subdirectories so the final path is the same as the loc key
                    final_path = settings[loc_key]
                elif settings[dt_key] == 'mmdd':
                    dynamic_key = '{:02d}{:02d}'.format(today.month, today.day)
                    # we assume that there are no subdirectories so the final path is the same as the loc key
                    final_path = settings[loc_key]
                else:
                    print('{} is set to climatology but date format from config_dict is {}'.format(settings[name_key],
                                                                                              settings[dt_key]))
                    sys.exit(0)

            elif not settings[clim_key]:
                print('non clim is happening')
                if settings[dt_key] == 'YYYYdoy':
                    dynamic_key = '{}{:03d}'.format(today.year, doy)
                    # Do a walk in case there are subdirectories with each file in there.
                    walk_obj = os.walk(settings[loc_key])
                    for path, subdir, files in walk_obj:
                        # print('nonclim path \n', path)
                        for file in files:
                            # todo - REALLY REALLY someday better error handling logging here please
                            # print('nonclim file name', file)
                            if file == settings[name_key].format(dynamic_key):
                                # print('nonclim final \n', path)
                                final_path = path
                elif settings[dt_key] == 'YYmmdd':
                    # get the last two digits of year
                    yr = str(today.year)
                    yr_2dig = yr[-2:]
                    dynamic_key = '{}{:02d}{:02d}'.format(yr_2dig, today.month, today.day)
                    # Do a walk in case there are subdirectories with each file in there.
                    walk_obj = os.walk(settings[loc_key])
                    for path, subdir, files in walk_obj:
                        # print('nonclim path \n', path)
                        durectory = os.path.split(path)[1]
                        if durectory == settings[name_key].format(dynamic_key):
                            final_path = os.path.split(path)[0]
                        else:
                            pass
                            for file in files:
                                # print('nonclim file name', file)
                                if file == settings[name_key].format(dynamic_key):
                                    # print('nonclim final \n', path)
                                    final_path = path
            else:
                print('Hey user, the format of the dt_fmt configuration you gave: {}, is not supported at '
                      'this time'.format(settings[dt_key]))
                sys.exit(0)

        elif self.config_dict['path_mode'] == 'aws':
            self.log.debug('PATH MODE is AWS - tony will be happy')
            # account for subdirs, called prefix in the cloud, if data is not climatology
            if settings[clim_key]:
                # for climatology then we expect a DOY format, and we assume no prefixes
                if settings[dt_key] == 'doy':
                    dynamic_key = '{:03d}'.format(doy)
                    final_path = settings[loc_key]
                else:
                    print('{} is set to climatology but date format from config_dict is {}'.format(settings[name_key],
                                                                                              settings[dt_key]))
                    sys.exit(0)
            # if its NOT a climatology add the year to the loc_dir (for now)
            elif settings[dt_key] == 'YYYYdoy':
                dynamic_key = '{}{:03d}'.format(today.year, doy)
                final_path = os.path.join(settings[loc_key], '{}'.format(today.year))
            else:
                print('Hey user, the format of the dt_fmt configuration you gave: {} is not supported at '
                      'this time'.format(settings[dt_key]))
                sys.exit(0)

        elif self.config_dict['path_mode'] == 'google':
            # TODO
            print('google cloud bucket is not implemented. Please use *aws* or *local* input args to run the model!!!')
            sys.exit(0)

        # print(self.config_dict['path_mode'])
        # print(settings[name_key].format(dynamic_key))
        fpath = os.path.join(final_path, settings[name_key].format(dynamic_key))
        return fpath


    def get_static_data(self, settings):
        """
        This gets static data, such as soil data sets.
        :param settings: data set characteristics from the configurations file
        :return:
        """
        fpath = settings['file_loc']

        if self.config_dict['path_mode'] == 'local':
            pass
        elif self.config_dict['path_mode'] == 'aws':
            # fpath = s3.open(fpath)
            pass
        elif self.config_dict['path_mode'] == 'google':
            # TODO
            print('google cloud bucket is not implemented. Returning as if the path were for a local')
            pass
        return fpath

    def make_s3_output_path(self):
        print(self.config_dict['out_root_prefix'])
        print(self.config_dict['region'])
        region = self.config_dict['region']
        today = date.today()
        print("Current date =", today)
        date_str=today.strftime("%m_%d_%Y")
        s3_output_path = self.config_dict['out_root_prefix'] + '/' + region + '/Run' + date_str + '/' + self.config_dict['tile']
        return(s3_output_path)


    def make_folder(self, folder_path):
        """"""
        # TODO - if user wants to put outputs in dir/subdir, how to handle with boto3???
        if self.config_dict['path_mode'] == 'local':
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
        elif self.config_dict['path_mode'] == 'aws':
            pass
            # if not os.path.exists(folder_path):
                # bucket= os.path.split(folder_path)[0]
                # KeyFileName = os.path.split(folder_path)[1]
                # s3 = boto3.client("s3")
                # with open(myfilename) as f:
                #     s3.
                #     s3.client.put_object(Bucket=bucket, Key=KeyFileName)





