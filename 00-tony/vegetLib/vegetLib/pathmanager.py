class PathManager:
    """
    This class creates the input paths for the dynamic and static data needed in the model.
    """


    # Todo - configuration non_std_inputs = True, accept parameters such as ndvi_fmt = userndvi_{}_vers1.tif
    #  where {} is the date and date_fmt YYYYmmdd or YYYYdd etc for a filepath. A dictionary can relate user-input
    #  name formatting to model std formats. Right now, non_std_inputs must be false and only 1 type of naming
    #  convention per paramter file is allowed.

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
        return fpath


    def get_static_data(self, settings):
        """
        This gets static data, such as soil data sets.
        :param settings: data set characteristics from the configurations file
        :return:
        """
        fpath = settings['file_loc']
        return fpath


    def s3_delete_local(self, outpath, from_file, bucket, prefix_no_slash):
        """
        This function will move the model outputs from a local folder to a cloud bucket.
        :param from_file: path the the local folder
        :param bucket: name of the cloud bucket = 'dev-et-data'
        :param prefix_no_slash: "folder" in cloud bucket  = 'v1DRB_outputs'
        :return:
        """

        objecta='{}/{}'.format(prefix_no_slash,outpath)
        s3 = boto3.client('s3')
        with open(from_file, "rb") as f:
            s3.upload_fileobj(f, bucket, objecta)
        os.remove(from_file)


