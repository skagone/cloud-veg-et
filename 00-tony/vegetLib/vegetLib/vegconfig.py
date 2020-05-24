class VegConfig:
    """"""

    # Attributes Here.

    # === Minimalist Beta Version Params ====
    start_year = None
    end_year = None
    start_day = None
    end_day = None
    sample_tiff = None
    shapefile = None

    # ==== Version 2.0 params ====
    #input/output
    in_root = None
    out_root = None

    precip_settings = {'crs': None, 'cols': None, 'rows': None, 'xres': None, 'yres': None, 'left': None, 'top': None}
    ndvi_settings = {'crs': None, 'cols': None, 'rows': None, 'xres': None, 'yres': None, 'left': None, 'top': None}
    pet_settings = {'crs': None, 'cols': None, 'rows': None, 'xres': None, 'yres': None, 'left': None, 'top': None}
    tmin_settings = {'crs': None, 'cols': None, 'rows': None, 'xres': None, 'yres': None, 'left': None, 'top': None}
    tavg_settings = {'crs': None, 'cols': None, 'rows': None, 'xres': None, 'yres': None, 'left': None, 'top': None}
    tmax_settings = {'crs': None, 'cols': None, 'rows': None, 'xres': None, 'yres': None, 'left': None, 'top': None}

    non_std_inputs = False


    def __init__(self, config_dictionary):
        for key, val in config_dictionary.items():
            setattr(self, key, val)

    def update_config(self, k, v, cfg_path):
        """"""
        setattr(self, k, v)
        with open(cfg_path, 'r') as cfg:
            file_dict = yaml.safe_load(cfg)
        file_dict['{}'.format(k)] = v
        with open(cfg_path, 'w') as w_file:
            w_file.write(yaml.dump(file_dict))
            # # old way
            # append_file.write('{}: {}\n'.format(k, v))
    def update_feature(selfs, k, v, cfg_path):
        """Used to update the config file when a configuration needs to be ADDED or APPENDED to something, like a
        dictionary that gets new entries
        This function does not update the attribute like update_config does. This function is used when an attribute
         has ALREADY been appended or updated in the code but we still need to update the configuration FILE"""
        with open(cfg_path, 'r') as cfg:
            file_dict = yaml.safe_load(cfg)
        # overprint the entries with the new config
        file_dict['{}'.format(k)] = v
        with open(cfg_path, 'w') as w_file:
            w_file.write(yaml.dump(file_dict))


