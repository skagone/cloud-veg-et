class RasterManager:
    """
    This class addresses all the data wrangling needs to get the data sets needed in the model
    to the extent nad resolution chosen by the user.
    The user needs to input a sample shapefile (polygon) and geotiff file that is used to
    get the attributes for processing extent and resolution, etc.
    """

    # TODO - get large extent tif files for testing the RasterManager warping functions
    # TODO - make overlapping of 2 dataset work
    # TODO - rasterio either clips based on a sample raster and shapefile, or based on user defined geo info
    # #TODO - https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html
    # Todo - test raster manager as a standalone module.

    out_root = None
    # --- geographic info for destination files ---
    crs = None
    cols = None
    rows = None
    xres = None
    yres = None
    # left geo coord 'e.g. Westernmost Longitude"
    left = None
    # top geo coord e.g highest Latitude extent
    top = None
    transform = [xres, 0.0, left, 0.0, yres, top]

    sample_tiff = None
    # can contain multiple features of interest
    shapefile = None

    temp_folder = None


    def __init__(self, config):
        self.config = config

        self.sample_tiff = config.sample_tiff
        self.shapefile = config.shapefile
        self.temp_folder = os.path.join(config.out_root, config.temp_folder)
        if not os.path.exists(self.temp_folder):
            os.makedirs(self.temp_folder)

        if self.sample_tiff == None or self.shapefile==None:
            print('Assuming the user entered values in the config for boundaries of the AOI not implemented at thsi time')
            sys.exit(0)

    # ----------- create output rasters -----------------
    def output_rasters(self, arr, outdir, outname):
        """
        This function creates geotiff files from the model output arrays.
        """

        outpath = os.path.join(outdir, outname)
        print('the outpath for file {} is {}'.format(outname, outpath))

        # get the geoinfo from sample tiff to output intermediate files
        ds = rasterio.open(self.sample_tiff)
        band1 = arr
        with rasterio.open(outpath, 'w', driver='GTiff', height=self.rows, width=self.cols,
                           count=1, dtype='float64', crs=self.crs, transform=self.transform) as wrast:
            wrast.write(band1, indexes=1)

        # TODO - Set an AWS Cloud flag in the config file to activate this function or not...
        # delete files created locally and put in bucket
        # PathManager.s3_delete_local(from_file, bucket, prefix_no_slash)


    def set_model_std_grid(self, feat=0):
        """Clips and crops a tiff to the extent of a feature in a shapefile
        :param feat: feat is  the feature id of the shapefile from like a GeoJSON)
        # https://rasterio.readthedocs.io/en/latest/topics/virtual-warping.html
        """
        print(self.shapefile)
        with fiona.open(self.shapefile, 'r') as shapefile:
            # todo - set up an error if user has shapefile with more than one feature.
            # shape = shapefile[0]['geometry']
            shapes = [feature["geometry"] for feature in shapefile]

            for feature in shapefile:
                # matching the FID of the given shapefile from a typical geoJSON (Not Ordered Dict nonsense)
                if feat == feature['id']:
                    shapes = [feature['geometry']]
        # print('This is the shape var:', shapes)

        with rasterio.open(self.sample_tiff, 'r') as src:
            out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
            out_meta = src.meta
        # once the image is cropped, the image metadata dictionary is updated with the cropped transform and bounds.
        out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform})

        self.crs = out_meta['crs']
        self.transform = out_meta['transform']
        self.left = self.transform[2]
        self.top = self .transform[5]
        self.cols = out_meta['width']
        self.rows = out_meta['height']
        self.xres = self.transform[0]
        self.yres = self.transform[4]
        # return out_meta

    def normalize_to_std_grid(self, inputs, resamplemethod = 'nearest'):
        """
        Uses rasterio virtual raster to standardize grids of different crs, resolution, boundaries based on  a shapefile geometry feature
        :param inputs: a list of (daily) raster input files for the water balance.
        :param outloc: output locations 'temp' for the virtual files
        :return: list of numpy arrays
        """
        outputs = []
        npy_outputs = []
        if resamplemethod == 'nearest':
            rs = Resampling.nearest
        else:
            print('only nearest neighbor resampling is supported at this time')
            sys.exit(0)

        for i, warpfile in enumerate(inputs):
            print('warpfile', warpfile)
            with rasterio.open(warpfile) as src:
                # create the virtual raster based on the standard rasterio attributes from the sample tiff and shapefile feature.
                with WarpedVRT(src, resampling=rs,
                               crs=self.crs,
                               transform=self.transform,
                               height=self.rows,
                               width=self.cols) as vrt:
                    data = vrt.read()
                    print(type(vrt))
                    # save the file as an enumerated tiff. reopen outside this loop with the outputs list
                    outwarp = os.path.join(self.temp_folder, 'temp_{}.tif'.format(i))
                    rio_shutil.copy(vrt, outwarp, driver='GTiff')
                    outputs.append(outwarp)

        # output each virtual file as a temporary .tif file in a temp folder somewhere in the outputs directory.
        # for each file in the temp directory read in the raster as a numpy array and return the list of numpy arrays
        # from this method for us in the rest of the code.
        for ow in outputs:
            with rasterio.open(ow, 'r') as src:
                arr = src.read(1)
                npy_outputs.append(arr)

        return npy_outputs


