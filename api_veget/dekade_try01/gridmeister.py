import os
import json

def _write_geojson(filename, coord_list):
    print(coord_list)
    polycs = []
    polycs.append(coord_list)
    geos = []
    for polyc in polycs:
        poly = {
            'type': 'Feature',
            'properties': {},
            'geometry': {
                'type': 'Polygon',
                'coordinates': [polyc]
            }
        }
        geos.append(poly)

    geometries = {
       'type': 'FeatureCollection',
        'features': geos,
    }

    with open(filename, 'w', encoding='utf-8') as f:
        json.dump(geometries, f, ensure_ascii=False, indent=4)
    print('WROTE filename {}'.format(filename))

def _write_shp(geojson_filename):
    shp_filename = geojson_filename.split('.json')[0] + '.shp'
    print(geojson_filename, shp_filename)
    cmd = 'ogr2ogr -f \"ESRI Shapefile\" {} {}'.format(shp_filename, geojson_filename)
    os.system(cmd)

def _make_extent_poly(extent):
    # (left(lon), bottom(lat) : right(lon), top(lat)) ->lon=x, lat=y
    ul_lon_lat = [extent[0], extent[3]]
    ur_lon_lat = [extent[2], extent[3]]
    lr_lon_lat = [extent[2], extent[1]]
    ll_lon_lat = [extent[0], extent[1]]
    lst = [ul_lon_lat, ur_lon_lat, lr_lon_lat, ll_lon_lat, ul_lon_lat]

    return lst

def _make_chip_poly(ul_lat, ul_lon, x_increment, y_increment):
    coord_list = []
    ul_lon_lat = [ul_lon, ul_lat]
    ur_lon_lat = [ul_lon + x_increment, ul_lat]
    lr_lon_lat = [ul_lon + x_increment, ul_lat - y_increment]
    ll_lon_lat = [ul_lon, ul_lat - y_increment]
    # print(ul_lon_lat)
    coord_list.append(ul_lon_lat)
    # print(ur_lon_lat)
    coord_list.append(ur_lon_lat)
    # print(lr_lon_lat)
    coord_list.append(lr_lon_lat)
    # print(ll_lon_lat)
    coord_list.append(ll_lon_lat)
    # print(ul_lon_lat)
    coord_list.append(ul_lon_lat)
    return coord_list

def make_filename(tile_name, chip_name, extension):
    filename = tile_name + '_' + chip_name + extension
    return(filename)

class GridMeister:
    """
    This class partitions a 10 degree tile into 2 degree chips
    """

    def __init__(self, tile_name, raster_extent, x_raster_res, y_raster_res):

        self.tile_name = tile_name
        self.chip_increment = 2
        self.extent = raster_extent
        self.width = x_raster_res
        self.height = y_raster_res
        self.xchip_increment = None
        self.ychip_increment = None


    def chip_list(self, max_pixels=500000):
        CHIP_LIST = []
        # box={'left': -78, 'bottom':36 , 'right': -72, 'top': 44}
        box = {'left': self.extent[0], 'bottom': self.extent[1],
               'right': self.extent[2], 'top': self.extent[3]}

        starting_lat = box['top']
        ending_lat = box['bottom']
        # print('lats', starting_lat, ending_lat)

        starting_lon = box['left']
        ending_lon = box['right']
        # print('lons', starting_lon, ending_lon)

        lat_cols = int(abs((starting_lat - ending_lat) / self.width))
        lon_rows = int(abs((starting_lon - ending_lon) / self.height))

        chip_pixels = lat_cols * lon_rows
        count = 1
        # print(chip_pixels, max_pixels)

        while chip_pixels > max_pixels:
            count += 1
            lat_cols = int(abs(((starting_lat - ending_lat)/count) / self.width))
            lon_rows = int(abs(((starting_lon - ending_lon)/count) / self.height))
            self.ychip_increment = abs(((starting_lat - ending_lat)/count))
            self.xchip_increment = abs(((starting_lon - ending_lon)/count))
            chip_pixels = lat_cols * lon_rows
            # print(chip_pixels, max_pixels)

        # print(self.xchip_increment, self.ychip_increment, lat_cols, lon_rows)
        # ===================================================================
        if not self.xchip_increment==None:
            lat = starting_lat
            while lat > ending_lat + 0.00000005: #-self.ychip_increment
                # print(lat)
                lon = starting_lon
                while lon < ending_lon:
                    chip = (lat, lon)
                    CHIP_LIST.append(chip)
                    lon = lon + self.xchip_increment
                lat = lat - self.ychip_increment
        else:
            CHIP_LIST = None

        return CHIP_LIST

    def create_chip_shp(self, ul_lat, ul_lon, out_location, unit_chip=False):
        print(ul_lat, ul_lon)
        if not unit_chip:
            coord_list = _make_chip_poly(ul_lat, ul_lon, self.xchip_increment, self.ychip_increment)
            print(coord_list)
            chip_name = 'chip' + str(round(ul_lat, 2)) + 'N' + str(round(ul_lon, 2)) + 'E'
            filename = '{}.json'.format(chip_name)
            print(filename)
        else:
            coord_list = _make_extent_poly(extent=self.extent)
            chip_name = 'chip' + self.tile_name
            filename = '{}.json'.format(chip_name)

        full_filename = os.path.join(out_location, filename)
        print('the fulll FILENAME: ', full_filename)
        _write_geojson(full_filename, coord_list)
        _write_shp(full_filename)

    def build_docker_run_bash(self, chip_list, optimize):
        print(chip_list)
        vols = '-v `pwd`/AOI:/home/veget/cloud-veg-et/api_veget/AOI'
        mycwd = os.getcwd()
        image_custom = mycwd.split('/')[-1]
        image = 'tbutzer/' + image_custom
        cmds=[]
        for chip_name in chip_list:
            filename = make_filename(self.tile_name, chip_name, '.shp')
            full_filename = self.aoi_dir + '/' + filename
            tile = filename.split('.shp')[0]
            cmd = 'docker run -i {} {} python3 api_veget.py -c running_config -s {}  {}'.format(vols,image,full_filename,tile)
            #cmd = 'docker run -i {} {} python3 bench_api_veget.py -c running_config -s {}  {}'.format(vols,image,full_filename,tile)
            if not optimize:
                cmd = 'docker run -i {} {} python3 api_veget.py -c running_config -s {} --optimize no  {}'.format(vols,image,full_filename,tile)
            print(cmd)
            logname = make_filename(self.tile_name, chip_name, '.log')
            # TODO -
            full_logname = './log' + '/' + logname
            print(full_logname)

            full_cmd = cmd + '  2>&1 | tee  ' + full_logname +'&'
            print(full_cmd)
            cmds.append(full_cmd)

        cmd_filename = 'cmd_runner_' + self.tile_name + '.sh'
        with open(cmd_filename, 'w') as f:
            for cmd in cmds:
                print(cmd)
                f.write(cmd+'\n')

        cmd_filename = 'test_cmd_runner_' + 'one' + '.sh'
        with open(cmd_filename, 'w') as f:
            cmd = cmds[0]
            print(cmd)
            f.write(cmd+'\n')

if __name__ == "__main__":

    print('======================')
    print('TESTING ZE GRIDMEISTER')
    print('======================')

    chip_output = r'D:\Users\gparrish\Desktop\gridmeistertest_Darin'
    # (left(lon), bottom(lat) : right(lon), top(lat) lon=x, lat=y
    exp_extent = (-76.7937822733839965, 38.0837012906591070,
                  -73.7104960956225455, 43.0419587927349596)
    darin_extent = (-77.4226093565661415, 38.3565890497327118, -73.2480170973858122, 42.7829967799980722)
    # (left, bottom: right, top)
    # tony_extent = (-78, 36, -72, 44)
    gm = GridMeister(tile_name='testtile', raster_extent=darin_extent,
                     x_raster_res=0.002310233679679207525, y_raster_res=0.002310233679679207525)
    # gm = GridMeister(tile_name='testtile', raster_extent=exp_extent,
    #                  x_raster_res=0.04166602942920882846, y_raster_res=0.04166602942920882846)
    #0.04166602942920882846
    lst = gm.chip_list()
    print('resulting chip list \n', lst)
    gm.create_chip_shp(ul_lat=None, ul_lon=None, out_location=chip_output, unit_chip=True)

    if lst == None:
        gm.create_chip_shp(ul_lat=None, ul_lon=None, out_location=chip_output, unit_chip=True)
    else:
        for l in lst:
            gm.create_chip_shp(l[0], l[-1], out_location=chip_output)