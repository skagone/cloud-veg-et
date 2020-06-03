from veget.vegetLib.vegetLib.veget import VegET
import os
import rasterio
from matplotlib import pyplot as plt


veggie = VegET(veget_config_path=r'D:\Users\gparrish\PycharmProjects\cloud-veg-et\veget\api_veget\sample_config', tile='40N-80E',
               shp=r'Z:\Projects\Cloud_Veg_ET\Data\shapefiles\DRB.shp')

veggie.run_veg_et()