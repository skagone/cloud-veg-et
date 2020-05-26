import fiona
import rasterio
import rasterio.mask

def test():
    with fiona.open("in_box.shp", "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]

    print('This is the shape var:' + shapes)