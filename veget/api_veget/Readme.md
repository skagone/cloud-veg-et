# API for veget

## The main function caller

## Synopsys


delaware:
	python3 api_veget.py -c delaware_config -s DRB.shp lilDRB_temp
    

## Errors

### Error 101

warpfile /vsis3/dev-et-data/NA_data_for_cloud/ETo_mosaic/pet366.tif
rasterio._err.CPLE_OpenFailedError: '/vsis3/dev-et-data/NA_data_for_cloud/ETo_mosaic/pet366.tif' does not exist in the file system, and is not recognized as a supported dataset name.

upyter-kagone@ip-10-12-68-72:~$ aws s3 ls dev-et-data/NA_data_for_cloud/ETo_mosaic/ | tail -4
2020-05-22 16:39:09       2733 pet364.tif.aux.xml
2020-05-22 16:39:09         83 pet365.tfw
2020-05-22 16:39:10 2814442690 pet365.tif
2020-05-22 16:39:50       2729 pet365.tif.aux.xml

### Workaround

- adjust run param to be day 364 - ask Steffi about pet366?