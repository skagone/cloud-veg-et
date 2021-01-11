# API for veget

## The main function caller

## Synopsys


delaware:
	python3 api_veget.py -c delaware_config -s DRB.shp lilDRB_temp

## TODO

1. XX Cleanup debug messages
2. XX add timings to input warp warp_one function
3. run in docker container and capture logs
4. analyze timings for average warp per file type
5. fix ETo pet366 problem

## Optimizations

1. use fsspec to write directly the outputs and time the difference between the local and s3copy and direct s3 writes
2. use zarray and zarr on the outputs
    

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


### Error 102

```
 start_dt = datetime.strptime("{}-{:03d}".format(self.start_year, self.start_day), '%Y-%j')
        print(start_dt)
        end_dt = datetime.strptime("{}-{:03d}".format(self.end_year, self.end_day), '%Y-%j')
        print(end_dt)
        time_interval = end_dt - start_dt
        num_days = time_interval.days
        print(num_days)

```

The above code only handles consecutive days

```
start_year: 2012
end_year: 2014
start_day: 335
end_day: 360
```

therefore the above does not work as advertised


### Git Stash Tutorial

https://www.youtube.com/watch?v=KLEDKgMmbBI

- git stash list
- git stash pop
