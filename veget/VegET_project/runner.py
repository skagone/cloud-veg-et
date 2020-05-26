from local_testing.VegET_model import VegConfig, PathManager, RasterManager, VegET
# from veget_model.car import Car

config_path = r'C:\WaterSmart\Projects\CloudVegET\local_testing\attributes.yml'

veggie = VegET(veget_config_path=config_path)


# run thyself!!!
veggie.run_veg_et()



