Introduction
==============
The VegET model is a distributed-parameter soil-water-balance model that uses available energy
and precipitation to model ET based on plant water-stress and remotely sensed phenology.

Soil Water Balance Modeling
=============================
Stefi

The VegET Model
=================

Model Algorithms
~~~~~~~~~~~~~~~~~~
Stefi

Model Inputs
~~~~~~~~~~~~~

Stefi/GELP

All water and energy-balance components of the VegET model
are input in the form of gridded datasets (GeoTiffs). The GeoTiff inputs
are stored either locally in directories, or within the cloud (e.g. AWS).

Precipitation
-----------------

Normalized Difference Vegetation Index (NDVI)
----------------------------------------------

Reference Evapotranspiration
-------------------------------

Soil Grids
------------

VegET Class
~~~~~~~~~~~~
Stefi/GELP

PathManager Class
~~~~~~~~~~~~~~~~~~~
GELP

The PathManager Class is a helper-class that handles the
management and tracking of model inputs at each time-step.

RasterManager Class
~~~~~~~~~~~~~~~~~~~~~~
GELP

The RasterManager Class is a helper-class that employs
open-source Python packages such as Rasterio to integrate
and standardize the gridded datasets that make up the inputs to
the VegET model.

Cloud Computing
~~~~~~~~~~~~~~~~~~~

Tony Butzer
