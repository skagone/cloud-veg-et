VegET Tutorial
===================

Blah blah

Datasets
---------------------

The datasets, etc.

Gridmet Precipitation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODIS NDVI
~~~~~~~~~~~~~

Gridmet ETo, MACA ETo
~~~~~~~~~~~~~~~~~~~~~~~~~

NCRS Soils
~~~~~~~~~~~

Creating a Configuration File
------------------------------

Non-Dynamic Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most common application of the VegET model involves selecting time-series of gridded data that have overlapping
temporal coverage, and running the soil water balance using the same datasources of gridded input throughout the run. If
the user wishes to change datasources mid-stream (e.g. switching from Gridmet ETo to NLDAS ETo after 2005), please refer
to the section below on creating dynamic configuration files.

The standard configuration file is created in the YAML (YAML Ain't Markup Language) format...

Dynamic Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some situations, you may want to run VegET
soil-water-balance with different datasources during different time periods.
A good example is setting up a VegET run based on historical gridded dataset
availability, and transitioning datasets later in the run to model soil-water-balance
under future climate regimes.

The primary configuration files remain the same in this case, with the addtion of a dynamic configuration file...


Case Study: Setting up a Model Run for the Mesilla Valley.
===========================================================

Using many of the datasets above, let us demonstrate how a user would run the


Case Study II: Setting up a Model Run for a Large Study Area, using future climate scenarios.
==================================================================================================

