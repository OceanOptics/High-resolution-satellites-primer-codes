## Case Study 3.2 Provisional Aquatic Reflectance from USGS: Coastal to Gulf Stream Rrs and Chl-a (chla)

This jupyter notebook processes the data and all it needs as input is file locations and spatial bounding boxes.

To run this notebook we recommend using Docker which can easily be pulled and run via the command

docker run -it -v <code dir>:/home/jovyan --rm -p 8888:8888 pangeo/pangeo-notebook:2021.05.15 jupyter lab --ip 0.0.0.0

And just replace <code dir> with the directory you would like the notebook to access.

If you want to build the environment yourself it is fairly basic and all you need are:
xarray matplotlib numpy unpackqa xrft rioxarray cartopy
