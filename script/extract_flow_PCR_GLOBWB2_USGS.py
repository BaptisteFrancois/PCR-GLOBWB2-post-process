
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import fiona
from shapely import geometry
from PIL import Image
from PIL.ExifTags import TAGS
from osgeo import gdal

import matplotlib.pyplot as plt

# Conversion factor from m^3/s to mm/day (x is the area of the basin in km^2)
m3s_to_mmday = lambda x: 86400 * 1000 / (x * 10**6)
m3s_to_cfs = 35.3147

# Path to the PCR-GLOBWB2 discharge output files
PCR_GLOBWB2_flow = '../data/PCR_GLOBWB2/discharge_monthAvg_output_1958-01-31_to_2015-12-31_zip.nc'

file2read = Dataset(PCR_GLOBWB2_flow, 'r')

# lat, lon = ... latitude and longitude of the outlet of the basins
latitude = file2read.variables['latitude'][:]
longitude = file2read.variables['longitude'][:]
time = file2read.variables['time'][:]
lon, lat = np.meshgrid(longitude, latitude)

# List of basin from the USGS

gages = pd.read_csv('../data/USGS/MetaData_USGS.csv', 
    dtype={'iD':str, 'area_km2':float, 'lat':float, 'lon':float, 'num_of_years':float}, 
    index_col='iD')


discharge = pd.DataFrame()

for i, gage in enumerate(gages.index):
    
    # Calculate the distance between the point and all grid cells
    distances = np.sqrt(
        (lat - gages.loc[gage, 'lat'])**2 
        + (lon - gages.loc[gage, 'lon'])**2)
    
    # Latitude and longitude of the closest PCR_GLOBWB2 grid cell
    id_closest = np.unravel_index(np.argmin(distances), distances.shape)
    
    # PCR_GLOBWB2 discharge in the closest grid cell (cfs)
    discharge = discharge.copy() # to get a defragmented dataframe
    discharge[gage] = file2read.variables['discharge'][:, id_closest[0],id_closest[1]] \
        * m3s_to_cfs
        
    print(f'{i+1}/{len(gages)}: {gage} done!')


dates = pd.date_range(start='1958-01-31', end='2015-12-31', freq='ME')
discharge.index = dates
discharge.index.name = 'Dates'
discharge.to_csv('../data/PCR_GLOBWB2/discharge_monthAvg_USGSgages_1958-01-31_to_2015-12-31.csv')


