
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

# Path to the PCR-GLOBWB2 discharge output files
PCR_GLOBWB2_flow = '../data/PCR_GLOBWB2/discharge_monthAvg_output_1958-01-31_to_2015-12-31_zip.nc'

file2read = Dataset(PCR_GLOBWB2_flow, 'r')

# lat, lon = ... latitude and longitude of the outlet of the basins
latitude = file2read.variables['latitude'][:]
longitude = file2read.variables['longitude'][:]
time = file2read.variables['time'][:]
lon, lat = np.meshgrid(longitude, latitude)

# List of basin IDs used in Maharjan et al. 2024
gages = np.loadtxt('../data/DPL-caravan_Maharjan_et_al_2024/GLOBAL_ALL_exclmissQ_area100to2k.txt',
                    dtype=str)

discharge = pd.DataFrame()

for i, gage in enumerate(gages):
    
    # Read the lat and lon of the outlet of the basin

    # The gage is located in United States
    if 'camels_' in gage:

        attributes = pd.read_csv(
            '../data/caravan/attributes/camels/attributes_other_camels.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])
        
    # The gage is located in Australia
    elif 'camelsaus_' in gage:
        
        attributes = pd.read_csv(
            '../data/caravan/attributes/camelsaus/attributes_other_camelsaus.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])
        
    # The gage is located in Brazil
    elif 'camelsbr_' in gage:

        attributes = pd.read_csv(
            '../data/caravan/attributes/camelsbr/attributes_other_camelsbr.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])
    
    # The gage is located in Chile
    elif 'camelscl_' in gage:

        attributes = pd.read_csv(
            '../data/caravan/attributes/camelscl/attributes_other_camelscl.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])
    
    # The gage is located in the United Kingdom
    elif 'camelsgb_' in gage:

        attributes = pd.read_csv(
            '../data/caravan/attributes/camelsgb/attributes_other_camelsgb.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])

    # The gage is located in the United States or Canada
    elif 'hysets_' in gage:

        attributes = pd.read_csv(
            '../data/caravan/attributes/hysets/attributes_other_hysets.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])

    # The gage is located in Central Europe
    elif 'lamah' in gage:

        attributes = pd.read_csv(
            '../data/caravan/attributes/lamah/attributes_other_lamah.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])
        
    # The gage is located in the Alps
    elif 'camelsch' in gage:
        
        attributes = pd.read_csv(
            '../data/caravan/attributes/camelsch/attributes_other_camelsch.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])
    
    
    # Calculate the distance between the point and all grid cells
    distances = np.sqrt(
        (lat - attributes.loc[gage, 'gauge_lat'])**2 
        + (lon - attributes.loc[gage, 'gauge_lon'])**2)
    
    # Latitude and longitude of the closest PCR_GLOBWB2 grid cell
    id_closest = np.unravel_index(np.argmin(distances), distances.shape)
    
    # PCR_GLOBWB2 discharge in the closest grid cell (mm/day)
    discharge = discharge.copy() # to get a defragmented dataframe
    discharge[gage] = file2read.variables['discharge'][:, id_closest[0],id_closest[1]] \
        * m3s_to_mmday(attributes.loc[gage, 'area'])
        
    print(f'{i+1}/{len(gages)}: {gage} done!')


dates = pd.date_range(start='1958-01-31', end='2015-12-31', freq='ME')
discharge.index = dates
discharge.index.name = 'Dates'
discharge.to_csv('../data/PCR_GLOBWB2/discharge_monthAvg_output_1958-01-31_to_2015-12-31.csv')



