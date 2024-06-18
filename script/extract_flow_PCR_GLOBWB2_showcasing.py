
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import fiona
from shapely import geometry
from PIL import Image
from PIL.ExifTags import TAGS
from osgeo import gdal
from itertools import product
import matplotlib.pyplot as plt

# Conversion factors
m3s_to_mmday_area_in_km2 = lambda x: 86400 * 1000 / (x *10**6)
m3s_to_mmday_area_in_mile2 = lambda x: 86400 * 1000 / (x * 2.58999 *10**6)
cfs_to_mmday_area_in_mile2 = lambda x: 0.0283168 * 86400 * 1000 / (x * 2.58999 * 10**6)

# Path to the PCR-GLOBWB2 discharge output files
PCR_GLOBWB2_flow = '../data/PCR_GLOBWB2/discharge_monthAvg_output_1958-01-31_to_2015-12-31_zip.nc'

file2read = Dataset(PCR_GLOBWB2_flow, 'r')

# lat, lon = ... latitude and longitude of the outlet of the basins
latitude = file2read.variables['latitude'][:]
longitude = file2read.variables['longitude'][:]
time = file2read.variables['time'][:]
lon, lat = np.meshgrid(longitude, latitude)

# List of basins we are showcasing for the 29-01-2024 meeting
# (name, latitude, longitude, area, area_unit)
gages = (
    ('Brazil', 'IGUATU - 3650545', -6.36, -39.3, 21770, 'km2'),
    ('Brazil', 'PEIXE GORDO - 3650549', -5.22, -38.19, 48200, 'km2'),
    ('Speyside', '97002', 58.51, -3.49, 412.8, 'km2'),
    ('RioGrande', 'USGS08330000', 35.08, -106.68, 14500, 'sqmile'),
    ('Pangani', '1dd1', -3.25, 37.25, 240, 'km2'),
    
)

discharge = pd.DataFrame()

for i, (gage_region, gage_id, gage_lat, gage_lon, gage_area, area_unit) in enumerate(gages):
    
    # Read the observed gage data
    if gage_region == 'Pangani':
        gage = pd.read_csv('../data/showcasing_case_studies/pangani_runoff_observed.csv',
                           index_col='Dates', parse_dates=True)
        gage = gage.resample('ME').mean()
        start = '1958-01-31'
        end = '2006-03-30'
        gage = gage.truncate(before=start, after=end)
        
        
    
    elif gage_region == 'RioGrande':
        gage = pd.read_csv('../data/showcasing_case_studies/08330000_RIO_GRANDE_AT_ALBUQUERQUE.csv',
                           index_col='Dates', parse_dates=True) * cfs_to_mmday_area_in_mile2(gage_area)
        gage = gage.resample('ME').mean()
        start = '1975-01-31'
        end = '2015-12-31'
        gage = gage.truncate(before=start, after=end)

    elif gage_region == 'Speyside':
        gage = pd.read_csv('../data/showcasing_case_studies/97002_gdf_formatted.csv', 
                           index_col='Dates', parse_dates=True) * m3s_to_mmday_area_in_km2(gage_area)
        gage = gage.resample('ME').mean()
        start = '1972-01-31'
        end = '2015-12-31'
        gage = gage.truncate(before=start, after=end)

    elif gage_region == 'Brazil':
        if gage_id == 'IGUATU - 3650545':
            gage = pd.read_csv('../data/showcasing_case_studies/3650645_month.csv', 
                               index_col='Dates', parse_dates=True) * m3s_to_mmday_area_in_km2(gage_area)
            gage = gage.resample('ME').mean()
            gage['ms-1'] = np.where(gage < 0, np.nan, gage)
            start = '1958-01-31'
            end = '1992-12-31'

        elif gage_id == 'PEIXE GORDO - 3650549':
            gage = pd.read_csv('../data/showcasing_case_studies/3650649_month.csv', 
                               index_col='Dates', parse_dates=True) * m3s_to_mmday_area_in_km2(gage_area)
            gage = gage.resample('ME').mean()
            gage['ms-1'] = np.where(gage < 0, np.nan, gage)
            start = '1961-01-31'
            end = '2015-12-31'
        gage = gage.truncate(before=start, after=end)

        
    if area_unit == 'km2':
        conversion_PCR = m3s_to_mmday_area_in_km2(gage_area)
    elif area_unit == 'sqmile':
        conversion_PCR = m3s_to_mmday_area_in_mile2(gage_area)
    

    # Calculate the distance between the point and all grid cells
    distances = np.sqrt( (lat - gage_lat)**2 + (lon - gage_lon)**2 )
    
    # Latitude and longitude of the closest PCR_GLOBWB2 grid cell
    id_closest = np.unravel_index(np.argmin(distances), distances.shape)
    
    # PCR_GLOBWB2 discharge in the closest grid cell (mm/day)
    discharge = discharge.copy() # to get a defragmented dataframe

    dx = (-1,0,1)
    dy = (-1,0,1)
    difference = 999999999999999999999999999.9
    for k, (dxx, dyy) in enumerate(product(dx, dy)):
        flow_cell = \
            file2read.variables['discharge'][:, id_closest[0]+dyy,id_closest[1]+dxx] * conversion_PCR
        flow_cell = pd.DataFrame(
            flow_cell.data, index=pd.date_range(start='1958-01-31', end='2015-12-31', freq='ME'))
        flow_cell = flow_cell.truncate(before=start, after=end)
        flow_cellAvg = flow_cell.mean()[0]
        if abs(flow_cellAvg - gage.mean().values[0]) < difference:
            difference = abs(flow_cellAvg - gage.mean().values[0])
            id_closest_final = (id_closest[0]+dyy, id_closest[1]+dxx)
        
    flow = file2read.variables['discharge'][:, id_closest_final[0],id_closest_final[1]] * conversion_PCR
    flow = pd.DataFrame(
        flow.data, index=pd.date_range(start='1958-01-31', end='2015-12-31', freq='ME'))
    flow = flow.truncate(before=start, after=end)
    #print(f'{i+1}/{len(gages)}: {gage} done!')
    
    #bias.append(np.mean(flow.values - gage.values))

    # Plot some figures
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    ax.plot(gage, label='observed')
    ax.plot(flow, label='PCR-GLOBWB 2.0')
    ax.grid()
    ax.set_ylabel('mm', fontsize=12)
    ax.set_title(f'{gage_region} -- Gage id:{gage_id}', fontsize=12)
    ax.legend(fontsize=12)
    plt.tight_layout()
    fig.savefig(f'../figures/showcasing/{gage_region}_{gage_id}_PCR-GLOBWB2.0.png')
    plt.show()

    fig,ax = plt.subplots(1, 1, figsize=(4,4))
    ax.plot(gage, flow, 'ko')
    ax.set_xlabel('Observed (mm)', fontsize=12)
    ax.set_ylabel('PCR-GLOBWB 2.0 (mm)', fontsize=12)
    ax.axline((0, 0), slope=1, color="r", linestyle='--')
    ax.grid()
    plt.tight_layout()
    fig.savefig(f'../figures/showcasing/{gage_region}_{gage_id}_PCR-GLOBWB2.0_scatter.png')
    plt.show()

    
#dates = pd.date_range(start='1958-01-31', end='2015-12-31', freq='ME')
#discharge.index = dates
#discharge.index.name = 'Dates'
#discharge.to_csv('../data/PCR_GLOBWB2/discharge_monthAvg_output_1958-01-31_to_2015-12-31.csv')



