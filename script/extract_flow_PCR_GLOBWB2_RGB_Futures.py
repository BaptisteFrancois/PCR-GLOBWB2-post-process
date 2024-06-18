
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

def nash(predictions, targets):
    return 1-(np.sum((targets-predictions)**2)/np.sum((targets-np.mean(targets))**2))

# Conversion factors
m3s_to_mmday_area_in_km2 = lambda x: 86400 * 1000 / (x *10**6)
m3s_to_mmday_area_in_mile2 = lambda x: 86400 * 1000 / (x * 2.58999 *10**6)
cfs_to_mmday_area_in_mile2 = lambda x: 0.0283168 * 86400 * 1000 / (x * 2.58999 * 10**6)
m3s_to_afd = 70.0452
cfs_to_afd = 1.98347

# Path to the PCR-GLOBWB2 discharge output files
PCR_GLOBWB2_flow = '../data/PCR_GLOBWB2/discharge_monthAvg_output_1958-01-31_to_2015-12-31_zip.nc'

file2read = Dataset(PCR_GLOBWB2_flow, 'r')

# lat, lon = ... latitude and longitude of the outlet of the basins
latitude = file2read.variables['latitude'][:]
longitude = file2read.variables['longitude'][:]
time = file2read.variables['time'][:]
lon, lat = np.meshgrid(longitude, latitude)

# # Flow data from these files is in AF/day
gages = (
    ('RGB', 'USGS08353000', 34.41, -106.85, 'USGS_08353000.csv'),
    ('RGB', 'USGS08362500', 32.88, -107.29, '08362500_RIO_GRANDE_BLW_CABALLO_DAM.csv'),
    ('RGB', 'USGS08361000', 33.15, -107.21, '08361000_RIO_GRANDE_BELOW_ELEPHANT_BUTTE_DAM_fill.csv'),
    ('RGB', 'USGS08317400', 35.62, -106.32, 'USGS_08317400.csv')
)

discharge = pd.DataFrame()

statistics = pd.DataFrame(columns=['min', 'q5', 'q25', 'q50', 'q75', 'q95', 'max', 'mean', 'std', 'nash'])

for i, (gage_region, gage_id, gage_lat, gage_lon, name_file) in enumerate(gages):
    
        
    if gage_id == 'USGS08362500':
        # Common period for the observed and simulated data (PCR and RGB Futures)
        start = '2001-01-31'
        end = '2013-12-31'
        
        gage = pd.read_csv('../data/RGB_Futures/{}'.format(name_file),
            index_col='date', parse_dates=True, usecols=['date', 'gauges'])
        gage = gage.resample('ME').mean()
        gage = gage.truncate(before=start, after=end)

        RGB_Futures = pd.read_csv(
            '../data/RGB_Futures/rio_grande_output_file_02152023_baseline_cost.csv',
            index_col='date', parse_dates=True, usecols=['date', 'caballo_release_AF', 'caballo_spill_AF'])
        RGB_Futures = RGB_Futures.resample('ME').mean()
        RGB_Futures = RGB_Futures.truncate(before=start, after=end)
        

    elif gage_id == 'USGS08361000':
        
        # Common period for the observed and simulated data (PCR and RGB Futures)
        start = '2001-01-31'
        end = '2015-12-31'
        
        gage = pd.read_csv('../data/RGB_Futures/{}'.format(name_file),
            index_col='date', parse_dates=True, usecols=['date', 'gauges'])
        gage = gage.resample('ME').mean()
        gage = gage.truncate(before=start, after=end)

        RGB_Futures = pd.read_csv(
            '../data/RGB_Futures/rio_grande_output_file_02152023_baseline_cost.csv',
            index_col='date', parse_dates=True, usecols=['date', 'elephant_butte_release_AF', 'elephant_butte_spill_AF'])
        RGB_Futures = RGB_Futures.resample('ME').mean()
        RGB_Futures = RGB_Futures.truncate(before=start, after=end)

    elif gage_id == 'USGS08317400':

        # Common period for the observed and simulated data (PCR and RGB Futures)
        start = '2001-01-31'
        end = '2015-12-31'

        gage = pd.read_csv('../data/RGB_Futures/{}'.format(name_file),
            index_col='Dates', parse_dates=True, usecols=['Dates', 'Flow_cfs']) * cfs_to_afd
        gage = gage.resample('ME').mean()
        gage = gage.truncate(before=start, after=end)

        RGB_Futures = pd.read_csv(
            '../data/RGB_Futures/rio_grande_output_file_02152023_baseline_cost.csv',
            index_col='date', parse_dates=True, usecols=['date', 'cochiti_release_AF', 'cochiti_spill_AF'])
        RGB_Futures = RGB_Futures.resample('ME').mean()
        RGB_Futures = RGB_Futures.truncate(before=start, after=end)

    elif gage_id == 'USGS08353000':
        # Common period for the observed and simulated data (PCR and RGB Futures)
        start = '1958-11-01'
        end = '2015-12-31'

        gage = pd.read_csv('../data/RGB_Futures/{}'.format(name_file),
            index_col='Dates', usecols=['Dates', 'Flow_cfs']) * cfs_to_afd
        gage.index = pd.to_datetime(gage.index)
        gage = gage.resample('ME').mean()
        gage = gage.truncate(before=start, after=end)

            
   

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
            file2read.variables['discharge'][:, id_closest[0]+dyy,id_closest[1]+dxx] * m3s_to_afd
        flow_cell = pd.DataFrame(
            flow_cell.data, index=pd.date_range(start='1958-01-31', end='2015-12-31', freq='ME'))
        flow_cell = flow_cell.truncate(before=start, after=end)
        flow_cellAvg = flow_cell.mean()[0]
        print(dxx, dyy, flow_cellAvg, gage.mean().values[0])
        print('difference is {}'.format(abs(flow_cellAvg - gage.mean().values[0])))

        if abs(flow_cellAvg - gage.mean().values[0]) < difference:
            difference = abs(flow_cellAvg - gage.mean().values[0])
            id_closest_final = (id_closest[0]+dyy, id_closest[1]+dxx)
        
    #flow = file2read.variables['discharge'][:, id_closest_final[0],id_closest_final[1]] * m3s_to_afd
    
    flow = file2read.variables['discharge'][:, id_closest[0],id_closest[1]] * m3s_to_afd

    flow = pd.DataFrame(
        flow.data, index=pd.date_range(start='1958-01-31', end='2015-12-31', freq='ME'))
    flow = flow.truncate(before=start, after=end)
    flow.index.name = 'Date'
    flow.columns = ['PCR_GLOBWB2']

    
    flow.to_csv('../data/RGB_Futures/PCR_GLOBWB2_{}.csv'.format(gage_id))
    #print(f'{i+1}/{len(gages)}: {gage} done!')
    
    #bias.append(np.mean(flow.values - gage.values))

    # Plot some figures
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    ax.plot(gage, 'k', label='Observed')
    ax.plot(flow, 'b', label='Competitor')
    if gage_id != 'USGS08353000':
        ax.plot(RGB_Futures.sum(axis=1),'r', label='TOVA')
    #ax.grid()
    ax.set_ylabel('AF/day', fontsize=12)
    ax.set_title(f'{gage_region} -- Gage id:{gage_id}', fontsize=12)
    ax.legend(fontsize=12)
    plt.tight_layout()
    fig.savefig(f'../figures/showcasing/{gage_region}_{gage_id}_PCR-GLOBWB2.0.png')
    #plt.show()
    plt.close()

    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.plot(gage.groupby(gage.index.month).mean(), 'k', label='Observed')
    ax.plot(flow.groupby(gage.index.month).mean(), 'b', label='Competitor')
    if gage_id != 'USGS08353000':
        ax.plot(RGB_Futures.groupby(gage.index.month).mean().sum(axis=1), 'r', label='TOVA')

    ax.set_xticks(range(1,13))
    ax.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
    ax.legend(fontsize=12)
    ax.set_ylabel('Water availability\n(Acre Feet/day)', fontsize=12)
    ax.set_title(f'{gage_region} -- Gage id:{gage_id}', fontsize=12)

    plt.tight_layout()
    fig.savefig(f'../figures/showcasing/{gage_region}_{gage_id}_PCR-GLOBWB2.0_IAM.png')
    #plt.show()
    plt.close()

    for year in range(2001, 2016):

        fig, ax = plt.subplots(1,1, figsize=(6,4))
        gage_y = gage.truncate(before=f'{year}-01-31', after=f'{year}-12-31')
        flow_y = flow.truncate(before=f'{year}-01-31', after=f'{year}-12-31')
        if gage_id != 'USGS08353000':
            RGB_Futures_y = RGB_Futures.truncate(before=f'{year}-01-31', after=f'{year}-12-31')
        
        ax.plot(gage_y, 'k', label='Observed')
        ax.plot(flow_y, 'b', label='Competitor')
        if gage_id != 'USGS08353000':
            ax.plot(RGB_Futures_y.sum(axis=1), 'r', label='TOVA')
        ax.legend(fontsize=12)
        ax.set_ylabel('Water availability\n(Acre Feet/day)', fontsize=12)
        ax.set_title(f'{gage_region} -- Gage id:{gage_id} -- Year {year}', fontsize=12)
        
        plt.tight_layout()
        fig.savefig(f'../figures/showcasing/{gage_region}_{gage_id}_PCR-GLOBWB2.0_IAM_{year}.png')
        #plt.show()
        plt.close()


    if gage_id != 'USGS08353000':
        
        fig,ax = plt.subplots(1, 1, figsize=(4,4))
        ax.plot(gage, flow, 'ko')
        ax.set_xlabel('Observed (AFD)', fontsize=12)
        ax.set_ylabel('PCR-GLOBWB 2.0 (AFD)', fontsize=12)
        ax.axline((0, 0), slope=1, color="r", linestyle='--')
        ax.grid()
        plt.tight_layout()
        fig.savefig(f'../figures/showcasing/{gage_region}_{gage_id}_PCR-GLOBWB2.0_scatter.png')
        #plt.show()
        plt.close()

        fig,ax = plt.subplots(1, 1, figsize=(4,4))
        ax.plot(gage, RGB_Futures.sum(axis=1), 'ko')
        ax.set_xlabel('Observed (AFD)', fontsize=12)
        ax.set_ylabel('RGB Futures (AFD)', fontsize=12)
        ax.axline((0, 0), slope=1, color="r", linestyle='--')
        ax.grid()
        plt.tight_layout()
        fig.savefig(f'../figures/showcasing/{gage_region}_{gage_id}_RGB_Futures_scatter.png')
        #plt.show()
        plt.close()


    # Statistics
    statistics.loc['USGS'] = [
        gage.min()[0], gage.quantile(0.05)[0], gage.quantile(0.25)[0], gage.quantile(0.5)[0],
        gage.quantile(0.75)[0], gage.quantile(0.95)[0], gage.max()[0], gage.mean()[0], gage.std()[0],
        -999.99
    ]
    
    if gage_id != 'USGS08353000':
        

        statistics.loc['RGB_Futures'] = [
            RGB_Futures.sum(axis=1).min(), RGB_Futures.sum(axis=1).quantile(0.05), 
            RGB_Futures.sum(axis=1).quantile(0.25), RGB_Futures.sum(axis=1).quantile(0.5),
            RGB_Futures.sum(axis=1).quantile(0.75), RGB_Futures.sum(axis=1).quantile(0.95),
            RGB_Futures.sum(axis=1).max(), RGB_Futures.sum(axis=1).mean(), 
            RGB_Futures.sum(axis=1).std(), nash(RGB_Futures.sum(axis=1), gage.values.flatten())
        ]

    statistics.loc['PCR_GLOB'] = [
        flow.min()[0], flow.quantile(0.05)[0], flow.quantile(0.25)[0], flow.quantile(0.5)[0],
        flow.quantile(0.75)[0], flow.quantile(0.95)[0], flow.max()[0], flow.mean()[0], flow.std()[0],
        nash(flow.values.flatten(), gage.values.flatten())
    ]
    
    
    statistics.to_csv('../data/RGB_Futures/statistics_{}_closest.csv'.format(gage_id))

#dates = pd.date_range(start='1958-01-31', end='2015-12-31', freq='ME')
#discharge.index = dates
#discharge.index.name = 'Dates'
#discharge.to_csv('../data/PCR_GLOBWB2/discharge_monthAvg_output_1958-01-31_to_2015-12-31.csv')



