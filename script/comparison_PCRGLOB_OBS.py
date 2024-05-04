

import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Read PCR GLOBWB 2.0 simulations
PCR = pd.read_csv(
    '../data/PCR_GLOBWB2/discharge_monthAvg_output_1958-01-31_to_2015-12-31.csv',
    index_col='Dates', parse_dates=True)
PCR = PCR.truncate(before='1981-01-31', after='2015-12-31')
# Annual total for the water year (starting in October)
PCR_annual = PCR.groupby(pd.Grouper(freq='YE-SEP')).sum()
# Remove first and last year beceause they are incomplete water year 
PCR_annual.drop(PCR_annual.index[-1], inplace=True)
PCR_annual.drop(PCR_annual.index[0], inplace=True)

# Large basins
PCR_LB = pd.read_csv(
    '../data/PCR_GLOBWB2/discharge_monthAvg_USGSgages_tillTX_1958-01-31_to_2015-12-31.csv',
    index_col='Dates', parse_dates=True)


def nash(predictions, targets):
    return 1-(np.sum((targets-predictions)**2)/np.sum((targets-np.mean(targets))**2))

gage_id = []   
bias = []
nse = []
rmse = []
area = []
corr = []

gage_id_LB = []
bias_LB = []
nse_LB = []
rmse_LB = []
area_LB = []
corr_LB = []

# Loop over the gages
for gage in PCR.columns:

    # The gage is located in United States
    if 'camels_' in gage:

        folder = 'camels'
        attributes = pd.read_csv(
            '../data/caravan/attributes/camels/attributes_other_camels.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])
        
    # The gage is located in Australia
    elif 'camelsaus_' in gage:
        
        folder = 'camelsaus'
        attributes = pd.read_csv(
            '../data/caravan/attributes/camelsaus/attributes_other_camelsaus.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])
        
    # The gage is located in Brazil
    elif 'camelsbr_' in gage:

        folder = 'camelsbr'
        attributes = pd.read_csv(
            '../data/caravan/attributes/camelsbr/attributes_other_camelsbr.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])
    
    # The gage is located in Chile
    elif 'camelscl_' in gage:

        folder = 'camelscl'
        attributes = pd.read_csv(
            '../data/caravan/attributes/camelscl/attributes_other_camelscl.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])
    
    # The gage is located in the United Kingdom
    elif 'camelsgb_' in gage:

        folder = 'camelsgb'
        attributes = pd.read_csv(
            '../data/caravan/attributes/camelsgb/attributes_other_camelsgb.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])

    # The gage is located in the United States or Canada
    elif 'hysets_' in gage:

        folder = 'hysets'
        attributes = pd.read_csv(
            '../data/caravan/attributes/hysets/attributes_other_hysets.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])

    # The gage is located in Central Europe
    elif 'lamah' in gage:

        folder = 'lamah'
        attributes = pd.read_csv(
            '../data/caravan/attributes/lamah/attributes_other_lamah.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])
        
    # The gage is located in the Alps
    elif 'camelsch' in gage:
        
        folder = 'camelsch'
        attributes = pd.read_csv(
            '../data/caravan/attributes/camelsch/attributes_other_camelsch.csv',
            index_col='gauge_id', usecols=['gauge_id', 'area', 'gauge_lat', 'gauge_lon'])


    # Read Observed Data
    obs = pd.read_csv('G:/Caravan/Caravan/timeseries/csv/{}/{}.csv'.format(folder, gage), 
                      index_col='date', parse_dates=True, usecols=['date', 'streamflow'])
    obs = obs.truncate(before='1981-01-31', after='2015-12-31')
    obs_monthAvg = obs.resample('ME').mean()
    
    # Annual total for the water year (starting in October)
    obs_annual = obs.groupby(pd.Grouper(freq='YE-SEP')).sum()
    # Remove first and last year beceause they are incomplete water year 
    obs_annual.drop(obs_annual.index[-1], inplace=True)
    obs_annual.drop(obs_annual.index[0], inplace=True)

    # Only select gage if data are found
    if True in np.isnan(PCR[gage].values):
        continue

    # ID of the gage
    gage_id.append(gage)

    # Bias, NSE, and RMSE
    bias.append(
        (np.mean(PCR[gage] - obs_monthAvg['streamflow'])) / 
        np.mean(obs_monthAvg['streamflow']) * 100)

    nse.append(nash(PCR[gage], obs_monthAvg['streamflow']))

    rmse.append(
        np.sqrt(np.mean((PCR[gage] - obs_monthAvg['streamflow'])**2)))
    
    # Correlation between annual mean
    corr.append(np.corrcoef(PCR_annual[gage], obs_annual['streamflow'])[0,1])


# Loop over the gages
for gage in PCR_LB.columns:

    # Only select gage if data are found
    if True in np.isnan(PCR_LB[gage].values):
        continue

    
    # Read Observed Data
    obs_LB = pd.read_csv('../data/USGS/USGS_{}.csv'.format(gage), 
                      index_col='Dates', parse_dates=True, usecols=['Dates', 'Flow_cfs'])
    obs_LB = obs_LB.truncate(before='1958-01-31', after='2015-12-31')
    
    if obs_LB['Flow_cfs'].shape[0]==0:
        continue
    
    obs_LB_monthAvg = obs_LB.resample('ME').mean()
    
    PCR_LB_plot = PCR_LB.truncate(before='{}-{}-{}'.format(
        obs_LB_monthAvg.index[0].year, obs_LB_monthAvg.index[0].month, obs_LB_monthAvg.index[0].day),
                             after='{}-{}-{}'.format(
        obs_LB_monthAvg.index[-1].year, obs_LB_monthAvg.index[-1].month, obs_LB_monthAvg.index[-1].day))

    # Annual total for the water year (starting in October)
    PCR_LB_annual = PCR_LB_plot.groupby(pd.Grouper(freq='YE-SEP')).sum()
    # Remove first and last year beceause they are incomplete water year 
    PCR_LB_annual.drop(PCR_LB_annual.index[-1], inplace=True)
    PCR_LB_annual.drop(PCR_LB_annual.index[0], inplace=True)

    # Annual total for the water year (starting in October)
    obs_LB_annual = obs_LB.groupby(pd.Grouper(freq='YE-SEP')).sum()
    # Remove first and last year beceause they are incomplete water year 
    obs_LB_annual.drop(obs_LB_annual.index[-1], inplace=True)
    obs_LB_annual.drop(obs_LB_annual.index[0], inplace=True)

    # ID of the gage
    gage_id_LB.append(gage)

    # Bias, NSE, and RMSE
    bias_LB.append(
        (np.mean(PCR_LB_plot[gage] - obs_LB_monthAvg['Flow_cfs'])) / 
        np.mean(obs_LB_monthAvg['Flow_cfs']) * 100)

    nse_LB.append(nash(PCR_LB_plot[gage], obs_LB_monthAvg['Flow_cfs']))

    rmse_LB.append(
        np.sqrt(np.mean((PCR_LB_plot[gage] - obs_LB_monthAvg['Flow_cfs'])**2)))
    
    # Correlation between annual mean
    corr_LB.append(np.corrcoef(PCR_LB_annual[gage], obs_LB_annual['Flow_cfs'])[0,1])
    
   



fig, ax = plt.subplots(1,1, figsize=(6,6))
ax.plot(np.sort(nse), 'k', label='Caravan (<2,000 km$^2$)')
ax.plot(np.sort(nse_LB), 'b', label='Basins (>2,000 km$^2$)')
ax.set_ylim(-1,1)
ax.grid()
ax.set_xlabel('Gages',fontsize=12)
ax.set_ylabel('NSE',fontsize=12)
ax.set_title('PCR-GLOBWB2.0',fontsize=12)
plt.tight_layout()
fig.savefig('../figures/NSE_PCRGLOBWB2.0.png')
plt.show()


fig, ax = plt.subplots(1,1, figsize=(6,6))
ax.plot(np.sort(corr), 'k', label='Caravan (<2,000 km$^2$)')
ax.plot(np.sort(corr_LB), 'b', label='Basins (>2,000 km$^2$)')
ax.set_ylim(-0.5,1)
ax.grid()
ax.set_xlabel('Gages',fontsize=12)
ax.set_ylabel('Correlation',fontsize=12)
ax.set_title('PCR-GLOBWB2.0',fontsize=12)
plt.tight_layout()
fig.savefig('../figures/Correlation_PCRGLOBWB2.0_w_LB.png')
plt.show()

