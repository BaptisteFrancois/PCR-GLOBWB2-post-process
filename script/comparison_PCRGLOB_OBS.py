

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


def nash(predictions, targets):
    return 1-(np.sum((targets-predictions)**2)/np.sum((targets-np.mean(targets))**2))

gage_id = []   
bias = []
nse = []
rmse = []
area = []
corr = []

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

    # Area of the basin
    area.append(attributes.loc[gage]['area'])

    # Bias, NSE, and RMSE
    bias.append(
        (np.mean(PCR[gage] - obs_monthAvg['streamflow'])) / 
        np.mean(obs_monthAvg['streamflow']) * 100)

    nse.append(nash(PCR[gage], obs_monthAvg['streamflow']))

    rmse.append(
        np.sqrt(np.mean((PCR[gage] - obs_monthAvg['streamflow'])**2)))
    
    # Correlation between annual mean
    corr.append(np.corrcoef(PCR_annual[gage], obs_annual['streamflow'])[0,1])
    
    

fig, ax = plt.subplots(1,1, figsize=(6,6))
ax.plot(np.sort(nse), 'k')
ax.set_ylim(-1,1)
ax.grid()
ax.set_xlabel('Gages',fontsize=12)
ax.set_ylabel('NSE',fontsize=12)
ax.set_title('PCR-GLOBWB2.0',fontsize=12)
plt.tight_layout()
fig.savefig('../figures/NSE_PCRGLOBWB2.0.png')
plt.show()


fig, ax = plt.subplots(1,1, figsize=(6,6))
ax.plot(np.sort(corr), 'k')
ax.set_ylim(-0.5,1)
ax.grid()
ax.set_xlabel('Gages',fontsize=12)
ax.set_ylabel('Correlation',fontsize=12)
ax.set_title('PCR-GLOBWB2.0',fontsize=12)
plt.tight_layout()
fig.savefig('../figures/Correlation_PCRGLOBWB2.0.png')
plt.show()

