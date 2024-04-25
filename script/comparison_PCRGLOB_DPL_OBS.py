
# ENV: 

import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Read PCR GLOBWB 2.0 simulations
PCR = pd.read_csv(
    '../data/PCR_GLOBWB2/discharge_monthAvg_output_1958-01-31_to_2015-12-31.csv',
    index_col='Dates', parse_dates=True)
PCR = PCR.truncate(before='2008-10-01', after='2014-09-30')
# Annual total for the water year (starting in October)
PCR_annual = PCR.groupby(pd.Grouper(freq='YE-SEP')).sum()




# Read DPL simulation
DPL_label_file = open('../data/DPL-caravan_Maharjan_et_al_2024/Global_tr_val_split_data.pkl', 'rb')
DPL_label = pickle.load(DPL_label_file)

## Validation simulation period is 2008-10-01 to 2014-09-30
DPL_flow = pd.read_pickle('../data/DPL-caravan_Maharjan_et_al_2024/G1_dPL_val_TS.pkl')


def nash(predictions, targets):
    return 1-(np.sum((targets-predictions)**2)/np.sum((targets-np.mean(targets))**2))

gage_id = []
num_gages = 0
area = []

bias_PCR_OBS = []
biasV_PCR_OBS = []
nse_PCR_OBS = []
rmse_PCR_OBS = []
corr_PCR_OBS = []

bias_DPL_OBS = []
biasV_DPL_OBS = []
nse_DPL_OBS = []
rmse_DPL_OBS = []
corr_DPL_OBS = []

rmse_by_month_PCR_OBS = np.zeros((12,283))
corr_by_month_PCR_OBS = np.zeros_like(rmse_by_month_PCR_OBS)
rmse_by_month_DPL_OBS = np.zeros_like(rmse_by_month_PCR_OBS)
corr_by_month_DPL_OBS = np.zeros_like(rmse_by_month_PCR_OBS)

DPL = pd.DataFrame()

count_basins_for_comparison = 0
# Loop over the gages
for gage in PCR.columns:

    # Look for the experiement in the DPL simulation that include the current gage
    dpl_sim_all_exp = pd.DataFrame()
    obs_dpl = pd.DataFrame()
    for exp in DPL_label['val'].keys():
        if gage in DPL_label['val'][exp]:
            index = [i for i, label in enumerate(DPL_label['val'][exp]) if label == gage][0]
            dpl_sim_all_exp[exp] = DPL_flow[int(exp)][index]['sim']


    # Only select gage if data are found
    if dpl_sim_all_exp.empty:
        continue

    # Average the DPL simulation across all available training sets
    DPL[gage] = dpl_sim_all_exp.mean(axis=1)
    DPL.index = pd.date_range(start='2008-10-01', end='2014-09-30')
    DPL_monthAvg = DPL.resample('ME').mean()
    DPL_annual = DPL.groupby(pd.Grouper(freq='YE-SEP')).sum()

    
    # We are counting the number of basins for which we compare PCR and DPL-HBV
    count_basins_for_comparison += 1

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
    obs = obs.truncate(before='2008-10-01', after='2014-09-30')
    obs_monthAvg = obs.resample('ME').mean()
    
    # Annual total for the water year (starting in October)
    obs_annual = obs.groupby(pd.Grouper(freq='YE-SEP')).sum()
    

    

    # Only select gage if data are found
    if True in np.isnan(PCR[gage].values):
        continue

    # ID of the gage
    gage_id.append(gage)

    # Area of the basin
    area.append(attributes.loc[gage]['area'])

    # PCR Performance
    # ---------------
    # Bias, NSE, RMSE and correlation between PCR and Observed
    bias_PCR_OBS.append(
        (np.mean(PCR[gage] - obs_monthAvg['streamflow'])) / 
        np.mean(obs_monthAvg['streamflow']) * 100)
    
    biasV_PCR_OBS.append(
        (np.std(PCR[gage])**2 - np.std(obs_monthAvg['streamflow'])**2) /
        np.std(obs_monthAvg['streamflow'])**2 * 100)

    nse_PCR_OBS.append(nash(PCR[gage], obs_monthAvg['streamflow']))

    rmse_PCR_OBS.append(
        np.sqrt(np.mean((PCR[gage] - obs_monthAvg['streamflow'])**2)))
    
    # Correlation between annual mean
    corr_PCR_OBS.append(np.corrcoef(PCR_annual[gage], obs_annual['streamflow'])[0,1])
    
    # DPL Performance
    # ---------------
    # Bias, NSE, RMSE and correlation between DPL and Observed
    bias_DPL_OBS.append(
        (np.mean(DPL_monthAvg[gage] - obs_monthAvg['streamflow'])) / 
        np.mean(obs_monthAvg['streamflow']) * 100)
    
    biasV_DPL_OBS.append(
        (np.std(DPL[gage])**2 - np.std(obs_monthAvg['streamflow'])**2) /
        np.std(obs_monthAvg['streamflow'])**2 * 100)

    nse_DPL_OBS.append(nash(DPL_monthAvg[gage], obs_monthAvg['streamflow']))

    rmse_DPL_OBS.append(
        np.sqrt(np.mean((DPL_monthAvg[gage] - obs_monthAvg['streamflow'])**2)))
    
    # Correlation between annual mean
    corr_DPL_OBS.append(np.corrcoef(DPL_annual[gage], obs_annual['streamflow'])[0,1])

    for im, month in enumerate(range(1,13)):

        PCR_by_month = PCR[PCR.index.month == month]
        DPL_by_month = DPL_monthAvg[DPL_monthAvg.index.month == month]
        obs_by_month = obs_monthAvg[obs_monthAvg.index.month == month]
        
        # PCR Performance
        rmse_by_month_PCR_OBS[im,num_gages] = (
            np.sqrt(np.mean((PCR_by_month[gage] - obs_by_month['streamflow'])**2)))
        
        corr_by_month_PCR_OBS[im,num_gages] = \
            np.corrcoef(PCR_by_month[gage], obs_by_month['streamflow'])[0,1]
        
        rmse_by_month_DPL_OBS[im,num_gages] = (
            np.sqrt(np.mean((DPL_by_month[gage] - obs_by_month['streamflow'])**2)))
        
        corr_by_month_DPL_OBS[im,num_gages] = \
            np.corrcoef(DPL_by_month[gage], obs_by_month['streamflow'])[0,1]
    
    num_gages += 1
    print(num_gages)

fig, ax = plt.subplots(1,1, figsize=(4,4))
ax.plot(np.sort(nse_PCR_OBS), 'b')
ax.plot(np.sort(nse_DPL_OBS), 'r')
ax.set_ylim(-1,1)
ax.grid()
ax.set_xlabel('Gages', fontsize=12)
ax.set_ylabel('NSE', fontsize=12)
ax.legend(['PCR-GLOBWB 2.0', 'DPL-HBV (Global Training)'])
plt.tight_layout()
fig.savefig('../figures/CDF_Nash_PCR_DPL_OBS.png')
#plt.show()
plt.close()


fig, ax = plt.subplots(1,1, figsize=(4,4))
ax.plot(np.sort(corr_PCR_OBS), 'b')
ax.plot(np.sort(corr_DPL_OBS), 'r')
ax.set_ylim(-0.5,1)
ax.grid()
ax.set_xlabel('Gages', fontsize=12)
ax.set_ylabel('Correlation', fontsize=12)
ax.legend(['PCR-GLOBWB 2.0', 'DPL-HBV (Global Training)'])
plt.tight_layout()
fig.savefig('../figures/CDF_CorrAnnual_PCR_DPL_OBS.png')
#plt.show()
plt.close()

fig, ax = plt.subplots(1,1, figsize=(4,4))
ax.plot(np.sort(bias_PCR_OBS), 'b')
ax.plot(np.sort(bias_DPL_OBS), 'r')
#ax.set_yscale('log')
ax.grid()
ax.set_xlabel('Gages', fontsize=12)
ax.set_ylabel('Bias (%)', fontsize=12)
ax.legend(['PCR-GLOBWB 2.0', 'DPL-HBV (Global Training)'])
plt.tight_layout()
fig.savefig('../figures/CDF_Bias_PCR_DPL_OBS.png')
#plt.show()
plt.close()

fig, ax = plt.subplots(1,1, figsize=(4,4))
ax.plot(np.sort(biasV_PCR_OBS), 'b')
ax.plot(np.sort(biasV_DPL_OBS), 'r')
#ax.set_yscale('log')
ax.grid()
ax.set_xlabel('Gages', fontsize=12)
ax.set_ylabel('Bias Variance (%)', fontsize=12)
ax.legend(['PCR-GLOBWB 2.0', 'DPL-HBV (Global Training)'])
plt.tight_layout()
fig.savefig('../figures/CDF_BiasV_PCR_DPL_OBS.png')
#plt.show()
plt.close()

fig, ax = plt.subplots(1,1, figsize=(4,4))
ax.plot(np.sort(rmse_PCR_OBS), 'b')
ax.plot(np.sort(rmse_DPL_OBS), 'r')
ax.set_yscale('log')
ax.grid()
ax.set_xlabel('Gages', fontsize=12)
ax.set_ylabel('RMSE (mm.day-1)', fontsize=12)
ax.legend(['PCR-GLOBWB 2.0', 'DPL-HBV (Global Training)'])
plt.tight_layout()
fig.savefig('../figures/CDF_RMSE_PCR_DPL_OBS.png')
#plt.show()
plt.close()

fig, ax = plt.subplots(1,1, figsize=(4,4))
ax.plot(nse_PCR_OBS, nse_DPL_OBS, 'xk')
ax.plot([-1,1], [-1,1], '--r')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.grid()
ax.set_xlabel('PCR-GLOBWB 2.0', fontsize=12)
ax.set_ylabel('DPL-HBV (Global Training)', fontsize=12)
ax.set_title('Nash Sutcliffe Efficiency', fontsize=12)
plt.tight_layout()
fig.savefig('../figures/Scatterplot_Nash_PCR_DPL_OBS.png')
#plt.show()
plt.close()

fig, ax = plt.subplots(1,1, figsize=(4,4))
ax.plot(corr_PCR_OBS, corr_DPL_OBS, 'xk')
ax.plot([-1,1], [-1,1], '--r')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.grid()
ax.set_xlabel('PCR-GLOBWB 2.0', fontsize=12)
ax.set_ylabel('DPL-HBV (Global Training)', fontsize=12)
ax.set_title('Correlation annual flow', fontsize=12)
plt.tight_layout()
fig.savefig('../figures/Scatterplot_CorrAnnual_PCR_DPL_OBS.png')
#plt.show()
plt.close()

fig, ax = plt.subplots(1,1, figsize=(4,4))
ax.plot(bias_PCR_OBS, bias_DPL_OBS, 'xk')
ax.axline((1, 1), slope=1, color="r", linestyle='--')
ax.grid()
ax.set_xlabel('PCR-GLOBWB 2.0', fontsize=12)
ax.set_ylabel('DPL-HBV (Global Training)', fontsize=12)
ax.set_title('Bias (%)', fontsize=12)
plt.tight_layout()
fig.savefig('../figures/Scatterplot_Bias_PCR_DPL_OBS.png')
#plt.show()
plt.close()

fig, ax = plt.subplots(1,1, figsize=(4,4))
ax.plot(biasV_PCR_OBS, biasV_DPL_OBS, 'xk')
#ax.plot([-1,1], [-1,1], '--r')
#ax.set_xlim(-1,1)
#ax.set_ylim(-1,1)
ax.grid()
ax.set_xlabel('PCR-GLOBWB 2.0', fontsize=12)
ax.set_ylabel('DPL-HBV (Global Training)', fontsize=12)
ax.set_title('Bias (%)', fontsize=12)
#ax.set_xscale('log')
#ax.set_yscale('log')
plt.tight_layout()
fig.savefig('../figures/Scatterplot_BiasV_PCR_DPL_OBS.png')
#plt.show()
plt.close()

fig, ax = plt.subplots(1,1, figsize=(4,4))
ax.plot(rmse_PCR_OBS, rmse_DPL_OBS, 'xk')
ax.axline((1, 1), slope=1, color="r", linestyle='--')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()
ax.set_xlabel('PCR-GLOBWB 2.0', fontsize=12)
ax.set_ylabel('DPL-HBV (Global Training)', fontsize=12)
ax.set_title('RMSE (mm.day-1)', fontsize=12)
plt.tight_layout()
fig.savefig('../figures/Scatterplot_RMSE_PCR_DPL_OBS.png')
#plt.show()
plt.close()


# 4x3 subplots
month_label=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
custom_lines = [Line2D([0], [0], color='red', lw=1),
                Line2D([0], [0], color='blue', linestyle='--', lw=1),
                ]

fig, axes = plt.subplots(4,3, figsize=(12,8))
for i, ax in enumerate(axes.flatten()):
    ax.plot(np.sort(rmse_by_month_PCR_OBS[i,:]), 'b')
    ax.plot(np.sort(rmse_by_month_DPL_OBS[i,:]), 'r')
    
    ax.grid()
    ax.set_ylabel('RMSE (mm.day-1)', fontsize=12)
    ax.set_xlabel('Gages', fontsize=10)
    ax.set_title(month_label[i], fontsize=12)
plt.tight_layout()
plt.figlegend(custom_lines, ['DPL-HBV', 'PCR-GLOBWB 2.0'],
    loc = 'lower center', ncol=2, bbox_to_anchor=(0.5,-0.004), fontsize=12)
plt.subplots_adjust(bottom=0.1)
fig.savefig('../figures/CDF_RMSE_month_PCR_DPL_OBS.png')
#plt.show()
plt.close()

fig, axes = plt.subplots(4,3, figsize=(12,8))
for i, ax in enumerate(axes.flatten()):
    ax.plot(np.sort(rmse_by_month_PCR_OBS[i,:]), 'b')
    ax.plot(np.sort(rmse_by_month_DPL_OBS[i,:]), 'r')
    ax.set_ylim(0,1.5)
    ax.grid()
    ax.set_ylabel('RMSE (mm.day-1)', fontsize=12)
    ax.set_xlabel('Gages', fontsize=10)
    ax.set_title(month_label[i], fontsize=12)
plt.tight_layout()
plt.figlegend(custom_lines, ['DPL-HBV', 'PCR-GLOBWB 2.0'],
    loc = 'lower center', ncol=2, bbox_to_anchor=(0.5,-0.004), fontsize=12)
plt.subplots_adjust(bottom=0.1)
fig.savefig('../figures/CDF_RMSE_month_PCR_DPL_OBS_zoom.png')
#plt.show()
plt.close()


fig, axes = plt.subplots(4,3, figsize=(12,8))
for i, ax in enumerate(axes.flatten()):
    ax.plot(np.sort(corr_by_month_PCR_OBS[i,:]), 'b')
    ax.plot(np.sort(corr_by_month_DPL_OBS[i,:]), 'r')
    ax.set_ylim(-1,1)
    ax.grid()
    ax.set_ylabel('Correlation', fontsize=12)
    ax.set_xlabel('Gages', fontsize=10)
    ax.set_title(month_label[i], fontsize=12)
plt.tight_layout()
plt.figlegend(custom_lines, ['DPL-HBV', 'PCR-GLOBWB 2.0'],
    loc = 'lower center', ncol=2, bbox_to_anchor=(0.5,-0.004), fontsize=12)
plt.subplots_adjust(bottom=0.1)
fig.savefig('../figures/CDF_Corr_month_PCR_DPL_OBS.png')
#plt.show()
plt.close()




fig, axes = plt.subplots(4,3, figsize=(12,8))
for i, ax in enumerate(axes.flatten()):
    ax.plot(rmse_by_month_PCR_OBS[i,:], rmse_by_month_DPL_OBS[i,:], 'xk')
    
    ax.grid()
    ax.set_ylabel('PCR-GLOBWB 2.0', fontsize=12)
    ax.set_xlabel('DPL-HBV', fontsize=10)
    ax.set_title(month_label[i], fontsize=12)
plt.suptitle('RMSE (mm.day-1)', fontsize=12)
plt.tight_layout()
fig.savefig('../figures/Scatterplot_RMSE_month_PCR_DPL_OBS.png')
#plt.show()
plt.close()


# Time series plot
id_plotted = ['camels_11522500']
for id_plotted in gage_id:
    iid = id_plotted.find('_')
    obs = pd.read_csv('G:/Caravan/Caravan/timeseries/csv/{}/{}.csv'.format(id_plotted[:iid], id_plotted), 
                        index_col='date', parse_dates=True, usecols=['date', 'streamflow'])
    obs = obs.truncate(before='2008-10-01', after='2014-09-30')
    obs_monthAvg = obs.resample('ME').mean()

    fig, ax = plt.subplots(1,1, figsize=(12,6))
    ax.plot(obs_monthAvg.index, obs_monthAvg['streamflow'], 'k', label='Observed')
    ax.plot(PCR.index, PCR[id_plotted], 'b', label='PCR-GLOBWB 2.0')
    ax.plot(DPL_monthAvg.index, DPL_monthAvg[id_plotted], 'r', label='DPL-HBV (Global Training)')

    ax.grid()
    ax.set_ylabel('Streamflow (mm.day-1)', fontsize=12)
    ax.set_title(id_plotted, fontsize=12)
    ax.legend()
    plt.tight_layout()
    fig.savefig('../figures/timeseries/TimeSeries_{}.png'.format(id_plotted))
    #plt.show()
    plt.close()
