
import numpy as np
import pandas as pd

import itertools 
import matplotlib.pyplot as plt

m3s_to_afd = 70.0452
cfs_to_afd = 1.98347

def rrv(supply, targets):
    reliability = 1 - np.sum(supply <= targets) / len(targets)
    
    failure = np.where(supply <= targets, 1, 0)
    bool_failure = supply[supply <= targets] == 1
    duration_failures = \
        [sum(g) for bool_failure, g in itertools.groupby(
         failure) if bool_failure]

    if len(duration_failures)==0:
        duration_failures = [0]        

    gap = np.where(supply <= targets, targets - supply, 0)
    supply_gap = \
        [sum(g) for bool_failure, g in itertools.groupby(
         gap) if bool_failure]
    
    if len(supply_gap)==0:  
        supply_gap = [0]
    
    #supply_gap = np.sum(np.where(supply <= targets, targets - supply, 0))
    #avg_supply_gap = supply_gap / len(duration_failures )

    return reliability, duration_failures, supply_gap

gages = (
    ('RGB', 'USGS08362500', 32.88, -107.29, '08362500_RIO_GRANDE_BLW_CABALLO_DAM.csv'),
    ('RGB', 'USGS08361000', 33.15, -107.21, '08361000_RIO_GRANDE_BELOW_ELEPHANT_BUTTE_DAM_fill.csv'),
    ('RGB', 'USGS08317400', 35.62, -106.32, 'USGS_08317400.csv')
)

data_center_demand = (0.05, 0.1, 0.15, 0.2, 0.25, 0.3)

discharge = pd.DataFrame()


with pd.ExcelWriter('../data/RGB_Futures/Performance_metrics_datacenter_RioGrande.xlsx') as writer:

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

            PCR = pd.read_csv('../data/RGB_Futures/PCR_GLOBWB2_{}.csv'.format(gage_id), index_col='Date',
                            parse_dates=True)
            PCR = PCR.resample('ME').mean()
            PCR = PCR.truncate(before=start, after=end)

            

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

            PCR = pd.read_csv('../data/RGB_Futures/PCR_GLOBWB2_{}.csv'.format(gage_id), index_col='Date',
                            parse_dates=True)
            PCR = PCR.resample('ME').mean()
            PCR = PCR.truncate(before=start, after=end)

        elif gage_id == 'USGS08317400':

            # Common period for the observed and simulated data (PCR and RGB Futures)
            start = '2001-01-31'
            end = '2015-12-31'

            gage = pd.read_csv('../data/RGB_Futures/{}'.format(name_file),
                index_col='Dates', parse_dates=True, usecols=['Dates', 'Flow_cfs']) * cfs_to_afd
            gage = gage.resample('ME').mean()
            gage = gage.truncate(before=start, after=end)
            gage.columns = ['gauges']

            RGB_Futures = pd.read_csv(
                '../data/RGB_Futures/rio_grande_output_file_02152023_baseline_cost.csv',
                index_col='date', parse_dates=True, usecols=['date', 'cochiti_release_AF', 'cochiti_spill_AF'])
            RGB_Futures = RGB_Futures.resample('ME').mean()
            RGB_Futures = RGB_Futures.truncate(before=start, after=end)

            PCR = pd.read_csv('../data/RGB_Futures/PCR_GLOBWB2_{}.csv'.format(gage_id), index_col='Date',
                            parse_dates=True)
            PCR = PCR.resample('ME').mean()
            PCR = PCR.truncate(before=start, after=end)

        performance = pd.DataFrame(columns=['demand','reliability', 'mean duration_failure', 
                                            'max duration_failure', 'avg_supply_gap(by episode)'])

        for k, demand_perct in enumerate(data_center_demand, start=1):
            

            demand = np.ones(gage['gauges'].shape) * demand_perct * gage['gauges'].mean()

            fig, ax = plt.subplots(3, 1, figsize=(8,12))

            ax[0].plot(gage.index, gage['gauges'], 'k', label='Observed')
            ax[0].plot(gage.index, demand, 'r', label='Demand')
            ax[0].fill_between(gage.index, gage['gauges'].values, demand, gage['gauges'].values <=demand,
                            color='red', alpha=0.5)

            ax[1].plot(RGB_Futures.index, RGB_Futures.sum(axis=1), 'k', label='RGB Futures')
            ax[1].plot(RGB_Futures.index, demand, 'r', label='Demand')
            ax[1].fill_between(RGB_Futures.index, RGB_Futures.sum(axis=1).values, demand, RGB_Futures.sum(axis=1).values <=demand,
                                color='red', alpha=0.5)
            
            ax[2].plot(PCR.index, PCR['PCR_GLOBWB2'], 'k', label='PCR GLOBWB2')
            ax[2].plot(PCR.index, demand, 'r', label='Demand')
            ax[2].fill_between(PCR.index, PCR['PCR_GLOBWB2'].values, demand, PCR['PCR_GLOBWB2'].values <=demand,
                                color='red', alpha=0.5)
            
            ax[0].set_title('Observed')
            ax[1].set_title('RGB Futures')
            ax[2].set_title('PCR GLOBWB2')

            for axx in ax:
                axx.legend()
                axx.grid()
                axx.set_ylabel('Flow [AF/day]')

            plt.suptitle('{} -- Demand = {}x avg Obs flow'.format(gage_id, demand_perct))

            plt.tight_layout()
            fig.savefig('../figures/sketch_datacenter_RioGrande/{}_{}.png'.format(gage_id, demand_perct))
            #plt.show()
            plt.close()

            # Calculate the supply
            rel_PCR, dur_PCR, supply_gap_PCR = rrv(PCR['PCR_GLOBWB2'].values, demand)
            rel_RBG, dur_RBG, supply_gap_RBG = rrv(RGB_Futures.sum(axis=1).values, demand)
            rel_OBS, dur_OBS, supply_gap_OBS = rrv(gage.values.ravel(), demand)

            performance.loc['OBS_{}'.format(k)] = \
                [demand_perct, rel_OBS, np.mean(dur_OBS), np.max(dur_OBS), np.max(supply_gap_OBS)]
            performance.loc['RGB_Futures_{}'.format(k)] = \
                [demand_perct, rel_RBG, np.mean(dur_RBG), np.max(dur_RBG), np.max(supply_gap_RBG)]
            performance.loc['PCR_{}'.format(k)] = \
                [demand_perct, rel_PCR, np.mean(dur_PCR), np.max(dur_PCR), np.max(supply_gap_PCR)]
            
            OBS_performance = performance.loc[performance.index.str.contains('OBS')]
            RGB_Futures_performance = performance.loc[performance.index.str.contains('RGB_Futures')]
            PCR_performance = performance.loc[performance.index.str.contains('PCR')]

            fig, ax = plt.subplots(2, 2, figsize=(8,6))

            ax[0,0].plot(OBS_performance['demand'], OBS_performance['reliability'], 'ko-', label='Observed')
            ax[0,0].plot(RGB_Futures_performance['demand'], RGB_Futures_performance['reliability'], 'r-', label='RGB Futures')
            ax[0,0].plot(PCR_performance['demand'], PCR_performance['reliability'], 'b--', label='PCR GLOBWB2')

            ax[0,1].plot(OBS_performance['demand'], OBS_performance['mean duration_failure'], 'ko-', label='Observed')
            ax[0,1].plot(RGB_Futures_performance['demand'], RGB_Futures_performance['mean duration_failure'], 'ro-', label='RGB Futures')
            ax[0,1].plot(PCR_performance['demand'], PCR_performance['mean duration_failure'], 'b--', label='PCR GLOBWB2')

            ax[1,0].plot(OBS_performance['demand'], OBS_performance['max duration_failure'], 'ko-', label='Observed')
            ax[1,0].plot(RGB_Futures_performance['demand'], RGB_Futures_performance['max duration_failure'], 'r-', label='RGB Futures')
            ax[1,0].plot(PCR_performance['demand'], PCR_performance['max duration_failure'], 'b--', label='PCR GLOBWB2')

            ax[1,1].plot(OBS_performance['demand'], OBS_performance['avg_supply_gap(by episode)'], 'ko-', label='Observed')
            ax[1,1].plot(RGB_Futures_performance['demand'], RGB_Futures_performance['avg_supply_gap(by episode)'], 'r-', label='RGB Futures')
            ax[1,1].plot(PCR_performance['demand'], PCR_performance['avg_supply_gap(by episode)'], 'b--', label='PCR GLOBWB2')

            for axx in ax.flatten():
                axx.grid()
                
                axx.set_xlabel('Demand (ratio of avg observed flow)')
                

            ax[0,0].set_ylabel('Reliability')
            ax[0,1].set_ylabel('Mean duration failure (months)')
            ax[1,0].set_ylabel('Max duration failure (months)')
            ax[1,1].set_ylabel('Max supply gap (AF)')

            ax[0,1].legend()

            plt.suptitle('{}'.format(gage_id))
            plt.tight_layout()
            fig.savefig('../figures/sketch_datacenter_RioGrande/{}_performance.png'.format(gage_id))
            #plt.show()
            plt.close()

        performance.to_excel(writer, sheet_name='{}'.format(gage_id))






