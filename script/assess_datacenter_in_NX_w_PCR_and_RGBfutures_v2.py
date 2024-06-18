# Conda: pcr_globwb2_env

import numpy as np
import pandas as pd

import itertools 
import matplotlib.pyplot as plt

m3s_to_afd = 70.0452
cfs_to_afd = 1.98347
mgd_to_afd = 3.06888

def rrv(supply, targets):
    reliability = 1 - np.sum(supply <= targets) / len(targets)
    
    failure = np.where(supply <= targets, 1, 0)
    bool_failure = supply <=targets #supply[supply <= targets] == 1
    duration_failures = \
        [sum(g) for bool_failure, g in itertools.groupby(failure) if bool_failure]

    if len(duration_failures)==0:
        duration_failures = [0]        

    gap = np.where(supply <= targets, targets - supply, 0)
    bool_gap = supply <= targets
    supply_gap = \
        [sum(g) for bool_gap, g in itertools.groupby(gap,lambda x: x>0) if bool_gap]
    
    if len(supply_gap)==0: 
        supply_gap = [0]

    #gap_episodes = [supply_gap[i] * 30.5 for i in range(len(duration_failures))]
    gap_episodes = [supply_gap[i] * 30.5 for i in range(len(supply_gap))]

    
    #supply_gap = np.sum(np.where(supply <= targets, targets - supply, 0))
    #avg_supply_gap = supply_gap / len(duration_failures )

    return reliability, duration_failures, supply_gap, gap_episodes

gages = (
    ('RGB', 'USGS08353000', 34.41, -106.85, 'USGS_08353000.csv'),
    ('RGB', 'USGS08362500', 32.88, -107.29, '08362500_RIO_GRANDE_BLW_CABALLO_DAM.csv'),
    ('RGB', 'USGS08361000', 33.15, -107.21, '08361000_RIO_GRANDE_BELOW_ELEPHANT_BUTTE_DAM_fill.csv'),
    ('RGB', 'USGS08317400', 35.62, -106.32, 'USGS_08317400.csv')
)

demand_scaling_factor = (0.8, 0.9, 1, 1.1, 1.2)
hyperscale_data_center = 0.55 # MGD

discharge = pd.DataFrame()


with pd.ExcelWriter('../data/RGB_Futures/Performance_metrics_hyperscale_datacenter_RioGrande.xlsx') as writer:

    for i, (gage_region, gage_id, gage_lat, gage_lon, name_file) in enumerate(gages):
        
            
        if gage_id == 'USGS08362500':
            # Common period for the observed and simulated data (PCR and TOVA)
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
            
            # Common period for the observed and simulated data (PCR and TOVA)
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

            # Common period for the observed and simulated data (PCR and TOVA)
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

        elif gage_id == 'USGS08353000':
            # Common period for the observed and simulated data (PCR and TOVA)
            start = '1958-11-01'
            end = '2015-12-31'

            gage = pd.read_csv('../data/RGB_Futures/{}'.format(name_file),
                index_col='Dates', usecols=['Dates', 'Flow_cfs']) * cfs_to_afd
            gage.index = pd.to_datetime(gage.index)
            gage = gage.resample('ME').mean()
            gage = gage.truncate(before=start, after=end)
            gage.columns = ['gauges']

            PCR = pd.read_csv('../data/RGB_Futures/PCR_GLOBWB2_{}.csv'.format(gage_id), index_col='Date',
                            parse_dates=True)
            PCR = PCR.resample('ME').mean()
            PCR = PCR.truncate(before=start, after=end)




        performance = pd.DataFrame(columns=['demand (MGD)','reliability', 'mean duration_failure', 
                                            'max duration_failure', 'max_supply_gap(by episode, AF)'])

        for k, scaling_factor in enumerate(demand_scaling_factor, start=1):
            

            demand = np.ones(gage['gauges'].shape) * scaling_factor * hyperscale_data_center * mgd_to_afd

            fig, ax = plt.subplots(3, 1, figsize=(8,12))

            ax[0].plot(gage.index, gage['gauges']+0.0000000000001, 'k', label='Observed')
            ax[0].plot(gage.index, demand, 'r', label='Demand')
            ax[0].fill_between(gage.index, gage['gauges'].values, demand, gage['gauges'].values<=demand,
                            color='red', alpha=0.5)

            if gage_id != 'USGS08353000':
                ax[1].plot(RGB_Futures.index, RGB_Futures.sum(axis=1)+0.0000000000001, 'k', label='TOVA')
                ax[1].plot(RGB_Futures.index, demand, 'r', label='Demand')
                ax[1].fill_between(RGB_Futures.index, RGB_Futures.sum(axis=1).values, demand, RGB_Futures.sum(axis=1).values <=demand,
                                    color='red', alpha=0.5)
            
            ax[2].plot(PCR.index, PCR['PCR_GLOBWB2']+0.0000000000001, 'k', label='Competitor')
            ax[2].plot(PCR.index, demand, 'r', label='Demand')
            ax[2].fill_between(PCR.index, PCR['PCR_GLOBWB2'].values, demand, PCR['PCR_GLOBWB2'].values <=demand,
                                color='red', alpha=0.5)
            
            ax[0].set_title('Observed')
            if gage_id != 'USGS08353000':
                ax[1].set_title('TOVA')
            ax[2].set_title('Competitor')

            for axx in ax:
                axx.legend()
                #axx.grid()
                axx.set_ylabel('Water availabilty\n(Acre Feet/day)')
                axx.set_ylim(bottom=0.2)
                axx.set_yscale('log')

            plt.suptitle('{} -- Demand = {} MGD'.format(gage_id, np.round(hyperscale_data_center*scaling_factor,2)))

            plt.tight_layout()
            fig.savefig('../figures/sketch_hyperscale_datacenter_RioGrande/{}_{}.png'.format(gage_id, scaling_factor))
            #plt.show()
            plt.close()

            
            # Calculate the supply
            rel_PCR, dur_PCR, ss, supply_gap_PCR = rrv(PCR['PCR_GLOBWB2'].values, demand)
            if gage_id != 'USGS08353000':
                rel_RBG, dur_RBG, ss, supply_gap_RBG = rrv(RGB_Futures.sum(axis=1).values, demand)
            rel_OBS, dur_OBS, ss, supply_gap_OBS = rrv(gage.values.ravel(), demand)

            performance.loc['OBS_{}'.format(k)] = \
                [scaling_factor*hyperscale_data_center, rel_OBS, np.mean(dur_OBS), np.max(dur_OBS), np.max(supply_gap_OBS)]
            if gage_id != 'USGS08353000':
                performance.loc['RGB_Futures_{}'.format(k)] = \
                    [scaling_factor*hyperscale_data_center, rel_RBG, np.mean(dur_RBG), np.max(dur_RBG), np.max(supply_gap_RBG)]
            performance.loc['PCR_{}'.format(k)] = \
                [scaling_factor*hyperscale_data_center, rel_PCR, np.mean(dur_PCR), np.max(dur_PCR), np.max(supply_gap_PCR)]
            
            OBS_performance = performance.loc[performance.index.str.contains('OBS')]
            if gage_id != 'USGS08353000':
                RGB_Futures_performance = performance.loc[performance.index.str.contains('RGB_Futures')]
            PCR_performance = performance.loc[performance.index.str.contains('PCR')]

            fig, ax = plt.subplots(2, 2, figsize=(8,6))

            ax[0,0].plot(OBS_performance['demand (MGD)'], OBS_performance['reliability']*100, 'ko-', label='Observed')
            if gage_id != 'USGS08353000':
                ax[0,0].plot(RGB_Futures_performance['demand (MGD)'], RGB_Futures_performance['reliability'], 'r-', label='TOVA')
            ax[0,0].plot(PCR_performance['demand (MGD)'], PCR_performance['reliability']*100, 'b--', label='Competitor')

            ax[0,1].plot(OBS_performance['demand (MGD)'], OBS_performance['mean duration_failure'], 'ko-', label='Observed')
            if gage_id != 'USGS08353000':
                ax[0,1].plot(RGB_Futures_performance['demand (MGD)'], RGB_Futures_performance['mean duration_failure'], 'ro-', label='TOVA')
            ax[0,1].plot(PCR_performance['demand (MGD)'], PCR_performance['mean duration_failure'], 'b--', label='Competitor')

            ax[1,0].plot(OBS_performance['demand (MGD)'], OBS_performance['max duration_failure'], 'ko-', label='Observed')
            if gage_id != 'USGS08353000':
                ax[1,0].plot(RGB_Futures_performance['demand (MGD)'], RGB_Futures_performance['max duration_failure'], 'r-', label='TOVA')
            ax[1,0].plot(PCR_performance['demand (MGD)'], PCR_performance['max duration_failure'], 'b--', label='Competitor')

            ax[1,1].plot(OBS_performance['demand (MGD)'], OBS_performance['max_supply_gap(by episode, AF)'], 'ko-', label='Observed')
            if gage_id != 'USGS08353000':
                ax[1,1].plot(RGB_Futures_performance['demand (MGD)'], RGB_Futures_performance['max_supply_gap(by episode, AF)'], 'r-', label='TOVA')
            ax[1,1].plot(PCR_performance['demand (MGD)'], PCR_performance['max_supply_gap(by episode, AF)'], 'b--', label='Competitor')

            for axx in ax.flatten():
                axx.grid()
                
                axx.set_xlabel('Demand (x 0.55 MGD)')
                

            ax[0,0].set_ylabel('Reliability (%)')
            ax[0,1].set_ylabel('Mean duration failure (months)')
            ax[1,0].set_ylabel('Max duration failure (months)')
            ax[1,1].set_ylabel('Max cumul supply gap (Acre Feet)')

            ax[0,1].legend()

            plt.suptitle('{}'.format(gage_id))
            plt.tight_layout()
            fig.savefig('../figures/sketch_hyperscale_datacenter_RioGrande/{}_performance.png'.format(gage_id))
            #plt.show()
            plt.close()
        
        performance.to_excel(writer, sheet_name='{}'.format(gage_id))






