import sys
import pandas as pd
import datetime as dt
import numpy as np


#%% ------------------------------------- ###
###            INPUT CHOICES              ###
### ------------------------------------- ###

# sys.argv.append(r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\Inputs.xlsx') #f_inputs = pd.read_excel('Inputs.xlsx')
# sys.argv.append(r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\Balmorel Results\csv\Base\PriceElectricityHourly_BaseBeforeNoEXP.csv')
# sys.argv.append(r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\Balmorel Results\csv\Base\PriceHeatHourly_BaseBeforeNoEXP.csv')
# sys.argv.append(r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\Balmorel Results\csv\Base\DemandHeatHourly_BaseBeforeNoEXP.csv')
# sys.argv.append(r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\Balmorel Results\csv\Base\FuelConsumptionYearly_BaseBeforeNoEXP.csv') #f_inputs = pd.read_excel('Inputs.xlsx')

f_inputs = pd.read_excel(sys.argv[1])

# Input times
years = f_inputs[f_inputs['parameter'] == 'years']['value'].iloc[0]
years = list(np.array(years.rstrip(']').lstrip('[').split(',')).astype(float))
years = [int(element) for element in years]
weeks = f_inputs[f_inputs['parameter'] == 'seasons']['value'].iloc[0]
weeks = list(np.array(weeks.rstrip(']').lstrip('[').split(',')).astype(float))
weeks = [int(element) for element in weeks] # Apparently needed
hours = f_inputs[f_inputs['parameter'] == 'hours']['value'].iloc[0]

# Scenario
# SC = 'Base'

#%% ------------------------------------- ###
###     Create Full Timeseries Indices    ###
### ------------------------------------- ###


### Rep period ratio, for correcting demand and prices
# Electricity and heat prices are factored by this to correct yearly costs,
# yearly demand (not hourly) is divided to correct for capacities
rat = 52*168 / (len(weeks)*hours) 

# Hourly resolution
t_res = int(168 / hours)


act_time = pd.to_datetime([])
for i in range(len(years)):
    for j in range(len(weeks)):
        
        # Hours to first week
        delta_start = dt.timedelta(days = (weeks[j]-1)*7) 
        
        # First hour
        start = dt.datetime(years[i], 1, 1, 0) + delta_start
        
        # Timeseries for actual week
        temp = pd.to_datetime([k for k in range(0, 168, t_res)], unit='h',
                              origin=pd.Timestamp(start))
        
        # Save to actual timeseries
        act_time = act_time.append(temp)


# Modelled time
mod_time = pd.to_datetime([k for k in range(0, len(weeks)*168*len(years), t_res)], unit='h',
                          origin=pd.Timestamp(dt.datetime(years[0], 1, 1)))

# Hourly resolution time
mod_time_hourly = pd.to_datetime([k for k in range(0, len(weeks)*168*len(years))], unit='h',
                                 origin=pd.Timestamp(dt.datetime(years[0], 1, 1)))



#%% ------------------------------------- ###
###             Create CSV's              ###
### ------------------------------------- ###

### How to: Copy + Paste object names into this script at respective parameters

# End date
end = mod_time[-1] + dt.timedelta(hours=t_res)


# Indices between years
points = len(years)
ind = int(len(mod_time) / points)

# Years + last year
year_pL = years + [years[-1]]


#%% --------------------- ###
###       Time Inputs     ###
### --------------------- ###
# file_single = pd.DataFrame([], columns=['object class name', 'object name', 'parameter names', 'alternative names',	
#                                  'parameter index', 'parameter values',	 
#                                  'unit'])

file = pd.DataFrame([], columns=['object class name', 'object name', 'parameter names', 'alternative names',	
                                 'parameter index', 'parameter values',	 
                                 'unit'])


file_rel_2D = pd.DataFrame([], columns=['relationship class name', 'object class name 1', 'object class name 2', 
                                        'object name 1', 'object name 2', 'parameter names', 'alternative names',	
                                        'parameter index', 'parameter values',	 
                                        'unit'])

# Model start and end
# file = file.append(pd.DataFrame([['model', 'simple', 'model_start', SC, mod_time[0], 1, 'date']], 
#                                 columns=file.columns), ignore_index=True)
#file = file.append(pd.DataFrame([['model', 'simple', 'model_end', SC, end, str(end).replace(' ', 'T'), '#']], 
#                                columns=file.columns), ignore_index=True)
#{"data": "2030-01-15T00:00:00"}::date_time


#%% --------------------- ###
###  Storages etc. at 0   ###
### --------------------- ###

# All storage-like nodes, except biomass potentials

# COPY+PASTE HERE
nodes = """DK1_CO2_STO
DK1_DistrictHeating
DK1_H2_STO_SMA
DK2_CO2_STO
DK2_DistrictHeating
DK2_H2_STO_SMA
EXCESS_Road
NO1_CO2_STO
NO1_DistrictHeating
NO1_H2_STO_SMA
NO2_CO2_STO
NO2_DistrictHeating
NO2_H2_STO_SMA
NO3_CO2_STO
NO3_DistrictHeating
NO3_H2_STO_SMA
NO4_CO2_STO
NO4_DistrictHeating
NO4_H2_STO_SMA
NO5_CO2_STO
NO5_DistrictHeating
NO5_H2_STO_SMA
SE1_CO2_STO
SE1_DistrictHeating
SE1_H2_STO_SMA
SE2_CO2_STO
SE2_DistrictHeating
SE2_H2_STO_SMA
SE3_CO2_STO
SE3_DistrictHeating
SE3_H2_STO_SMA
SE4_CO2_STO
SE4_DistrictHeating
SE4_H2_STO_SMA""".split('\n')

vals = [0, 'NaN']
index = [mod_time[0]-dt.timedelta(hours=t_res), mod_time[0]]
SC = 'Base'
# Create data
for node in nodes:
    for i in range(2):   
        file = file.append(pd.DataFrame([['node', node, 
                                          'fix_node_state', SC, 
                                          index[i], vals[i], 
                                          'MWh']], 
                                        columns=file.columns), ignore_index=True)
        

             
#%% --------------------- ###
### Biomass Availability  ###
### --------------------- ###

f_fuel = pd.read_csv(sys.argv[5])
idxF = (f_fuel['FFF'] == 'WOODCHIPS') | (f_fuel['FFF'] == 'STRAW') 
idxY = f_fuel['Y'] == years[0]
idxC = (f_fuel['C'] == 'DENMARK') | (f_fuel['C'] == 'NORWAY') | (f_fuel['C'] == 'SWEDEN')

# fuels = f_fuel[idxF & idxY & idxA].groupby(['Y', 'AAA', 'FFF', 'TECH_TYPE', 'UNITS'], as_index=False)
fuels = f_fuel[idxF & idxY & idxC].groupby(['Y', 'RRR','FFF','UNITS'], as_index=False) # Same emission factor for all
fuels = fuels.aggregate(np.sum)
fuels['FFF'] = fuels['FFF'].replace('STRAW', 'Straw')
fuels['FFF'] = fuels['FFF'].replace('WOODCHIPS', 'Wood')


cols = file_rel_2D.columns

index = []
for i in range(len(years)): 
    index.append(mod_time[i*ind]-dt.timedelta(hours=t_res))  
    index.append(mod_time[i*ind])   
index.append(end)


bio = f_inputs[f_inputs['category'] == 'biomass potential']


# Create data
for node in np.unique(bio['parameter']):
    
    for SC in np.unique(bio['scenario']):
        # Get node and scenario
        temp = bio[bio['parameter'] == node]
        temp = temp[temp['scenario'] == SC]
        
        # Values (Assuming no more than two years!)
        val1 = np.int32(temp['value'].iloc[0])
        val2 = np.int32(temp['value'].iloc[1])
        
        # print(node,val1)
        # Subtract Balmorel consumption
        try:
            val1_bal = round(fuels[(fuels['RRR'] == node[:3]) & (fuels['FFF'] == node[4:])]['Val'].values[0]*1e6)
        except IndexError:
            val1_bal = 0
        # print(val1_bal)
        vals = [val1, 'NaN', val2, 'NaN']
        pro_vals = ['NaN', 0, 'NaN', 0]
        
        for i in range(len(years)*2):   
            if type(vals[i])==str:
                val = vals[i] 
            else:
                val = vals[i] - val1_bal
                # val = vals[i]/rat
                val = np.int32(val)
                print('Biomass potential', node, val, 'MWh')
                
            file = file.append(pd.DataFrame([['node', node, 'fix_node_state', SC, index[i], val, 'MWh']], 
                                            columns=file.columns), ignore_index=True)
    
            # file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__to_node'], 
            #                                               cols[1] : ['unit'], 
            #                                               cols[2] : ['node'], 
            #                                               cols[3] : ['bio_infinite_producer'],
            #                                               cols[4] : [node], 
            #                                               cols[5] : ['fix_unit_flow'], 
            #                                               cols[6] : [SC], 
            #                                               cols[7] : [index[i]], 
            #                                               cols[8] : [pro_vals[i]], 
            #                                               cols[9] : ['MWh']}))
            
        # Last value
        # file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__to_node'], 
        #                                               cols[1] : ['unit'], 
        #                                               cols[2] : ['node'], 
        #                                               cols[3] : ['bio_infinite_producer'],
        #                                               cols[4] : [node], 
        #                                               cols[5] : ['fix_unit_flow'], 
        #                                               cols[6] : [SC], 
        #                                               cols[7] : [end], 
        #                                               cols[8] : [pro_vals[i]], 
        #                                               cols[9] : ['MWh']}))
    






#%% --------------------- ###
###      Fuel demand      ###
### --------------------- ###

# """
# Assuming hourly demands for a representative year
# """

# dems = f_inputs[f_inputs['category'] == 'demand']

# index = []
# index_act = []
# for i in range(len(years)):
#     index.append(mod_time[i*ind])
#     index_act.append(act_time[i*ind])
    
# # Enddates
# index.append(end)
# index_act.append(act_time[i*ind]+dt.timedelta(hours=t_res))


# for i in range(len(np.unique(dems['parameter']))): 
    
#     for SC in np.unique(dems['scenario']):
#         # Select tech
#         temp = dems[dems['parameter'] == np.unique(dems['parameter'])[i]] 

#         # Select scenario
#         temp = temp[temp['scenario'] == SC]
#         for j in range(len(index)):
            
#             if len(temp) != 0:
#                 idx = temp['index'] == index_act[j].year
            
#                 name = temp['parameter'][idx].iloc[0]
#                 val = temp['value'][idx].iloc[0]
        
                
#                 file = file.append(pd.DataFrame([['node', name, 'demand', 
#                                                   'Base', index[j], 
#                                                   val, 'MWh']], 
#                                                 columns=file.columns), ignore_index=True)
    
#%% --------------------- ###
###      Yrly Demand      ###
### --------------------- ###

# Hourly demands 
dems = f_inputs[f_inputs['category'] == 'demand']

# Time index
index = []
index_act = []
for i in range(len(years)):
    index.append(mod_time[i*ind])
    index_act.append(act_time[i*ind])
index = [mod_time[0]-dt.timedelta(hours=t_res), mod_time[0], mod_time[-1]]


# Filter for year
idx = dems['index'] == index_act[0].year


dems_filt = dems[idx]
# Create data
for i in range(len(dems_filt)):
    node = dems_filt['parameter'].values[i]
    # print(node)
    
    # Extract value
    # val = dems_filt['value'].values[i]*8736 / rat # THIS SHOULD BE CORRECTED TO A FULL YEAR
    if node.find('H2') == -1:
        # print(node)
        val = dems_filt['value'].values[i]*8736 # THIS SHOULD BE CORRECTED TO A FULL YEAR
        val = np.int32(val)
        
        vals = [0, 'NaN', val]
        
        # Extract scenario
        SC = dems_filt['scenario'].values[i]
        # print(val*len(weeks)*hours)
        for j in range(3):   
            file = file.append(pd.DataFrame([['node', node, 
                                              'fix_node_state', SC, 
                                              index[j], vals[j], 
                                              'MWh']], 
                                            columns=file.columns), ignore_index=True)
    else:
        # print(node)
        val = dems_filt['value'].values[i] * 13 # Hourly demand assumed (Flow x 13)
        # val = np.int32(val)
        # print(dems_filt['value'].values[i])
        vals = [val, val]
        # print(val)
        # Extract scenario
        SC = dems_filt['scenario'].values[i]
        # print(val*len(weeks)*hours)
        for j in range(2):    
            file = file.append(pd.DataFrame([['node', node, 
                                              'demand', SC, 
                                              index[j+1], vals[j], 
                                              'MWh']], 
                                            columns=file.columns), ignore_index=True)
            
        
#%% --------------------- ###
###        El Prices      ###
### --------------------- ###


# COPY+PASTE HERE
rels = pd.Series("""DK1_El,DK1
DK2_El,DK2
NO1_El,NO1
NO2_El,NO2
NO3_El,NO3
NO4_El,NO4
NO5_El,NO5
SE1_El,SE1
SE2_El,SE2
SE3_El,SE3
SE4_El,SE4
DK1_SellingElectricity,DK1
DK2_SellingElectricity,DK2
NO1_SellingElectricity,NO1
NO2_SellingElectricity,NO2
NO3_SellingElectricity,NO3
NO4_SellingElectricity,NO4
NO5_SellingElectricity,NO5
SE1_SellingElectricity,SE1
SE2_SellingElectricity,SE2
SE3_SellingElectricity,SE3
SE4_SellingElectricity,SE4""".split('\n')).str.split(',')

# Heat areas - HARDCODED FOR NOW - maybe make dictionary for ESB, CPH and so on?
# AAA = """DK1_Large
# DK2_Medium
# NO1_A2""".split('\n')

# path = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure' + '/Balmorel Results/csv/'

# f_elPrice = pd.read_csv(path+SC+'/PriceElectricityHourly_'+SC+'.csv')
f_elPrice = pd.read_csv(sys.argv[2])

SC = 'Base'
tariffs = f_inputs[f_inputs['category'] == 'tariff']

# Make Balmorel Weeks
BalmWeeks = pd.Series(['S' + str(Bweek) for Bweek in weeks])
idx = BalmWeeks.str.len() < 3
BalmWeeks[idx] = BalmWeeks[idx].str.replace('S', 'S0')


# Make electricity prices
# for SC in np.unique(f_elPrice['SC']):
#     idxel = f_elPrice['SC'] == SC
for i in range(len(years)):
    idxel0 = f_elPrice['Y'] == years[i]
    for j in range(len(BalmWeeks)):
        index = mod_time[i*ind+j*hours:i*ind+j*hours + hours]    
        idxel1 = f_elPrice['SSS'] == BalmWeeks[j]
        for r in range(len(rels)):
            idxel2 = f_elPrice['RRR'] == rels[r][1]
            
            val = f_elPrice[idxel0 & idxel1 & idxel2]['Val'].iloc[::t_res].round(1)
            
            # Include tariff
            tar = tariffs[(tariffs['parameter'] == rels[r][1]) & (tariffs['index'] == years[i])]['value'].values[0]
            
            
            
            if rels[r][0][-2:] == 'El':
                val = val + tar # Production requires tariff
                
                cols = file_rel_2D.columns
                file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__from_node']*len(index), 
                                                              cols[1] : ['unit']*len(index), 
                                                              cols[2] : ['node']*len(index), 
                                                              cols[3] : [rels[r][0]]*len(index), 
                                                              cols[4] : ['ElProduced']*len(index), 
                                                              cols[5] : ['fuel_cost']*len(index), 
                                                              cols[6] : [SC]*len(index), 
                                                              cols[7] : index, 
                                                              cols[8] : val, 
                                                              cols[9] : ['euro/MWh']*len(index)}),
                                                  ignore_index=True)
            else:
                
                cols = file_rel_2D.columns
                file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__from_node']*len(index), 
                                                              cols[1] : ['unit']*len(index), 
                                                              cols[2] : ['node']*len(index), 
                                                              cols[3] : [rels[r][0]]*len(index), 
                                                              cols[4] : [rels[r][1]]*len(index), 
                                                              cols[5] : ['fuel_cost']*len(index), 
                                                              cols[6] : [SC]*len(index), 
                                                              cols[7] : index, 
                                                              cols[8] : -val, 
                                                              cols[9] : ['euro/MWh']*len(index)}),
                                                  ignore_index=True)
                

                
#%% --------------------- ###
###     Biomass Prices    ###
### --------------------- ###

### Not a time series! So doesn't work

# # COPY+PASTE HERE
# rels = pd.Series("""DK1_ST_HTL,DK1_Straw
# DK1_ST_Pyro,DK1_Straw
# DK1_ST_TG-FT,DK1_Straw
# DK1_ST_TG-FT-H2,DK1_Straw
# DK1_ST_TG-MeOH,DK1_Straw
# DK1_WO_HTL,DK1_Wood
# DK1_WO_Pyro,DK1_Wood
# DK1_WO_TG-FT,DK1_Wood
# DK1_WO_TG-FT-H2,DK1_Wood
# DK1_WO_TG-MeOH,DK1_Wood
# DK2_ST_HTL,DK2_Straw
# DK2_ST_Pyro,DK2_Straw
# DK2_ST_TG-FT,DK2_Straw
# DK2_ST_TG-FT-H2,DK2_Straw
# DK2_ST_TG-MeOH,DK2_Straw
# DK2_WO_HTL,DK2_Wood
# DK2_WO_Pyro,DK2_Wood
# DK2_WO_TG-FT,DK2_Wood
# DK2_WO_TG-FT-H2,DK2_Wood
# DK2_WO_TG-MeOH,DK2_Wood
# NO1_ST_HTL,NO1_Straw
# NO1_ST_Pyro,NO1_Straw
# NO1_ST_TG-FT,NO1_Straw
# NO1_ST_TG-FT-H2,NO1_Straw
# NO1_ST_TG-MeOH,NO1_Straw
# NO1_WO_HTL,NO1_Wood
# NO1_WO_Pyro,NO1_Wood
# NO1_WO_TG-FT,NO1_Wood
# NO1_WO_TG-FT-H2,NO1_Wood
# NO1_WO_TG-MeOH,NO1_Wood
# NO2_ST_HTL,NO2_Straw
# NO2_ST_Pyro,NO2_Straw
# NO2_ST_TG-FT,NO2_Straw
# NO2_ST_TG-FT-H2,NO2_Straw
# NO2_ST_TG-MeOH,NO2_Straw
# NO2_WO_HTL,NO2_Wood
# NO2_WO_Pyro,NO2_Wood
# NO2_WO_TG-FT,NO2_Wood
# NO2_WO_TG-FT-H2,NO2_Wood
# NO2_WO_TG-MeOH,NO2_Wood
# NO3_ST_HTL,NO3_Straw
# NO3_ST_Pyro,NO3_Straw
# NO3_ST_TG-FT,NO3_Straw
# NO3_ST_TG-FT-H2,NO3_Straw
# NO3_ST_TG-MeOH,NO3_Straw
# NO3_WO_HTL,NO3_Wood
# NO3_WO_Pyro,NO3_Wood
# NO3_WO_TG-FT,NO3_Wood
# NO3_WO_TG-FT-H2,NO3_Wood
# NO3_WO_TG-MeOH,NO3_Wood
# NO4_ST_HTL,NO4_Straw
# NO4_ST_Pyro,NO4_Straw
# NO4_ST_TG-FT,NO4_Straw
# NO4_ST_TG-FT-H2,NO4_Straw
# NO4_ST_TG-MeOH,NO4_Straw
# NO4_WO_HTL,NO4_Wood
# NO4_WO_Pyro,NO4_Wood
# NO4_WO_TG-FT,NO4_Wood
# NO4_WO_TG-FT-H2,NO4_Wood
# NO4_WO_TG-MeOH,NO4_Wood
# NO5_ST_HTL,NO5_Straw
# NO5_ST_Pyro,NO5_Straw
# NO5_ST_TG-FT,NO5_Straw
# NO5_ST_TG-FT-H2,NO5_Straw
# NO5_ST_TG-MeOH,NO5_Straw
# NO5_WO_HTL,NO5_Wood
# NO5_WO_Pyro,NO5_Wood
# NO5_WO_TG-FT,NO5_Wood
# NO5_WO_TG-FT-H2,NO5_Wood
# NO5_WO_TG-MeOH,NO5_Wood
# SE1_ST_HTL,SE1_Straw
# SE1_ST_Pyro,SE1_Straw
# SE1_ST_TG-FT,SE1_Straw
# SE1_ST_TG-FT-H2,SE1_Straw
# SE1_ST_TG-MeOH,SE1_Straw
# SE1_WO_HTL,SE1_Wood
# SE1_WO_Pyro,SE1_Wood
# SE1_WO_TG-FT,SE1_Wood
# SE1_WO_TG-FT-H2,SE1_Wood
# SE1_WO_TG-MeOH,SE1_Wood
# SE2_ST_HTL,SE2_Straw
# SE2_ST_Pyro,SE2_Straw
# SE2_ST_TG-FT,SE2_Straw
# SE2_ST_TG-FT-H2,SE2_Straw
# SE2_ST_TG-MeOH,SE2_Straw
# SE2_WO_HTL,SE2_Wood
# SE2_WO_Pyro,SE2_Wood
# SE2_WO_TG-FT,SE2_Wood
# SE2_WO_TG-FT-H2,SE2_Wood
# SE2_WO_TG-MeOH,SE2_Wood
# SE3_ST_HTL,SE3_Straw
# SE3_ST_Pyro,SE3_Straw
# SE3_ST_TG-FT,SE3_Straw
# SE3_ST_TG-FT-H2,SE3_Straw
# SE3_ST_TG-MeOH,SE3_Straw
# SE3_WO_HTL,SE3_Wood
# SE3_WO_Pyro,SE3_Wood
# SE3_WO_TG-FT,SE3_Wood
# SE3_WO_TG-FT-H2,SE3_Wood
# SE3_WO_TG-MeOH,SE3_Wood
# SE4_ST_HTL,SE4_Straw
# SE4_ST_Pyro,SE4_Straw
# SE4_ST_TG-FT,SE4_Straw
# SE4_ST_TG-FT-H2,SE4_Straw
# SE4_ST_TG-MeOH,SE4_Straw
# SE4_WO_HTL,SE4_Wood
# SE4_WO_Pyro,SE4_Wood
# SE4_WO_TG-FT,SE4_Wood
# SE4_WO_TG-FT-H2,SE4_Wood
# SE4_WO_TG-MeOH,SE4_Wood""".split('\n')).str.split(',')

# # Heat areas - HARDCODED FOR NOW - maybe make dictionary for ESB, CPH and so on?
# # AAA = """DK1_Large
# # DK2_Medium
# # NO1_A2""".split('\n')

# #path = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure'

# #f_BioP = pd.read_csv(path+'Inputs.csv')

# bio_p = f_inputs[f_inputs['category'] == 'biomass price']

# # Make biomass price
# cols = file_rel_2D.columns
# for r in range(len(rels)):
#     idx = bio_p['parameter'] == rels[r][1].split('_')[1]
    
#     file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__from_node']*len(bio_p[idx]), 
#                                                   cols[1] : ['unit']*len(bio_p[idx]), 
#                                                   cols[2] : ['node']*len(bio_p[idx]), 
#                                                   cols[3] : [rels[r][0]]*len(bio_p[idx]), 
#                                                   cols[4] : [rels[r][1]]*len(bio_p[idx]), 
#                                                   cols[5] : ['fuel_cost']*len(bio_p[idx]), 
#                                                   cols[6] : bio_p[idx]['scenario'], 
#                                                   cols[7] : bio_p[idx]['index'], 
#                                                   cols[8] : bio_p[idx]['value'], 
#                                                   cols[9] : ['euro/MWh']*len(bio_p[idx])}), ignore_index=True)

#%% --------------------- ###
###      Heat Prices      ###
### --------------------- ###

# COPY+PASTE HERE
rels = pd.Series("""DK1_SellingWasteheat,DK1_PotWasteheat
DK2_SellingWasteheat,DK2_PotWasteheat
SE4_SellingWasteheat,SE4_PotWasteheat
NO2_SellingWasteheat,NO2_PotWasteheat
NO1_SellingWasteheat,NO1_PotWasteheat
SE3_SellingWasteheat,SE3_PotWasteheat
NO3_SellingWasteheat,NO3_PotWasteheat
NO4_SellingWasteheat,NO4_PotWasteheat
NO5_SellingWasteheat,NO5_PotWasteheat
SE2_SellingWasteheat,SE2_PotWasteheat
SE1_SellingWasteheat,SE1_PotWasteheat""".split('\n')).str.split(',')

# Heat areas - HARDCODED FOR NOW - maybe make dictionary for ESB, CPH and so on?
AAA = pd.Series("""DK1_Large
DK2_Large
SE4_large
NO2_A2
NO1_A3
SE3_large
NO3_A3
NO4_A2
NO5_A2
SE2_medium
SE1_medium""".split('\n')).str.split('_')

# f_HPrice = pd.read_csv(path+SC+'/PriceHeatHourly_'+SC+'.csv')
f_HPrice = pd.read_csv(sys.argv[3])

# Make heat prices
# cols = file.columns
cols = file_rel_2D.columns
# for SC in np.unique(f_HPrice['SC']):
#     idxH = f_HPrice['SC'] == SC
for i in range(len(years)):
    idxH0 = f_HPrice['Y'] == years[i]
    for j in range(len(BalmWeeks)):
        index = mod_time[i*ind+j*hours:i*ind+j*hours + hours]  
        idxH1 = f_HPrice['SSS'] == BalmWeeks[j]
        for r in range(len(rels)):
            idxH2 = f_HPrice['RRR'] == AAA[r][0]
            idxH3 = f_HPrice['AAA'] == AAA.str.join('_')[r]
            
            vals = f_HPrice[idxH0 & idxH1 & idxH2 & idxH3]['Val'].iloc[::t_res].round(1)
            # vals = f_HPrice[idxH0 & idxH1 & idxH2 & idxH3]['Val'].iloc[::t_res].round(1) / (1 + 1/0.99) # For connection_flow_cost
            
            
            # print()
            # print(index)
            # print(units[r])
            # print(units[r].replace('Heat', 'PotWasteheat'))
            if len(vals[vals <= 0]) >= 1:
                # print(vals[vals < 0])
                vals[vals <= 0] = 1
                # print(vals[vals < 0])
                print(np.unique(vals))
            
            # file = file.append(pd.DataFrame({cols[0] : ['connection']*len(index), 
            #                                               cols[1] : [conns[r]]*len(index), 
            #                                               cols[2] : ['connection_flow_cost']*len(index), 
            #                                               cols[3] : [SC]*len(index), 
            #                                               cols[4] : index, 
            #                                               cols[5] : -vals, 
            #                                               cols[6] : ['euro/MWh']*len(index)}),
            #                     ignore_index=True)

            file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__from_node']*len(index), 
                                                          cols[1] : ['unit']*len(index), 
                                                          cols[2] : ['node']*len(index), 
                                                          cols[3] : [rels[r][0]]*len(index), 
                                                          cols[4] : [rels[r][1]]*len(index), 
                                                          cols[5] : ['fuel_cost']*len(index), 
                                                          cols[6] : [SC]*len(index), 
                                                          cols[7] : index, 
                                                          cols[8] : -vals, 
                                                          cols[9] : ['euro/MWh']*len(index)}))
            
#%% --------------------- ###
###      Heat Demand      ###
### --------------------- ###

### connection_capacity

# COPY+PASTE HERE
rels = pd.Series("""DK1_SellingWasteheat,DK1_DistrictHeating
DK2_SellingWasteheat,DK2_DistrictHeating
SE4_SellingWasteheat,SE4_DistrictHeating
NO2_SellingWasteheat,NO2_DistrictHeating
NO1_SellingWasteheat,NO1_DistrictHeating
SE3_SellingWasteheat,SE3_DistrictHeating
NO3_SellingWasteheat,NO3_DistrictHeating
NO4_SellingWasteheat,NO4_DistrictHeating
NO5_SellingWasteheat,NO5_DistrictHeating
SE2_SellingWasteheat,SE2_DistrictHeating
SE1_SellingWasteheat,SE1_DistrictHeating""".split('\n')).str.split(',')

# Heat areas - HARDCODED FOR NOW - maybe make dictionary for ESB, CPH and so on?
AAA = pd.Series("""DK1_Large
DK2_Large
SE4_large
NO2_A2
NO1_A3
SE3_large
NO3_A3
NO4_A2
NO5_A2
SE2_medium
SE1_medium""".split('\n')).str.split('_')


# path = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure' + '/Balmorel Results/csv/'

# f_Hdem = pd.read_csv(path+SC+'/DemandHeatHourly_'+SC+'.csv')
f_Hdem = pd.read_csv(sys.argv[4])

# fnew = f_Hdem.groupby(['Y', 'C', 'RRR', 'AAA', 'SSS', 'TTT', 'UNITS', 'SC'], as_index=False)
fnew = f_Hdem[f_Hdem['VARIABLE_CATEGORY'] == 'EXOGENOUS']
fnew = fnew.groupby(['Y', 'C', 'RRR', 'AAA', 'SSS', 'TTT', 'UNITS'], as_index=False)
fnew = fnew.aggregate(np.sum)
# pd.pivot_table(f_Hdem, columns=['Y', 'C', 'RRR', 'AAA', 'SSS', 'UNITS'], index='TTT',
#          values='Val')

# Make heat demands
cols = file_rel_2D.columns
# for SC in np.unique(fnew['SC']):
#     idxH = fnew['SC'] == SC
for i in range(len(years)):
    idxH0 = fnew['Y'] == years[i]
    for j in range(len(BalmWeeks)):
        index = mod_time[i*ind+j*hours:i*ind+j*hours + hours]    
        idxH1 = fnew['SSS'] == BalmWeeks[j]
        for r in range(len(rels)):
            idxH2 = fnew['RRR'] == AAA[r][0]
            idxH3 = fnew['AAA'] == AAA.str.join('_')[r]
            
            # file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['connection__to_node']*len(index), 
            #                                               cols[1] : ['connection']*len(index), 
            #                                               cols[2] : ['node']*len(index), 
            #                                               cols[3] : [rels[r][0]]*len(index), 
            #                                               cols[4] : [rels[r][1]]*len(index), 
            #                                               cols[5] : ['connection_capacity']*len(index), 
            #                                               cols[6] : [SC]*len(index), 
            #                                               cols[7] : index, 
            #                                               cols[8] : fnew[idxH0 & idxH1 & idxH2 & idxH3]['Val'].iloc[::t_res].round(1), 
            #                                               cols[9] : ['euro/MWh']*len(index)}))

            file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__to_node']*len(index), 
                                                          cols[1] : ['unit']*len(index), 
                                                          cols[2] : ['node']*len(index), 
                                                          cols[3] : [rels[r][0]]*len(index), 
                                                          cols[4] : [rels[r][1]]*len(index), 
                                                          cols[5] : ['unit_capacity']*len(index), 
                                                          cols[6] : [SC]*len(index), 
                                                          cols[7] : index, 
                                                          cols[8] : fnew[idxH0 & idxH1 & idxH2 & idxH3]['Val'].iloc[::t_res].round(1), 
                                                          cols[9] : ['euro/MWh']*len(index)}))

#%% --------------------- ###
### Heat Production Price ###
### --------------------- ###

# COPY+PASTE HERE
units = """DK1_Heat
DK2_Heat
SE4_Heat
NO2_Heat
NO1_Heat
SE3_Heat
NO3_Heat
NO4_Heat
NO5_Heat
SE2_Heat
SE1_Heat""".split('\n')

# Heat areas - HARDCODED FOR NOW - maybe make dictionary for ESB, CPH and so on?
AAA = pd.Series("""DK1_Large
DK2_Large
SE4_large
NO2_A2
NO1_A3
SE3_large
NO3_A3
NO4_A2
NO5_A2
SE2_medium
SE1_medium""".split('\n')).str.split('_')

# f_HPrice = pd.read_csv(path+SC+'/PriceHeatHourly_'+SC+'.csv')
f_HPrice = pd.read_csv(sys.argv[3])

# Make heat prices
cols = file_rel_2D.columns
# for SC in np.unique(f_HPrice['SC']):
#     idxH = f_HPrice['SC'] == SC
for i in range(len(years)):
    idxH0 = f_HPrice['Y'] == years[i]
    for j in range(len(BalmWeeks)):
        index = mod_time[i*ind+j*hours:i*ind+j*hours + hours]  
        idxH1 = f_HPrice['SSS'] == BalmWeeks[j]
        for r in range(len(units)):
            idxH2 = f_HPrice['RRR'] == AAA[r][0]
            idxH3 = f_HPrice['AAA'] == AAA.str.join('_')[r]
            
            # file = file.append(pd.DataFrame({cols[0] : ['unit']*len(index), 
            #                                               cols[1] : [units[r]]*len(index), 
            #                                               cols[2] : ['fuel_cost']*len(index), 
            #                                               cols[3] : [SC]*len(index), 
            #                                               cols[4] : index, 
            #                                               cols[5] : f_HPrice[idxH0 & idxH1 & idxH2 & idxH3]['Val'].iloc[::t_res].round(1), 
            #                                               cols[6] : ['euro/MWh']*len(index)}),
            #                     ignore_index=True)
                
            # Get prices and set negative to zero
            vals = f_HPrice[idxH0 & idxH1 & idxH2 & idxH3]['Val'].iloc[::t_res].round(1) # Heat prices from Balmorel areas
            # vals = np.array([25]*len(index)) 
            # print()
            # print(index)
            # print(units[r])
            # print(units[r].replace('Heat', 'PotWasteheat'))
            if len(vals[vals <= 25]) >= 1:
                # print(vals[vals < 0])
                vals[vals <= 25] = 25 # IF HEAT IS PRODUCED, IT IS HIGH TEMPERATURE HEAT 25 €/MWh from Chmielarz2021
                
                # print(vals[vals < 0])
            

            file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__from_node']*len(index), 
                                                          cols[1] : ['unit']*len(index), 
                                                          cols[2] : ['node']*len(index), 
                                                          cols[3] : [units[r]]*len(index), 
                                                          cols[4] : ['HeatProduced']*len(index), 
                                                          cols[5] : ['fuel_cost']*len(index), 
                                                          cols[6] : [SC]*len(index), 
                                                          cols[7] : index, 
                                                          cols[8] : vals, 
                                                          cols[9] : ['euro/MWh']*len(index)}))
            # cols = file.columns
            # file = file.append(pd.DataFrame({cols[0] : ['connection']*len(index), 
            #                                               cols[1] : [units[r]]*len(index), 
            #                                               cols[2] : ['connection_flow_cost']*len(index), 
            #                                               cols[3] : [SC]*len(index), 
            #                                               cols[4] : index, 
            #                                               cols[5] : vals / (1 + 1/0.99), 
            #                                               cols[6] : ['euro/MWh']*len(index)}),
            #                    ignore_index=True)


#%% --------------------- ###
###     CO2 Potentials    ###
### --------------------- ###
f_fuel = pd.read_csv(sys.argv[5])

AAA = pd.Series("""DK1_Large
DK2_Large
SE4_large
NO2_A2
NO1_A3
SE3_large
NO3_A3
NO4_A2
NO5_A2
SE2_medium
SE1_medium""".split('\n'))


idxF = (f_fuel['FFF'] == 'WOODCHIPS') | (f_fuel['FFF'] == 'MUNIWASTE') | (f_fuel['FFF'] == 'WOODWASTE') | (f_fuel['FFF'] == 'STRAW') | (f_fuel['FFF'] == 'BIOGAS') | (f_fuel['FFF'] == 'WOODPELLETS') | (f_fuel['FFF'] == 'WOOD')
idxY = f_fuel['Y'] == years[0]
# idxY = f_fuel['Y'] == f_fuel['Y']
idxA = f_fuel['AAA'] != f_fuel['AAA']
for A in AAA:
    idxA = idxA | (f_fuel['AAA'] == A)

# fuels = f_fuel[idxF & idxY & idxA].groupby(['Y', 'AAA', 'FFF', 'TECH_TYPE', 'UNITS'], as_index=False)
fuels = f_fuel[idxF & idxY & idxA].groupby(['Y', 'RRR','UNITS'], as_index=False) # Same emission factor for all
fuels = fuels.aggregate(np.sum)

# COPY+PASTE HERE
nodes = """DK1_CO2_PS
DK2_CO2_PS
SE4_CO2_PS
NO2_CO2_PS
NO1_CO2_PS
SE3_CO2_PS
NO3_CO2_PS
NO4_CO2_PS
NO5_CO2_PS
SE2_CO2_PS
SE1_CO2_PS""".split('\n')

index = [mod_time[0]-dt.timedelta(hours=t_res), mod_time[0]]
SC = 'Base'
# Create data
for node in nodes:
    # Get yearly fuel consumption if possible
    try:
        val = fuels[fuels['RRR'] == node[:3]]['Val'].values[0]
    except IndexError:
        val = 0
    
    EM = 94.6 # Emission factor in kg/GJ
    val = val * 3600000 * EM / 1000 # TWh to GJ to ton CO2
    vals = [val, 'NaN'] 
    
    # print(node)
    print('%s = %d tCO2'%(node, val))
    for i in range(2):    
        file = file.append(pd.DataFrame([['node', 
                                          node, 
                                          'fix_node_state', 
                                          SC, 
                                          index[i], 
                                          vals[i],
                                          'tCO2']], 
                                        columns=file.columns), ignore_index=True)



#%% --------------------- ###
###          Save         ###
### --------------------- ###

# Time series
file.to_csv('timeseries.csv', index=False)

# Single values
# file_single.to_csv('single_values.csv', index=False)

# 2D Relationships
file_rel_2D.to_csv('timeseries_rel_2D.csv', index=False)
