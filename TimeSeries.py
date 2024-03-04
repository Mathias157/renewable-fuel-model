import sys
import pandas as pd
import datetime as dt
import numpy as np
import json

#%% ------------------------------------- ###
###            INPUT CHOICES              ###
### ------------------------------------- ###

sys.argv.append(r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\Inputs.xlsx') #f_inputs = pd.read_excel('Inputs.xlsx')
sys.argv.append(r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\Balmorel Results\csv\Base\PriceElectricityHourly_ItalianCase.csv')
sys.argv.append(r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\Balmorel Results\csv\Base\PriceHeatHourly_ItalianCase.csv')
sys.argv.append(r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\Balmorel Results\csv\Base\DemandHeatHourly_ItalianCase.csv')
sys.argv.append(r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\Balmorel Results\csv\Base\FuelConsumptionYearly_ItalianCase.csv') #f_inputs = pd.read_excel('Inputs.xlsx')

f_inputs = pd.read_excel(sys.argv[1])

# Input times
years = f_inputs[f_inputs['parameter'] == 'years']['value'].iloc[0]
years = list(np.array(years.rstrip(']').lstrip('[').split(',')).astype(float))
years = [int(element) for element in years]
weeks = f_inputs[f_inputs['parameter'] == 'seasons']['value'].iloc[0]
weeks = list(np.array(weeks.rstrip(']').lstrip('[').split(',')).astype(float))
weeks = [int(element) for element in weeks] # Apparently needed
hours = f_inputs[f_inputs['parameter'] == 'hours']['value'].iloc[0]


### Geographical Scope
reg_all = f_inputs[f_inputs['parameter'] == 'all_regions']['value'].iloc[0]
reg_all = list(np.array(reg_all.rstrip(']').lstrip('[').split(',')).astype(str))
reg_all = [element.replace(' ', '') for element in reg_all]
reg_cho = f_inputs[f_inputs['parameter'] == 'regions']['value'].iloc[0]
reg_cho = list(np.array(reg_cho.rstrip(']').lstrip('[').split(',')).astype(str))
reg_cho = [element.replace(' ', '') for element in reg_cho]

# Balmorel regions
reg_map = f_inputs[f_inputs['parameter'] == 'balmorel_mapping']['value'].iloc[0]
reg_map = reg_map.replace(':', '" : "')
reg_map = reg_map.replace(',', '" , "')
reg_map = reg_map.replace('{', '{"')
reg_map = reg_map.replace('}', '"}')
reg_map = reg_map.replace(' ', '')
reg_map = json.loads(reg_map)
    
# Balmorel areas
aaa_map = f_inputs[f_inputs['parameter'] == 'balmorel_areas']['value'].iloc[0]
aaa_map = aaa_map.replace(':', '" : "')
aaa_map = aaa_map.replace(',', '" , "')
aaa_map = aaa_map.replace('{', '{"')
aaa_map = aaa_map.replace('}', '"}')
aaa_map = aaa_map.replace(' ', '')
aaa_map = json.loads(aaa_map)
    
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
nodes = """Z01_CO2_STO
Z01_DistrictHeating
Z01_H2_STO_SMA
Z02_CO2_STO
Z02_DistrictHeating
Z02_H2_STO_SMA
EXCESS_Road
Z03_CO2_STO
Z03_DistrictHeating
Z03_H2_STO_SMA
Z04_CO2_STO
Z04_DistrictHeating
Z04_H2_STO_SMA
Z05_CO2_STO
Z05_DistrictHeating
Z05_H2_STO_SMA
Z06_CO2_STO
Z06_DistrictHeating
Z06_H2_STO_SMA
Z07_CO2_STO
Z07_DistrictHeating
Z07_H2_STO_SMA
Z08_CO2_STO
Z08_DistrictHeating
Z08_H2_STO_SMA
Z09_CO2_STO
Z09_DistrictHeating
Z09_H2_STO_SMA
Z10_CO2_STO
Z10_DistrictHeating
Z10_H2_STO_SMA
Z11_CO2_STO
Z11_DistrictHeating
Z11_H2_STO_SMA""".split('\n')

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
        try:
            val1 = np.int32(temp['value'].iloc[0])
            val2 = np.int32(temp['value'].iloc[1])
        except ValueError:
            val1 = 0
            val2 = 0
            
        # print(node,val1)
        # Subtract Balmorel consumption
        try:
            val1_bal = round(fuels[(fuels['RRR'] == node[:3]) & (fuels['FFF'] == node[4:])]['Val'].values[0]*1e6)
        except IndexError:
            val1_bal = 0
        # print(val1_bal)
        vals = [val1, 'NaN', val2, 'NaN']
        pro_vals = ['NaN', 0, 'NaN', 0]
        
        
        ### Only add value if in chosen regions
        if node.split('_')[0] in reg_cho:
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
                
        ### If not chosen, set to zero
        else:
            file = file.append(pd.DataFrame([['node', node, 'fix_node_state', SC, index[0], 0, 'MWh']],
                                            columns=file.columns), ignore_index=True)
            file = file.append(pd.DataFrame([['node', node, 'fix_node_state', SC, index[2], 0, 'MWh']],
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
indexH2 = [mod_time[0]-dt.timedelta(hours=t_res), mod_time[0]+dt.timedelta(hours=24), mod_time[-1]]
# Day-time demand
H2Profile = pd.read_excel(r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\H2Profile.xlsx')
# indexH2 = []
# for i in range(28*3):
#     indexH2.append(mod_time[0]+dt.timedelta(hours=8*i))
# indexH2.append(mod_time[-1])    
    
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
        
        
        ### Only add value if chosen region NB: Set custom made demand nodes in if statement
        if (node.split('_')[0] in reg_cho):
            for j in range(3):   
                file = file.append(pd.DataFrame([['node', node, 
                                                  'fix_node_state', SC, 
                                                  index[j], vals[j], 
                                                  'MWh']], 
                                                columns=file.columns), ignore_index=True)
        else:
            for j in range(3):
                file = file.append(pd.DataFrame([['node', node, 
                                                  'fix_node_state', SC, 
                                                  index[j], 0, 
                                                  'MWh']], 
                                                columns=file.columns), ignore_index=True)

    else:
        # print(node)
        val = dems_filt['value'].values[i] * 13 * 1.03448275862069 # Hourly demand assumed (Flow x 13, accounting for 0 first 2 days)
        # val = np.int32(val)

        # print(dems_filt['value'].values[i])
        vals = [0, val * (672/648), val * (672/648)] # Accounting for 0 the first day
        # vals = [0, val * (672/(16*28)), val * (672/(16*28))]*28 + [0] # Accounting for 0 in the night
        # print(val)
        # Extract scenario
        SC = dems_filt['scenario'].values[i]
        # print(val*len(weeks)*hours)
        
        
        ### Only add value if chosen region NB: Set custom made demand nodes, such as EU_H2 and FIN_H2 here!!! 
        if (node.split('_')[0] in reg_cho) | (node == 'EU_H2') | (node == 'FIN_H2'):
            for j in range(len(mod_time)):    
                file = file.append(pd.DataFrame([['node', node, 
                                                  'demand', SC, 
                                                  mod_time[j], val * H2Profile.values[j,0], 
                                                  'MWh']], 
                                                columns=file.columns), ignore_index=True)
        
        else:
            mod_times = [0, -1]
            for j in range(2):
                file = file.append(pd.DataFrame([['node', node, 
                                                  'demand', 'Base', 
                                                  mod_time[mod_times[j]], 0, 
                                                  'MWh']], 
                                                columns=file.columns), ignore_index=True)
                    
#%% --------------------- ###
###        El Prices      ###
### --------------------- ###


# COPY+PASTE HERE
rels = pd.Series("""Z01_El,Z01
Z02_El,Z02
Z03_El,Z03
Z04_El,Z04
Z05_El,Z05
Z06_El,Z06
Z07_El,Z07
Z08_El,Z08
Z09_El,Z09
Z10_El,Z10
Z11_El,Z11
Z01_SellingElectricity,Z01
Z02_SellingElectricity,Z02
Z03_SellingElectricity,Z03
Z04_SellingElectricity,Z04
Z05_SellingElectricity,Z05
Z06_SellingElectricity,Z06
Z07_SellingElectricity,Z07
Z08_SellingElectricity,Z08
Z09_SellingElectricity,Z09
Z10_SellingElectricity,Z10
Z11_SellingElectricity,Z11""".split('\n')).str.split(',')

# Heat areas - HARDCODED FOR NOW - maybe make dictionary for ESB, CPH and so on?
# AAA = """Z01_Large
# Z02_Medium
# Z03_A2""".split('\n')

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
            idxel2 = f_elPrice['RRR'] == reg_map[rels[r][1]]
            
            val = f_elPrice[idxel0 & idxel1 & idxel2]['Val'].iloc[::t_res].round(1)
            
            # Include tariff
            tar = tariffs[(tariffs['parameter'] == rels[r][1]) & (tariffs['index'] == years[i])]['value'].values[0]
            
            
            
            if rels[r][0][-2:] == 'El':
                        
                ### Only add value if in chosen regions
                if rels[r][0].split('_')[0] in reg_cho:
                    val = val + tar # ion requires tariff
                    index2 = index
                else:
                    val = pd.Series([1000, 1000])
                    index2 = index[[0, -1]] 
                
                cols = file_rel_2D.columns
                file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__from_node']*len(index2), 
                                                              cols[1] : ['unit']*len(index2), 
                                                              cols[2] : ['node']*len(index2), 
                                                              cols[3] : [rels[r][0]]*len(index2), 
                                                              cols[4] : ['ElProduced']*len(index2), 
                                                              cols[5] : ['fuel_cost']*len(index2), 
                                                              cols[6] : [SC]*len(index2), 
                                                              cols[7] : index2, 
                                                              cols[8] : val, 
                                                              cols[9] : ['euro/MWh']*len(index2)}),
                                                  ignore_index=True)
            else:
                
                ### Only add value if in chosen regions
                if rels[r][0].split('_')[0] in reg_cho:
                    index2 = index
                else:
                    val = pd.Series([-1000, -1000])
                    index2 = index[[0, -1]] 
                    
                    
                cols = file_rel_2D.columns
                file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__from_node']*len(index2), 
                                                              cols[1] : ['unit']*len(index2), 
                                                              cols[2] : ['node']*len(index2), 
                                                              cols[3] : [rels[r][0]]*len(index2), 
                                                              cols[4] : [rels[r][1]]*len(index2), 
                                                              cols[5] : ['fuel_cost']*len(index2), 
                                                              cols[6] : [SC]*len(index2), 
                                                              cols[7] : index2, 
                                                              cols[8] : -val, 
                                                              cols[9] : ['euro/MWh']*len(index2)}),
                                                  ignore_index=True)
                

                
#%% --------------------- ###
###     Biomass Prices    ###
### --------------------- ###

### Not a time series! So doesn't work

# # COPY+PASTE HERE
# rels = pd.Series("""Z01_ST_HTL,Z01_Straw
# Z01_ST_Pyro,Z01_Straw
# Z01_ST_TG-FT,Z01_Straw
# Z01_ST_TG-FT-H2,Z01_Straw
# Z01_ST_TG-MeOH,Z01_Straw
# Z01_WO_HTL,Z01_Wood
# Z01_WO_Pyro,Z01_Wood
# Z01_WO_TG-FT,Z01_Wood
# Z01_WO_TG-FT-H2,Z01_Wood
# Z01_WO_TG-MeOH,Z01_Wood
# Z02_ST_HTL,Z02_Straw
# Z02_ST_Pyro,Z02_Straw
# Z02_ST_TG-FT,Z02_Straw
# Z02_ST_TG-FT-H2,Z02_Straw
# Z02_ST_TG-MeOH,Z02_Straw
# Z02_WO_HTL,Z02_Wood
# Z02_WO_Pyro,Z02_Wood
# Z02_WO_TG-FT,Z02_Wood
# Z02_WO_TG-FT-H2,Z02_Wood
# Z02_WO_TG-MeOH,Z02_Wood
# Z03_ST_HTL,Z03_Straw
# Z03_ST_Pyro,Z03_Straw
# Z03_ST_TG-FT,Z03_Straw
# Z03_ST_TG-FT-H2,Z03_Straw
# Z03_ST_TG-MeOH,Z03_Straw
# Z03_WO_HTL,Z03_Wood
# Z03_WO_Pyro,Z03_Wood
# Z03_WO_TG-FT,Z03_Wood
# Z03_WO_TG-FT-H2,Z03_Wood
# Z03_WO_TG-MeOH,Z03_Wood
# Z04_ST_HTL,Z04_Straw
# Z04_ST_Pyro,Z04_Straw
# Z04_ST_TG-FT,Z04_Straw
# Z04_ST_TG-FT-H2,Z04_Straw
# Z04_ST_TG-MeOH,Z04_Straw
# Z04_WO_HTL,Z04_Wood
# Z04_WO_Pyro,Z04_Wood
# Z04_WO_TG-FT,Z04_Wood
# Z04_WO_TG-FT-H2,Z04_Wood
# Z04_WO_TG-MeOH,Z04_Wood
# Z05_ST_HTL,Z05_Straw
# Z05_ST_Pyro,Z05_Straw
# Z05_ST_TG-FT,Z05_Straw
# Z05_ST_TG-FT-H2,Z05_Straw
# Z05_ST_TG-MeOH,Z05_Straw
# Z05_WO_HTL,Z05_Wood
# Z05_WO_Pyro,Z05_Wood
# Z05_WO_TG-FT,Z05_Wood
# Z05_WO_TG-FT-H2,Z05_Wood
# Z05_WO_TG-MeOH,Z05_Wood
# Z06_ST_HTL,Z06_Straw
# Z06_ST_Pyro,Z06_Straw
# Z06_ST_TG-FT,Z06_Straw
# Z06_ST_TG-FT-H2,Z06_Straw
# Z06_ST_TG-MeOH,Z06_Straw
# Z06_WO_HTL,Z06_Wood
# Z06_WO_Pyro,Z06_Wood
# Z06_WO_TG-FT,Z06_Wood
# Z06_WO_TG-FT-H2,Z06_Wood
# Z06_WO_TG-MeOH,Z06_Wood
# Z07_ST_HTL,Z07_Straw
# Z07_ST_Pyro,Z07_Straw
# Z07_ST_TG-FT,Z07_Straw
# Z07_ST_TG-FT-H2,Z07_Straw
# Z07_ST_TG-MeOH,Z07_Straw
# Z07_WO_HTL,Z07_Wood
# Z07_WO_Pyro,Z07_Wood
# Z07_WO_TG-FT,Z07_Wood
# Z07_WO_TG-FT-H2,Z07_Wood
# Z07_WO_TG-MeOH,Z07_Wood
# Z08_ST_HTL,Z08_Straw
# Z08_ST_Pyro,Z08_Straw
# Z08_ST_TG-FT,Z08_Straw
# Z08_ST_TG-FT-H2,Z08_Straw
# Z08_ST_TG-MeOH,Z08_Straw
# Z08_WO_HTL,Z08_Wood
# Z08_WO_Pyro,Z08_Wood
# Z08_WO_TG-FT,Z08_Wood
# Z08_WO_TG-FT-H2,Z08_Wood
# Z08_WO_TG-MeOH,Z08_Wood
# Z09_ST_HTL,Z09_Straw
# Z09_ST_Pyro,Z09_Straw
# Z09_ST_TG-FT,Z09_Straw
# Z09_ST_TG-FT-H2,Z09_Straw
# Z09_ST_TG-MeOH,Z09_Straw
# Z09_WO_HTL,Z09_Wood
# Z09_WO_Pyro,Z09_Wood
# Z09_WO_TG-FT,Z09_Wood
# Z09_WO_TG-FT-H2,Z09_Wood
# Z09_WO_TG-MeOH,Z09_Wood
# Z10_ST_HTL,Z10_Straw
# Z10_ST_Pyro,Z10_Straw
# Z10_ST_TG-FT,Z10_Straw
# Z10_ST_TG-FT-H2,Z10_Straw
# Z10_ST_TG-MeOH,Z10_Straw
# Z10_WO_HTL,Z10_Wood
# Z10_WO_Pyro,Z10_Wood
# Z10_WO_TG-FT,Z10_Wood
# Z10_WO_TG-FT-H2,Z10_Wood
# Z10_WO_TG-MeOH,Z10_Wood
# Z11_ST_HTL,Z11_Straw
# Z11_ST_Pyro,Z11_Straw
# Z11_ST_TG-FT,Z11_Straw
# Z11_ST_TG-FT-H2,Z11_Straw
# Z11_ST_TG-MeOH,Z11_Straw
# Z11_WO_HTL,Z11_Wood
# Z11_WO_Pyro,Z11_Wood
# Z11_WO_TG-FT,Z11_Wood
# Z11_WO_TG-FT-H2,Z11_Wood
# Z11_WO_TG-MeOH,Z11_Wood""".split('\n')).str.split(',')

# # Heat areas - HARDCODED FOR NOW - maybe make dictionary for ESB, CPH and so on?
# # AAA = """Z01_Large
# # Z02_Medium
# # Z03_A2""".split('\n')

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
rels = pd.Series("""Z01_SellingWasteheat,Z01_PotWasteheat
Z02_SellingWasteheat,Z02_PotWasteheat
Z03_SellingWasteheat,Z03_PotWasteheat
Z04_SellingWasteheat,Z04_PotWasteheat
Z05_SellingWasteheat,Z05_PotWasteheat
Z06_SellingWasteheat,Z06_PotWasteheat
Z07_SellingWasteheat,Z07_PotWasteheat
Z08_SellingWasteheat,Z08_PotWasteheat
Z09_SellingWasteheat,Z09_PotWasteheat
Z10_SellingWasteheat,Z10_PotWasteheat
Z11_SellingWasteheat,Z11_PotWasteheat""".split('\n')).str.split(',')

# Heat areas - HARDCODED FOR NOW - maybe make dictionary in inputs file for Z01 = DK1 and so on?
# AAA = pd.Series("""DK1_Large
# DK2_Large
# SE4_large
# NO2_A2
# NO1_A3
# SE3_large
# NO3_A3
# NO4_A2
# NO5_A2
# SE2_medium
# SE1_medium""".split('\n')).str.split('_')

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
            idxH2 = f_HPrice['RRR'] == reg_map[rels[r][0][:3]]
            a0 = aaa_map[rels[r][0][:3]].replace(' ', '').lower()
            idxH3 = f_HPrice['AAA'].str.lower() == a0
            
            vals = f_HPrice[idxH0 & idxH1 & idxH2 & idxH3]['Val'].iloc[::t_res].round(1)
            # vals = f_HPrice[idxH0 & idxH1 & idxH2 & idxH3]['Val'].iloc[::t_res].round(1) / (1 + 1/0.99) # For connection_flow_cost, where cost is price*(flow in + flow out)
            
            # Make sure it's not profitable to just produce heat and sell it
            if len(vals[vals <= 0]) >= 1:
                # print(vals[vals < 0])
                vals[vals <= 0] = 1
                # print(vals[vals < 0])
                # print(np.unique(vals))
            
            
            ### Only add value if in chosen regions
            if (rels[r][0].split('_')[0] in reg_cho) & (a0 != 'none'):
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
            
            ### If area not chosen or no heat area
            else:
                file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__from_node']*2, 
                                                              cols[1] : ['unit']*2, 
                                                              cols[2] : ['node']*2, 
                                                              cols[3] : [rels[r][0]]*2, 
                                                              cols[4] : [rels[r][1]]*2, 
                                                              cols[5] : ['fuel_cost']*2, 
                                                              cols[6] : [SC]*2, 
                                                              cols[7] : index[[0, -1]], 
                                                              cols[8] : np.zeros(2), 
                                                              cols[9] : ['euro/MWh']*2}))
                
                
#%% --------------------- ###
###      Heat Demand      ###
### --------------------- ###

### connection_capacity

# COPY+PASTE HERE
rels = pd.Series("""Z01_SellingWasteheat,Z01_DistrictHeating
Z02_SellingWasteheat,Z02_DistrictHeating
Z03_SellingWasteheat,Z03_DistrictHeating
Z04_SellingWasteheat,Z04_DistrictHeating
Z05_SellingWasteheat,Z05_DistrictHeating
Z06_SellingWasteheat,Z06_DistrictHeating
Z07_SellingWasteheat,Z07_DistrictHeating
Z08_SellingWasteheat,Z08_DistrictHeating
Z09_SellingWasteheat,Z09_DistrictHeating
Z10_SellingWasteheat,Z10_DistrictHeating
Z11_SellingWasteheat,Z11_DistrictHeating""".split('\n')).str.split(',')

# Heat areas - HARDCODED FOR NOW - maybe make dictionary for ESB, CPH and so on?
# AAA = pd.Series("""DK1_Large
# DK2_Large
# SE4_large
# NO2_A2
# NO1_A3
# SE3_large
# NO3_A3
# NO4_A2
# NO5_A2
# SE2_medium
# SE1_medium""".split('\n')).str.split('_')


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
            idxH2 = fnew['RRR'] == reg_map[rels[r][0][:3]]
            a0 = aaa_map[rels[r][0][:3]].replace(' ', '').lower()
            idxH3 = fnew['AAA'].str.lower() == a0

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


            if (rels[r][0].split('_')[0] in reg_cho) & (a0 != 'none'):
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
                
            else:
                file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__to_node']*2, 
                                                              cols[1] : ['unit']*2, 
                                                              cols[2] : ['node']*2, 
                                                              cols[3] : [rels[r][0]]*2, 
                                                              cols[4] : [rels[r][1]]*2, 
                                                              cols[5] : ['unit_capacity']*2, 
                                                              cols[6] : [SC]*2, 
                                                              cols[7] : index[[0, -1]], 
                                                              cols[8] : np.zeros(2), 
                                                              cols[9] : ['euro/MWh']*2}))
                


#%% --------------------- ###
### Heat Production Price ###
### --------------------- ###

# COPY+PASTE HERE
units = """Z01_Heat
Z02_Heat
Z03_Heat
Z04_Heat
Z05_Heat
Z06_Heat
Z07_Heat
Z08_Heat
Z09_Heat
Z10_Heat
Z11_Heat""".split('\n')

# Heat areas - HARDCODED FOR NOW - maybe make dictionary for ESB, CPH and so on?
# AAA = pd.Series("""DK1_Large
# DK2_Large
# SE4_large
# NO2_A2
# NO1_A3
# SE3_large
# NO3_A3
# NO4_A2
# NO5_A2
# SE2_medium
# SE1_medium""".split('\n')).str.split('_')

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
            idxH2 = f_HPrice['RRR'] == reg_map[units[r][:3]]
            idxH3 = f_HPrice['AAA'] == aaa_map[units[r][:3]]
                
            # Get prices and set negative to zero
            # vals = f_HPrice[idxH0 & idxH1 & idxH2 & idxH3]['Val'].iloc[::t_res].round(1) # Heat prices from Balmorel areas
            vals = np.array([25]*len(index)) ### Just 25 €/MWh in general
            # print()
            # print(index)
            # print(units[r])
            # print(units[r].replace('Heat', 'PotWasteheat'))
            if len(vals[vals <= 25]) >= 1:
                # print(vals[vals < 0])
                vals[vals <= 25] = 25 # IF HEAT IS PRODUCED, IT IS HIGH TEMPERATURE HEAT 25 €/MWh from Chmielarz2021
                
                # print(vals[vals < 0])
            
            if units[r].split('_')[0] in reg_cho:
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
            else:
                file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__from_node']*2, 
                                                              cols[1] : ['unit']*2, 
                                                              cols[2] : ['node']*2, 
                                                              cols[3] : [units[r]]*2, 
                                                              cols[4] : ['HeatProduced']*2, 
                                                              cols[5] : ['fuel_cost']*2, 
                                                              cols[6] : [SC]*2, 
                                                              cols[7] : index[[0, -1]], 
                                                              cols[8] : [999, 999], 
                                                              cols[9] : ['euro/MWh']*2}))                
            
            
            

#%% --------------------- ###
###     CO2 Potentials    ###
### --------------------- ###
f_fuel = pd.read_csv(sys.argv[5])

# AAA = pd.Series("""Z01_Large
# Z02_Large
# Z11_large
# Z04_A2
# Z03_A3
# Z10_large
# Z05_A3
# Z06_A2
# Z07_A2
# Z09_medium
# Z08_medium""".split('\n'))


idxF = (f_fuel['FFF'] == 'WOODCHIPS') | (f_fuel['FFF'] == 'MUNIWASTE') | (f_fuel['FFF'] == 'WOODWASTE') | (f_fuel['FFF'] == 'STRAW') | (f_fuel['FFF'] == 'BIOGAS') | (f_fuel['FFF'] == 'WOODPELLETS') | (f_fuel['FFF'] == 'WOOD')
idxY = f_fuel['Y'] == years[0]
# idxY = f_fuel['Y'] == f_fuel['Y']
idxA = f_fuel['AAA'] != f_fuel['AAA']
for r in reg_all:
    idxA = idxA | (f_fuel['AAA'] == aaa_map[r])

### Only take potentials from area (idxA)
fuels = f_fuel[idxF & idxY & idxA].groupby(['Y', 'RRR','UNITS'], as_index=False) # Same emission factor for all
fuels = fuels.aggregate(np.sum)

# COPY+PASTE HERE
nodes = """Z01_CO2_PS
Z02_CO2_PS
Z03_CO2_PS
Z04_CO2_PS
Z05_CO2_PS
Z06_CO2_PS
Z07_CO2_PS
Z08_CO2_PS
Z09_CO2_PS
Z10_CO2_PS
Z11_CO2_PS""".split('\n')

index = [mod_time[0]-dt.timedelta(hours=t_res), mod_time[0]]
SC = 'Base'
# Create data
for node in nodes:
    # Get yearly fuel consumption if possible
    try:
        val = fuels[fuels['RRR'] == reg_map[node[:3]]]['Val'].values[0]
    except IndexError:
        val = 0
    
    EM = 94.6 # Emission factor in kg/GJ
    val = val * 3600000 * EM / 1000 # TWh to GJ to ton CO2
    vals = [val, 'NaN'] 
    
    if not(node[:3] in reg_cho):
        val = 0
        index = [mod_time[0], mod_time[-1]+dt.timedelta(hours=t_res)]
        vals = [0, 0]
    
    # print(node)
    print('CO2 Point Source Potential %s = %d tCO2'%(node, val))
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
