import sys
import pandas as pd
import datetime as dt
import numpy as np


#%% ----------------------------- ###
###       Functions & Read        ###
### ----------------------------- ###

sys.argv.append(r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\Inputs.xlsx') #f_inputs = pd.read_excel('Inputs.xlsx')

# Read files
f_inputs = pd.read_excel(sys.argv[1])

# Convert input array format - DOESN'T WORK FOR SOME REASON??? 
def inp_array(inp_array):
    """
    Converts input array format to float array

    Parameters
    ----------
    inp_array : String
        Example: '[2030, 2050]'.

    Returns
    -------
    out_array : Array
        Example: np.array([2030, 2050]).

    """
    temp = inp_array.rstrip(']').lstrip('[')
    out_array = pd.Series(temp.split(',')).astype(float)
    return out_array

# File for output
file = pd.DataFrame([], columns=['object class name', 'object name', 'parameter names', 'alternative names',	
                                 'parameter index', 'parameter values',	 
                                 'unit'])

file_rel_3D = pd.DataFrame([], columns=['relationship class name', 'object class name 1', 'object class name 2',
                                        'object class name 3', 'object name 1', 'object name 2', 'object name 3',
                                        'parameter names', 'alternative names', 'parameter index', 'parameter values',	 
                                        'unit'])

file_rel_2D = pd.DataFrame([], columns=['relationship class name', 'object class name 1', 'object class name 2',
                                        'object name 1', 'object name 2', 'parameter names', 'alternative names', 
                                        'parameter index', 'parameter values', 'unit'])


#%% ----------------------------- ###
###        Temporal Scope         ###
### ----------------------------- ###


# Input times
years = f_inputs[f_inputs['parameter'] == 'years']['value'].iloc[0]
years = np.array(years.rstrip(']').lstrip('[').split(',')).astype(float)
#years = inp_array(years) # DOESN'T WORK!! Can't load pandas module into func
weeks = f_inputs[f_inputs['parameter'] == 'seasons']['value'].iloc[0]
weeks = np.array(weeks.rstrip(']').lstrip('[').split(',')).astype(float)
#weeks = inp_array(weeks) # DOESN'T WORK!! Can't load pandas module into func
hours = f_inputs[f_inputs['parameter'] == 'hours']['value'].iloc[0]
SC = f_inputs[f_inputs['parameter'] == 'hours']['scenario'].iloc[0]

rat = 52*168 / (len(weeks)*hours) 


### ---- Model Parameters ---- ###
start_time = dt.datetime(int(years[0]),1,1)
end_time = start_time + dt.timedelta(hours=len(weeks)*168*len(years))

# Save
st_string = '{"data" : "%s"}::date_time'%dt.datetime.strftime(start_time, '%Y-%m-%dT%H:%M:%S')
file = file.append(pd.DataFrame([['model', 'simple', 'model_start', SC, 
                                  '', st_string,'date']], columns=file.columns),
                   ignore_index=True)

et_string = '{"data" : "%s"}::date_time'%dt.datetime.strftime(end_time, '%Y-%m-%dT%H:%M:%S')
file = file.append(pd.DataFrame([['model', 'simple', 'model_end', SC, 
                                  '', et_string,'date']], columns=file.columns),
                   ignore_index=True)



### ---- Temporal Block Parameters ---- ###
block_end = end_time - start_time
block_end = block_end.days*24 + block_end.seconds/3600

# Weeks block end and resolution
file = file.append(pd.DataFrame([['temporal_block', 'weeks', 'block_end', SC, 
                                  '', '{"data": "%sh"}::duration'%str(int(block_end)),'duration']], columns=file.columns),
                   ignore_index=True)

file = file.append(pd.DataFrame([['temporal_block', 'weeks', 'resolution', SC, 
                                  '', '{"data": "%sh"}::duration'%str(int(168/hours)),'duration']], columns=file.columns),
                   ignore_index=True)

# Investment block end and resolution
file = file.append(pd.DataFrame([['temporal_block', 'investment', 'resolution', SC, 
                                  '', '{"data": "%sh"}::duration'%str(int(block_end/len(years))),'duration']], columns=file.columns),
                   ignore_index=True)

#%% ----------------------------- ###
###      Geographical Scope       ###
### ----------------------------- ###

all_reg = f_inputs[f_inputs['parameter'] == 'all_regions']['value'].iloc[0]
all_reg = pd.Series(all_reg.rstrip(']').lstrip('[').split(',')).str.replace(' ','')
reg = f_inputs[f_inputs['parameter'] == 'regions']['value'].iloc[0]
reg = pd.Series(reg.rstrip(']').lstrip('[').split(',')).str.replace(' ','')

#%% --------------------- ###
###         CAPEX         ###
### --------------------- ###
Tidx = f_inputs['index'] == years[0]

CAPEX = f_inputs[(f_inputs['category'] == 'CAPEX') & Tidx]
FOM = f_inputs[(f_inputs['category'] == 'Fixed OM') & Tidx]
LifeT = f_inputs[(f_inputs['category'] == 'Lifetime') & Tidx]

reg = f_inputs[f_inputs['parameter'] == 'regions']['value'].iloc[0]
reg = pd.Series(reg.rstrip(']').lstrip('[').split(',')).str.replace(' ','')
reg[11] = 'NOS'
reg[12] = 'DKS'

# Index
# index = act_time[::ind].year
# index_mod = mod_time[::ind]

for r in range(len(reg)):
    # for n in range(len(CAPEX)):
    for n in range(len(CAPEX)):
        tech = reg[r] + '_' + CAPEX['parameter'].iloc[n]
        tech = tech.replace('__', '_') # Takes care if no index
        
        # True if not North Sea areas
        noNOS = (reg[r] != 'NOS') 
        noDKS = (reg[r] != 'DKS') 
        # True if offshore electrolyser
        OFFAEC = tech.find('Electrolyser_AECOFF') != -1 
        # OFFAECDK = tech.find('Electrolyser_AECOFFDK') != -1 
        if (noNOS & noDKS & ~OFFAEC) | (~(noNOS & noDKS) & OFFAEC):

            if (tech.split('_')[1] == 'SMR') & ((tech.split('_')[0][:2] != 'NO') | (tech.split('_')[0][:3] == 'NO1')):
                # Only allow SMR in Norwegian regions with gas infrastructure
                pass
            else:
                if ('STO' in tech) & ~('CON' in tech):
                    cl = 'node'
                    param = 'storage_investment_cost'
                elif ('STO' in tech) & ('CON' in tech):
                    cl = 'connection'
                    param = 'connection_investment_cost'
                elif 'TRAN' in tech:
                    cl = 'connection'
                    param = 'connection_investment_cost'
                else:
                    cl = 'unit'
                    param = 'unit_investment_cost'
        
                SC = CAPEX['scenario'].iloc[n]
        
                # data_time = CAPEX['index'].iloc[n] == years[0]
                val = CAPEX['value'].iloc[n]
                
                # Find fixed o&m and lifetime if exist
                idx1a = FOM['parameter'] == CAPEX['parameter'].iloc[n]
                idx1b = LifeT['parameter'] == CAPEX['parameter'].iloc[n]
                
                
                try:
                    idx2a = FOM['index'] == years[0]
                    fOM = FOM[idx1a & idx2a]['value'].values[0]
                    
                except IndexError: 
                    # If there's no value for the modeled year
                    fOM = 0
                    
                try:
                    idx2b = LifeT['index'] == years[0]
                    LifeTval = LifeT[idx1b & idx2b]['value'].values[0]
                    
                    
                    ### Annualisation
                    AnnCosts = np.zeros(LifeTval+1)
                    AnnCosts[0] = val
                    AnnCosts[1:] = fOM
                    
                    # Discount factor (disc rate of 4%)
                    disc_fac = np.array([1/(1+0.04)**i for i in range(LifeTval+1)])
                    
                    AnnCosts_disc = disc_fac*AnnCosts
                    
                    CRF = 0.04*(1+0.04)**LifeTval / ((1+0.04)**LifeTval - 1)
                    NPV = np.sum(AnnCosts_disc)
                    
                    # Annualised cost
                    ann_cost = NPV*CRF
                    # ann_cost = NPV*CRF / rat
                    ann_cost = np.int32(ann_cost)
                    
                    # Pipeline Annualised Costs
                    if reg[r] == 'DKS':
                        ann_cost = ann_cost + 10775.753
                    elif reg[r] == 'NOS':
                        ann_cost = ann_cost + 6887.605
                    
                    
                except IndexError:            
                    # If there's no value for the modeled year
                    LifeTval = 0
                    
                # Only add if 
                try:
                    # Check calc
                    print('%s Ann. cost = %0.0f'%(tech, ann_cost))
                    # print('%s NPV = %0.0f'%(tech, np.sum(AnnCosts_disc)))
                    # print('%s NPV via annualised cost = %0.0f'%(tech, np.sum(np.ones(LifeTval)*CRF*NPV*disc_fac[1:])))
                    
                    file = file.append(pd.DataFrame([[cl, tech, param, SC, 
                                                      '', ann_cost, 'EURO15/MW']], 
                                                    columns=file.columns), ignore_index=True)
                    
                except IndexError:
                    pass
            
                
      
    



#%% ----------------------------- ###
###         Efficiencies          ###
### ----------------------------- ###

### Assumes one year in simulation!

### fix_ratio_out_in
cols = file_rel_3D.columns

eta = f_inputs[f_inputs['category'] == 'efficiency']

idx = eta['index'] == years[0]
eta = eta[idx]

eta_names = eta['parameter'].str.rstrip(']').str.lstrip('[').str.split(', ')

reg = f_inputs[f_inputs['parameter'] == 'regions']['value'].iloc[0]
reg = pd.Series(reg.rstrip(']').lstrip('[').split(',')).str.replace(' ','')

for r in range(len(reg)):
    for n in range(len(eta)):
        tech = reg[r] + '_' + eta_names.iloc[n][0]
        tech = tech.replace('__', '_') # Takes care if no index
        
        # Nodetype
        ntype = {'Electricity' : reg[r],
                 'Straw' : reg[r] + '_Straw',
                 'Wood' : reg[r] + '_Wood',
                 'Gasoline' : reg[r] + '_Gasoline',
                 'Jetfuel' : reg[r] + '_Kerosene',
                 'NH3' : reg[r] + '_NH3',
                 'H2' : reg[r] + '_H2',
                 'DH' : reg[r] + '_PotWasteheat',
                 'CO2' : reg[r] + '_CO2',
                 'CO2_PS' : reg[r] + '_CO2_PS',
                 'CO2_STO' : reg[r] + '_CO2_STO',
                 'H2_STO_SMA' : reg[r] + '_H2_STO_SMA',
                 'Diesel' : reg[r] + '_Diesel',
                 'MeOH' : reg[r] + '_eMethanol',
                 'Bio-Oil' : reg[r] + '_Bio-Oil'}
        
        SC = eta['scenario'].iloc[n]
        
        node1 = ntype[eta_names.iloc[n][1]]
        node2 = ntype[eta_names.iloc[n][2]]
        val = np.round(eta['value'].iloc[n], 3)
        
        if ((tech.find('CCU_DAC') != -1) | (tech.find('SOEC') != -1) | (tech.find('CCU_PS') != -1)) & (node2.find('PotWasteheat') != -1):
            param = 'fix_ratio_in_out_unit_flow'
            node1 = ntype[eta_names.iloc[n][2]]
            node2 = ntype[eta_names.iloc[n][1]]
            val = round(1/val,3)
            N = 2
        elif ((tech.find('TG-FT') != -1) | (tech.find('TG-MeOH') != -1)) & ((node2.find('Straw') == -1) & (node2.find('Wood') == -1) & (node2.find('H2') == -1)):
            param = 'fix_ratio_out_out_unit_flow'
            N = 1
        else:
            param = 'fix_ratio_out_in_unit_flow'
            N = 1            

        SCs = [SC, 'ExclHighProcessHeat']        
        for i in range(N):
            if i == 1:
                node1 = node1.replace(reg[r]+'_PotWasteheat', 'HeatProduced')
                
            print(tech, node1, node2, SCs[i], val)
            file_rel_3D = file_rel_3D.append(pd.DataFrame({cols[0] : ['unit__node__node'], 
                                                          cols[1] : ['unit'], 
                                                          cols[2] : ['node'], 
                                                          cols[3] : ['node'],
                                                          cols[4] : [tech],
                                                          cols[5] : [node1],
                                                          cols[6] : [node2],
                                                          cols[7] : [param], 
                                                          cols[8] : [SCs[i]], 
                                                          cols[9] : '', 
                                                          cols[10] : [val], 
                                                          cols[11] : ['frac']}), ignore_index=True)

#%% --------------------- ###
###    Biomass Prices     ###
### --------------------- ###


# COPY+PASTE HERE
rels = pd.Series("""DK1_ST_HTL,DK1_Straw
DK1_ST_TG-FT,DK1_Straw
DK1_ST_TG-MeOH,DK1_Straw
DK1_WO_HTL,DK1_Wood
DK1_WO_TG-FT,DK1_Wood
DK1_WO_TG-MeOH,DK1_Wood
DK2_ST_HTL,DK2_Straw
DK2_ST_TG-FT,DK2_Straw
DK2_ST_TG-MeOH,DK2_Straw
DK2_WO_HTL,DK2_Wood
DK2_WO_TG-FT,DK2_Wood
DK2_WO_TG-MeOH,DK2_Wood
NO1_ST_HTL,NO1_Straw
NO1_ST_TG-FT,NO1_Straw
NO1_ST_TG-MeOH,NO1_Straw
NO1_WO_HTL,NO1_Wood
NO1_WO_TG-FT,NO1_Wood
NO1_WO_TG-MeOH,NO1_Wood
NO2_ST_HTL,NO2_Straw
NO2_ST_TG-FT,NO2_Straw
NO2_ST_TG-MeOH,NO2_Straw
NO2_WO_HTL,NO2_Wood
NO2_WO_TG-FT,NO2_Wood
NO2_WO_TG-MeOH,NO2_Wood
NO3_ST_HTL,NO3_Straw
NO3_ST_TG-FT,NO3_Straw
NO3_ST_TG-MeOH,NO3_Straw
NO3_WO_HTL,NO3_Wood
NO3_WO_TG-FT,NO3_Wood
NO3_WO_TG-MeOH,NO3_Wood
NO4_ST_HTL,NO4_Straw
NO4_ST_TG-FT,NO4_Straw
NO4_ST_TG-MeOH,NO4_Straw
NO4_WO_HTL,NO4_Wood
NO4_WO_TG-FT,NO4_Wood
NO4_WO_TG-MeOH,NO4_Wood
NO5_ST_HTL,NO5_Straw
NO5_ST_TG-FT,NO5_Straw
NO5_ST_TG-MeOH,NO5_Straw
NO5_WO_HTL,NO5_Wood
NO5_WO_TG-FT,NO5_Wood
NO5_WO_TG-MeOH,NO5_Wood
SE1_ST_HTL,SE1_Straw
SE1_ST_TG-FT,SE1_Straw
SE1_ST_TG-MeOH,SE1_Straw
SE1_WO_HTL,SE1_Wood
SE1_WO_TG-FT,SE1_Wood
SE1_WO_TG-MeOH,SE1_Wood
SE2_ST_HTL,SE2_Straw
SE2_ST_TG-FT,SE2_Straw
SE2_ST_TG-MeOH,SE2_Straw
SE2_WO_HTL,SE2_Wood
SE2_WO_TG-FT,SE2_Wood
SE2_WO_TG-MeOH,SE2_Wood
SE3_ST_HTL,SE3_Straw
SE3_ST_TG-FT,SE3_Straw
SE3_ST_TG-MeOH,SE3_Straw
SE3_WO_HTL,SE3_Wood
SE3_WO_TG-FT,SE3_Wood
SE3_WO_TG-MeOH,SE3_Wood
SE4_ST_HTL,SE4_Straw
SE4_ST_TG-FT,SE4_Straw
SE4_ST_TG-MeOH,SE4_Straw
SE4_WO_HTL,SE4_Wood
SE4_WO_TG-FT,SE4_Wood
SE4_WO_TG-MeOH,SE4_Wood
DK1_ST_TG-FT-H2,DK1_Straw
DK1_WO_TG-FT-H2,DK1_Wood
DK2_ST_TG-FT-H2,DK2_Straw
DK2_WO_TG-FT-H2,DK2_Wood
NO1_ST_TG-FT-H2,NO1_Straw
NO1_WO_TG-FT-H2,NO1_Wood
NO2_ST_TG-FT-H2,NO2_Straw
NO2_WO_TG-FT-H2,NO2_Wood
NO3_ST_TG-FT-H2,NO3_Straw
NO3_WO_TG-FT-H2,NO3_Wood
NO4_ST_TG-FT-H2,NO4_Straw
NO5_WO_TG-FT-H2,NO5_Wood
NO5_ST_TG-FT-H2,NO5_Straw
SE1_ST_TG-FT-H2,SE1_Straw
SE1_WO_TG-FT-H2,SE1_Wood
SE2_ST_TG-FT-H2,SE2_Straw
SE2_WO_TG-FT-H2,SE2_Wood
SE3_ST_TG-FT-H2,SE3_Straw
SE3_WO_TG-FT-H2,SE3_Wood
SE4_ST_TG-FT-H2,SE4_Straw
SE4_WO_TG-FT-H2,SE4_Wood""".split('\n')).str.split(',', expand=True)

# Heat areas - HARDCODED FOR NOW - maybe make dictionary for ESB, CPH and so on?
# AAA = """DK1_Large
# DK2_Medium
# NO1_A2""".split('\n')

# path = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure' + '/Balmorel Results/csv/'

# f_elPrice = pd.read_csv(path+SC+'/PriceElectricityHourly_'+SC+'.csv')
# f_bioPrice = pd.read_csv(sys.argv[2])
f_bioPrice = f_inputs[f_inputs['category'] == 'biomass price']

Tidx = years[0] == f_bioPrice['index']

# Make biomass prices
cols = file_rel_2D.columns
for SC in np.unique(f_bioPrice['scenario']):
    SCidx = f_bioPrice['scenario'] == SC
    # print(f_bioPrice[Tidx & SCidx])
    for r in range(len(rels)):
        Bidx = f_bioPrice['parameter'] == (rels.loc[r, 1].split('_')[1])
        
        file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__from_node'], 
                                                      cols[1] : ['unit'], 
                                                      cols[2] : ['node'], 
                                                      cols[3] : [rels.loc[r, 0]], 
                                                      cols[4] : [rels.loc[r, 1]], 
                                                      cols[5] : ['fuel_cost'], 
                                                      cols[6] : [f_bioPrice[Tidx & SCidx & Bidx]['scenario'].values[0]], 
                                                      cols[7] : '', 
                                                      cols[8] : [np.round(f_bioPrice[Tidx & SCidx & Bidx]['value'].values[0], 2)], 
                                                      cols[9] : ['euro/MWh']}))
        
    # End with bioimport price
    Bidx = f_bioPrice['parameter'] == 'Wood Import'
    wood_rels = rels[1].str.find('Wood') != -1
    wood_nodes = np.unique(rels[wood_rels].iloc[:,1])
    
    for r in range(len(wood_nodes)):
        file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__to_node'], 
                                                      cols[1] : ['unit'], 
                                                      cols[2] : ['node'], 
                                                      cols[3] : ['WoodImport'], 
                                                      cols[4] : [wood_nodes[r]], 
                                                      cols[5] : ['fuel_cost'], 
                                                      cols[6] : [f_bioPrice[Tidx & SCidx & Bidx]['scenario'].values[0]], 
                                                      cols[7] : '', 
                                                      cols[8] : [np.round(f_bioPrice[Tidx & SCidx & Bidx]['value'].values[0], 2)], 
                                                      cols[9] : ['euro/MWh']}))


#%% ----------------------------- ###
###         Variable OM           ###
### ----------------------------- ###

### Assumes one year in simulation!

### vom_cost
cols = file_rel_2D.columns

## Variable OM
VOM = f_inputs[f_inputs['category'] == 'VOM']

idx = VOM['index'] == years[0]
VOM = VOM[idx]

VOM_names = VOM['parameter'].str.rstrip(']').str.lstrip('[').str.split(', ', expand=True)

reg = f_inputs[f_inputs['parameter'] == 'regions']['value'].iloc[0]
reg = pd.Series(reg.rstrip(']').lstrip('[').split(',')).str.replace(' ','')

## Distribution costs - HARDCODED FOR NOW
fuel_dist = {'Bio-Oil' : 37.91722785/7.47*3.6,
             'MeOH' : 72.52489468/7.47*3.6,
             'Diesel' : 28.74482946/7.47*3.6,
             'Gasoline' : 34.83998799/7.47*3.6,
             'Jetfuel' : 2.166299727/7.47*3.6,
             'H2' : 0,
             'CO2' : 0,
             'NH3' : 0}

## Distribution costs - HARDCODED FOR NOW
fuel_dist = {'Bio-Oil' : 15,
             'MeOH' : 15,
             'Diesel' : 15,
             'Gasoline' : 15,
             'Jetfuel' : 15,
             'H2' : 0,
             'CO2' : 0,
             'NH3' : 15}

for r in range(len(reg)):
    for n in range(len(VOM)):
        tech = reg[r] + '_' + VOM_names.iloc[n][0]
        tech = tech.replace('__', '_') # Takes care if no index
        
        # Nodetype
        ntype = {'Electricity' : reg[r],
                 'Straw' : reg[r] + '_Straw',
                 'Wood' : reg[r] + '_Wood',
                 'Gasoline' : reg[r] + '_Gasoline',
                 'Jetfuel' : reg[r] + '_Kerosene',
                 'NH3' : reg[r] + '_NH3',
                 'H2' : reg[r] + '_H2',
                 'DH' : reg[r] + '_PotWasteheat',
                 'CO2' : reg[r] + '_CO2',
                 'CO2_PS' : reg[r] + '_CO2_PS',
                 'CO2_STO' : reg[r] + '_CO2_STO',
                 'Diesel' : reg[r] + '_Diesel',
                 'MeOH' : reg[r] + '_eMethanol',
                 'Bio-Oil' : reg[r] + '_Bio-Oil'}
        
        # Scenario
        SC = VOM['scenario'].iloc[n]
        
        # Node
        product = VOM_names.iloc[n][1]
        node1 = ntype[product]
        
        # VOM + Distribution Cost
        val = VOM['value'].iloc[n]
        
        try:
            # print(product, val)
            if tech.find('MeOH-Upgrade') == -1: # Don't double count distribution
                val = val + fuel_dist[product]
            # val = val + 15 # 15
            # print(val)
        except KeyError:
            # If no distribution cost
            pass
        
        # Round
        val = np.round(val, 2)
        
        if tech.find('Electrolyser_AECOFF') == -1:
            print(tech)
            print(val)
            # print(tech, node1, node2, SC, VOM['value'].iloc[n])
            
            file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__to_node'], 
                                                          cols[1] : ['unit'], 
                                                          cols[2] : ['node'], 
                                                          cols[3] : [tech],
                                                          cols[4] : [node1],
                                                          cols[5] : ['vom_cost'], 
                                                          cols[6] : [SC], 
                                                          cols[7] : '', 
                                                          cols[8] : [val], 
                                                          cols[9] : ['EURO15/MWh']}), ignore_index=True)


#%% ----------------------------- ###
###         Fossil Fuels          ###
### ----------------------------- ###
fossils = f_inputs[f_inputs['category'] == 'Fossil fuel price']
emi = f_inputs[f_inputs['category'] == 'Fossil fuel emissions']

reg = f_inputs[f_inputs['parameter'] == 'regions']['value'].iloc[0]
reg = pd.Series(reg.rstrip(']').lstrip('[').split(',')).str.replace(' ','')

for SC in np.unique(fossils['scenario']):
    Tidx = fossils['index'] == years[0]
    SCidx = SC == fossils['scenario']
    CO2tax = f_inputs[(f_inputs['category'] == 'CO2 Tax') & (f_inputs['index'] == years[0]) & (f_inputs['scenario'] == SC)]['value'].values[0] / 1000 # €15 / kg

    for fuel in np.unique(fossils['parameter']):
        
        fossil = fossils[(fossils['parameter'] == fuel) & SCidx & Tidx]['value'].values[0]
        em = emi[emi['parameter'] == fuel]['value'].values[0]
        
        # Fossil fuel price:
        FP = fossil + CO2tax*em
        
        # Dictionary
        dic_r = {'Jetfuel' : '_Jetfuel',
               'Diesel' : '_Road',
               'Heavy fuel oil' : '_Shipping'}
        
        dic_u = {'Jetfuel' : 'Fossil_Jetfuel',
               'Diesel' : 'Fossil_Diesel',
               'Heavy fuel oil' : 'Fossil_FuelOil'}
        
        for r in range(len(reg)):
            
            print(reg[r] + dic_r[fuel])
            print(FP)
            print(dic_u[fuel])
            FP = np.round(FP, 2)
            
            file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__to_node'], 
                                                          cols[1] : ['unit'], 
                                                          cols[2] : ['node'], 
                                                          cols[3] : [dic_u[fuel]],
                                                          cols[4] : [reg[r] + dic_r[fuel]],
                                                          cols[5] : ['vom_cost'], 
                                                          cols[6] : [SC], 
                                                          cols[7] : '', 
                                                          cols[8] : [FP], 
                                                          cols[9] : ['EURO15/MWh']}), ignore_index=True)

    
    # emi = fi[cat == 'Fossil fuel emissions']
    
    # emi[emi['parameter']=='Jetfuel']['value'].values[0]

#%% ----------------------------- ###
###         Fossil Share          ###
### ----------------------------- ###
dems = f_inputs[f_inputs['category'] == 'demand']
fos_share = f_inputs[f_inputs['category'] == 'fossil share']

Tidx = fos_share['index'] == years[0]
Tidx2 = dems['index'] == years[0]

dict_tech = {'Fossil_Jetfuel' : '_Jetfuel',
             'Fossil_Diesel' : '_Road',
             'Fossil_FuelOil' : '_Shipping'}

for SC in np.unique(fos_share['scenario']):
    SCidx = dems['scenario'] == SC
    for tech in np.unique(fos_share['parameter']):
        techidx = tech == fos_share['parameter']
        
        fos = fos_share[Tidx & techidx]
        
        for r in reg:
            Ridx = dems['parameter'].str.find(r) != -1
            techidx2 = dems['parameter'] == r + dict_tech[tech]
            
            # Retrieve hourly demand
            dem = dems[Tidx2 & SCidx & Ridx & techidx2]['value'].values[0]
            
            # Retrieve fossil fraction
            fs = fos_share[techidx & Tidx]['value'].values[0]
            
            if np.isnan(fs):
                val = 999999999999999999
            else:
                val = np.round(fs*dem*13, 2) # 4 week scale
            
            # Save
            # file = file.append(pd.DataFrame([['unit', tech, 
            #                                   'fix_unit_flow', SC, 
            #                                   '', val, 'MWh/h']], 
            #                                 columns=file.columns), ignore_index=True)
    
    
            file_rel_2D = file_rel_2D.append(pd.DataFrame({cols[0] : ['unit__to_node'], 
                                                          cols[1] : ['unit'], 
                                                          cols[2] : ['node'], 
                                                          cols[3] : [tech],
                                                          cols[4] : [r + dict_tech[tech]],
                                                          cols[5] : ['fix_unit_flow'], 
                                                          cols[6] : [SC], 
                                                          cols[7] : '', 
                                                          cols[8] : [val], 
                                                          cols[9] : ['MWh/h']}), ignore_index=True)


#%% ----------------------------- ###
###              Save             ###
### ----------------------------- ###

file.to_csv('scalars.csv', index=False)

file_rel_2D.to_csv('scalars_rel_2D.csv', index=False)

file_rel_3D.to_csv('scalars_rel_3D.csv', index=False)