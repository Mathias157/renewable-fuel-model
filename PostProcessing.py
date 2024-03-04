import sys
import os
import pandas as pd
import datetime as dt
import numpy as np

SC = 'Base_Export'
# p =
path_i = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\.spinetoolbox\items\investedobjects\output'
path = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\InvestmentInput'
path_inputs = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure'
# path_i = r'C:\Users\mberos\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Dokumenter\Codes\Electrolysis Infrastructure\.spinetoolbox\items\investedobjects\output'
# path = r'C:\Users\mberos\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Dokumenter\Codes\Electrolysis Infrastructure\InvestmentInput'
# path_inputs = r'C:\Users\mberos\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Dokumenter\Codes\Electrolysis Infrastructure'

#%% ----------------------------- ###
###          Investments          ###
### ----------------------------- ###


# Read files
if SC == '':
    f_inv = pd.read_excel(path_i + '\\unit_investments.xlsx')
else:
    f_inv = pd.read_excel(path_i + '\\' + SC + '\\unit_investments.xlsx')

# Filter away linear programming artifacts
idx = f_inv['Value'] < 1e-9
f_inv = f_inv[~idx]

# Scenario
try:
    SC2 = f_inv['Scenario'].iloc[0].split('__')[0]

    # Year
    yr = f_inv['Time'].iloc[0][:4]

    # Input file
    fi = pd.read_excel(path_inputs + '\\Inputs.xlsx', sheet_name='Inputs')

    # Lifetimes
    T = fi[(fi['category'] == 'Lifetime') & (int(yr) == fi['index'])]

    # Change param to spineinput
    f_inv['Object class'] = f_inv['Parameter'].str.replace('s_invested', '')
    f_inv['Object class'] = f_inv['Object class'].str.replace(
        'storage', 'node')

    f_inv['Parameter'] = f_inv['Parameter'].str.replace('units', 'fix_units')
    f_inv['Parameter'] = f_inv['Parameter'].str.replace(
        'connections', 'fix_connections')
    f_inv['Parameter'] = f_inv['Parameter'].str.replace(
        'storages', 'fix_storages')

    f_inv['Scenario short'] = f_inv['Scenario'].str.split('__', expand=True)[0] + '_INV'

    # Add for readability
    f_inv['Region'] = f_inv['Unit'].str.split('_', expand=True)[0]
    f_inv['Tech'] = f_inv['Unit'].str.replace('TRAN', 'TRA')
    f_inv['Tech'] = f_inv['Tech'].str[4:]
    f_inv['Lifetime'] = 50  # Pipeline lifetime

    for tech in np.unique(f_inv['Tech']):
        idx = f_inv['Tech'] == tech

        if not('Pipe' in tech):
            f_inv.loc[idx, 'Lifetime'] = T[T['parameter']
                                           == tech]['value'].values[0]

    # 2030 Investments decommisioning
    f_inv['New Value'] = """=IF(AND(J2=20,LEFT(D2,4)="2030"),0,E2)"""

    if not(os.path.exists(path + '\\' + SC)):
        os.makedirs(path + '\\' + SC)
        f_inv.to_csv(path + '\\' + SC + '\\' + yr + '_investments.csv', index=False)
    else:
        f_inv.to_csv(path + '\\' + SC + '\\' + yr + '_investments.csv', index=False)

except IndexError:
    print('No investments')


#%% ----------------------------- ###
###         Demand Flows          ###
### ----------------------------- ###
path = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\.spinetoolbox\items\demandflows\output'
# path = r'C:\Users\mberos\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\.spinetoolbox\items\demandflows\output'

# Read files
if SC != '':
    path = path + '\\' + SC 

f = pd.read_csv(path + '\DemandFlows.csv')    
# try:
    # 
# except FileNotFoundError:

# Resolution
t1 = dt.datetime.strptime(f['Time'][1], '%Y-%m-%dT%H:%M:%S')
t0 = dt.datetime.strptime(f['Time'][0], '%Y-%m-%dT%H:%M:%S')



f.loc[f['Direction'] == 'from_node', 'Value'] = - \
    f.loc[f['Direction'] == 'from_node', 'Value']
# yr = int(f['Time'][0][:4])
yr = f['Time'].str[:4]

f['Year'] = yr
fyear = f.groupby(['Unit', 'Node', 'Scenario', 'Year', 'Value'], as_index=False)
fyear = fyear.aggregate({'Value': np.sum})

unit2 = fyear['Unit'].copy()
fuel = fyear['Unit'].copy()

# Do operations
unit2 = unit2.str.split('_', expand=True)
node2 = fyear['Node'].copy()

# Jetfuel
tran_idx = unit2[0] == 'TRAN'
jet_idx = node2.str.find('Jetfuel') != -1
foss_idx = unit2[0] == 'Fossil'
fuel.loc[foss_idx] = 'Fossil'
fuel.loc[jet_idx & ~foss_idx & ~
         tran_idx] = fuel.loc[jet_idx & ~foss_idx].str[4:]
fuel.loc[jet_idx & ~foss_idx &
         tran_idx] = unit2.loc[jet_idx & ~foss_idx & tran_idx][2]

# Shipping and Road
tran_idx = (unit2[1] == 'Shipping') | (unit2[1] == 'Road')
# fuel.loc[~jet_idx & ~foss_idx & ~tran_idx] = unit2[1]
# fuel.loc[~jet_idx & ~foss_idx & tran_idx] = unit2[2]
fuel.loc[~jet_idx & ~foss_idx] = unit2[2]

# Region
reg = node2.str[:3]

# Excess
exc_idx = unit2[0] == 'EXCESS'
fuel.loc[exc_idx] = 'EXCESS'

# Type2

# Save
fyear['Value'] = fyear['Value']*168 # Assuming 7 day resolution
fyear['Reg'] = reg
fyear['Type'] = fuel
fyear['Type2'] = """=IF(ISNUMBER(FIND(F2,G2,1)), "TRAN", G2)"""
fyear['DemType'] = """=MID(B2,FIND("_", B2,1)+1,30)"""
fyear['SC'] = """=LEFT(C2,FIND("__",C2,1)-1)"""
# fyear['Year'] = yr
fyear['C'] = reg.str[:2]
# fyear.to_csv(path + '\\' + SC + '\DemandFlowsYearly_'+str(yr)+'.csv', index=False)
fyear.to_csv(path + '\DemandFlowsYearly'+'.csv', index=False)


#%% ----------------------------- ###
###          Other Flows          ###
### ----------------------------- ###
path = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\.spinetoolbox\items\fuelflows\output'
# path = r'C:\Users\mberos\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Dokumenter\Codes\Electrolysis Infrastructure\.spinetoolbox\items\fuelflows\output'

# Read files
if SC != '':
    path = path + '\\' + SC 

f = pd.read_csv(path + '\FuelFlows.csv')    
# try:
    # 
# except FileNotFoundError:

# Resolution
t1 = dt.datetime.strptime(f['Time'][1], '%Y-%m-%dT%H:%M:%S')
t0 = dt.datetime.strptime(f['Time'][0], '%Y-%m-%dT%H:%M:%S')


# yr = int(f['Time'][0][:4])
yr = f['Time'].str[:4]

f['Year'] = yr
# f['Value'] = f['Value'] * 4 # If 4h resolution
fyear = f.groupby(['Unit', 'Node', 'Direction', 'Scenario', 'Year', 'Value'], as_index=False)
fyear = fyear.aggregate({'Value': sum})

unit2 = fyear['Unit'].copy()

idx = unit2.str.find('Dem') != -1

# Do operations
# unit2 = unit2.str.split('_', expand=True)
# unit2[1][unit2[1] == 'ST'] = unit2[2][unit2[1] == 'ST']
# unit2[1][unit2[1] == 'WO'] = unit2[2][unit2[1] == 'WO']
# unit2[1][unit2[1] == 'Electrolyser'] = unit2[2][unit2[1] == 'Electrolyser']
# unit2[1][unit2[1] == 'CCU'] = unit2[1][unit2[1] == 'CCU'] + '-' + unit2[2][unit2[1] == 'CCU']
node2 = fyear['Node'].copy()

# Region
reg = node2.str[:3]


# # Save
fyear['Value'] = fyear['Value'] 
fyear['Reg'] = reg
# fyear['Type'] = node2.str[4:]
# fyear.loc[fyear['Type'] == '', 'Type'] = 'Electricity'
# fyear.loc[fyear['Type'] == 'oduced', 'Type'] = 'Electricity Sink'
# fyear.loc[fyear['Type'] == 'Produced', 'Type'] = 'Heat Sink'
# fyear.loc[fyear['Type'] == 'SS_Road', 'Type'] = 'Excess Road'
# fyear.loc[fyear['Type'] == 'atGas', 'Type'] = 'Natural Gas'
# fyear.loc[fyear['Type'] == '2', 'Type'] = 'H2'
# fyear['DemandType'] = ''
# fyear.loc[idx, 'DemandType'] = unit2[3][idx]
# fyear['Tech'] = unit2[1]
# fyear['Fuel'] = unit2[2]
# fyear.loc[(unit2[0] == 'TRAN') & (unit2[2].str.find('-') != -1), 'Tech'] = 'TRAN'
# fyear.loc[(unit2[0] == 'TRAN') & (unit2[2].str.find('-') == -1), 'Tech'] = unit2[2][(unit2[0] == 'TRAN') & (unit2[2].str.find('-') == -1)]
fyear['SC'] = fyear['Scenario'].str.split('__',expand=True)[0]
# fyear['Year'] = yr
fyear['C'] = reg.str[:2]
fyear['Tech'] = unit2.str[4:]
fyear['NodeType'] = node2.str[4:]
fyear.loc[fyear['NodeType'] == '', 'NodeType'] = 'Electricity'
# # fyear.to_csv(path + '\\' + SC + '\DemandFlowsYearly_'+str(yr)+'.csv', index=False)
fyear.to_csv(path + '\FuelFlowsYearly'+'.csv', index=False)



#%% ----------------------------- ###
###       Transport Flows         ###
### ----------------------------- ###

units = fyear['Unit'].copy()
units = units.str.split('_', expand=True)

nodes = fyear['Node'].copy()
nodes = nodes.str.split('_', expand=True)

idx = (units[0] == 'TRAN') & (units[3].str[4:] != 'Dem')

# print(fyear.loc[idx, ['Unit']])
f = fyear.loc[idx, ['Unit', 'Year', 'Direction', 'Reg', 'SC', 'Value']]
f['Type'] = nodes[idx][1]

tran_reg = nodes[idx][1].copy() # Placeholder
tran_reg[f['Type'] == 'H2'] = units[3][idx & (f['Type'] == 'H2')]
tran_reg[f['Type'] != 'H2'] = units[2][idx & (f['Type'] != 'H2')]

f['TranReg'] = tran_reg

### Coordinates
f['Coords'] = ''

coord = pd.read_excel('TransportGIS.xlsx')
for reg in tran_reg.unique():
    idx = f['TranReg'] == reg
    
    f['Coords'][idx] = coord['wkt_geom'][coord['Name'] == reg].values[0]

f = f.groupby(list(f.columns[:5])+list(f.columns[6:]), as_index=False)
f = f.aggregate({'Value': sum})

f.to_csv(path + '\TransportYearly.csv', index=False)
# f = pd.DataFrame(fyearf)


#%% ----------------------------- ###
###          Fuel Costs           ###
### ----------------------------- ###
# path = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\.spinetoolbox\items\heat&el\output'

# # ftime = pd.read_csv(path_inputs + '\\' + 'Timeseries.csv')
# # fhe = pd.read_excel(path + '\\' + SC + '\\heat&el.xlsx' )

# ### Heat and el
# idx1 = ftime['parameter names'] == 'connection_flow_cost'
# flows = np.unique(fhe['Unit'])

# fyear2 = pd.DataFrame({'Unit' : '',
#                        'Country' : '',
#                        'Region' : '',
#                        'Scenario' : '',
#                        'Type' : '',
#                        'Value' : '',
#                        'Unit' : ''},index=[0]).drop(0)
# cols = ['Unit', 'Country', 'Region', 'Scenario', 'Type', 'Value', 'Unit']

# for flow in flows[flows != 'Unit']:

#     # Check if heat or electricity
#     if flow[-2:] == 'El':
#         flow_cost = flow.replace('El', '')
#         cost_type = 'Electricity Costs'
#     elif flow[-4:] == 'Heat':
#         flow_cost = flow.replace('Heat', 'SellingWasteheat')
#         cost_type = 'Heat Costs'
#     else:
#         flow_cost = 'None'

#     idx2 = ftime['object name'] == flow_cost

#     # print(ftime[idx1 & idx2])

#     prices = np.array(fhe[fhe['Unit'] == flow]['Value'])
#     conn_flow = np.array(ftime[idx1 & idx2]['parameter values'])


#     if len(conn_flow) > 0:
#         SC = fhe[fhe['Unit'] == flow]['Scenario'].values[0]
#         obj = ftime[idx1 & idx2]['object name'].values[0]
#         fuel_cost_yr = np.sum(prices*conn_flow)
#         # fig = plt.figure()
#         # plt.plot(prices)
#         # plt.plot(conn_flow)
#         # plt.plot(prices*conn_flow)
#         # plt.legend(['Price','Flow', 'Total cost of %s'%obj])

#         fyear2 = fyear2.append(pd.DataFrame({cols[0] : [obj],
#                                              cols[1] : [obj[:2]],
#                                              cols[2] : [obj[:3]],
#                                              cols[3] : [SC],
#                                              cols[4] : cost_type,
#                                              cols[5] : fuel_cost_yr / 1e6,
#                                              cols[6] : ['Meuro']}),
#                                ignore_index=True)

