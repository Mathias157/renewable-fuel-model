import sys
import os
import pandas as pd
import datetime as dt
import numpy as np

SC = 'H2B_Export'
# p =
path_i = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\.spinetoolbox\items\investedobjects\output'
path = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\InvestmentInput'
path_inputs = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure'
# path_i = r'C:\Users\mberos\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Dokumenter\Codes\Electrolysis Infrastructure\.spinetoolbox\items\investedobjects\output'
# path = r'C:\Users\mberos\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Dokumenter\Codes\Electrolysis Infrastructure\InvestmentInput'
# path_inputs = r'C:\Users\mberos\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Dokumenter\Codes\Electrolysis Infrastructure'

if SC != '':
    pp = '\\' + SC
else:
    pp = ''
#%% ----------------------------- ###
###           Heat & El           ###
### ----------------------------- ###

# COPY+PASTE HERE
conns = """DK1_PotWasteheat
DK2_PotWasteheat
NO1_PotWasteheat
NO2_PotWasteheat
NO3_PotWasteheat
NO4_PotWasteheat
NO5_PotWasteheat
SE1_PotWasteheat
SE2_PotWasteheat
SE3_PotWasteheat
SE4_PotWasteheat""".split('\n')

# Heat areas - HARDCODED FOR NOW - maybe make dictionary for ESB, CPH and so on?
AAA = pd.Series("""DK1_Large
DK2_Large
NO1_A3
NO2_A2
NO3_A3
NO4_A2
NO5_A2
SE1_medium
SE2_medium
SE3_large
SE4_large""".split('\n'))

dict_A = {conns[i]: AAA[i] for i in range(len(AAA))}

p = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\.spinetoolbox\items\heat&el\output'
# p = r'C:\Users\mberos\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Dokumenter\Codes\Electrolysis Infrastructure\.spinetoolbox\items\heat&el\output'

# f = pd.read_excel(p + '\\' + SC + '\\heat&el.xlsx')
f1 = pd.read_excel(p + pp + '\\heat&el.xlsx').dropna()
fi = pd.read_excel(path_inputs + '\\Inputs.xlsx', sheet_name='Inputs')

SC_array = f1['Scenario'].str.split('__').str[0]
SCs = np.unique(SC_array)
SCs = SCs[SCs != 'Scenario']
for SC in SCs:
    yrs = f1['Time'].str[:4]
    yrs = yrs[yrs != 'Time']
    for yr in np.unique(yrs):
        
        f = f1[(SC_array == SC) & (yrs == yr)]
        
        if len(f) > 0:
            # Get seasons
            S = fi[(fi['parameter'] == 'seasons') & (
                fi['category'] == 'temp_scope')]['value'].values[0]
            S = S.rstrip(']').lstrip('[')
            S = S.replace(' ', '')
            S = S.split(',')
            
            # Adapt time
            times = np.unique(f['Time'][f['Time'] != 'Time'])
            week_length = int(len(times) / len(S))
            
            
            ### Create GAMS file for heat
            fDH = 'DH commands'
            
            fh = pd.DataFrame([], columns=['S', 'T'], index=[])
            cols = fh.columns
            
            hours = 'T' + pd.Series(np.arange(1, 168+1, (168/week_length)).astype(int)).astype(str)
            hours[hours.str.len() == 2] = hours[hours.str.len()
                                                == 2].str.replace('T', 'T00')
            hours[hours.str.len() == 3] = hours[hours.str.len()
                                                == 3].str.replace('T', 'T0')
            
            
            # hours = hours[::int(168/week_length)]
            for i in range(len(S)):
                w = S[i]
                if len(w) == 1:
                    w = 'S0' + w
                else:
                    w = 'S' + w
                    
                fh = fh.append(pd.DataFrame({cols[0]: ['%s . ' % w]*week_length,
                                             cols[1]: hours}),
                               ignore_index=True)
            
            node_idx = (f['Node'].str.find('PotWasteheat') != -1)
            for node in np.unique(f['Node'][node_idx]):
                # Specific node index
                idx1 = node == f['Node']
            
                # Wasteheat sell
                idx3 = f['Unit'].str.find('SellingWasteheat') != -1
            
                # Heat produced
                idx4 = (f['Unit'].str.find('Heat') != -1)
            
                # Make sum: Positive demand, negative sold
                val1 = np.array(f['Value'][idx1 & idx3])
                val2 = np.array(f['Value'][idx1 & idx4])
                # Make sum, negative convention
                # val = np.sum(val2 - val1) # HIGH HEAT OUT OF SCOPE
                val = -np.sum(val1)*168/week_length # Use this for no high heat
                
                ## HIGH TEMPERATURE BYPASS
                # IF HEAT IS REQUIRED IT IS HIGH TEMPERATURE PROCESS HEAT!
                # Therefore, 25 €/MWh is used for heat production in general
                # See Chmielarz2021 and average of IND-HT area prices in Balmorel
                if val > 0:
                    val = 0 
                if val != 0:
                    sign = val / abs(val)
                else:
                    sign = 1
                string = "GKFX('%s','%s','GNR_EH_HEAT') = GKFX('%s','%s','GNR_EH_HEAT') + %d ;" % (yr, dict_A[node], yr, dict_A[node], -val/8736) # CHANGED TO GKFX, check that it's negative!
                fDH = fDH + '\n' + string
                
                
                # for i in range(len(S)):
                #     date_start = times[week_length*i]
                #     date_stop = times[week_length*i + week_length-1]
            
                #     # Time index
                #     idx2 = (f['Time'] >= date_start) & (f['Time'] <= date_stop)
            
            
                #     # Create data - variation should be positive, sum decides if it's demand or production
                #     if sign > 0:
                #         # If more electricity is used than sold
                #         val1 = -np.array(f['Value'][idx1 & idx2 & idx3]) # Negative sell convention
                #         val2 = np.array(f['Value'][idx1 & idx2 & idx4])  # Positive demand
                #     else:
                #         # If more electricity is sold than used
                #         val1 = np.array(f['Value'][idx1 & idx2 & idx3])  # Positive sell convention
                #         val2 = -np.array(f['Value'][idx1 & idx2 & idx4]) # Negative demand
                #     # val = val1 + val2
                #     val = val1*168/week_length # NO TAKING INTO ACCOUNT USED HEAT
            
                #     fh.loc[len(hours)*i:(i+1)*len(hours)-1, dict_A[node]] = val
            
            if not(os.path.exists(p + '\\' + SC)):
                os.makedirs(p + '\\' + SC)
            
            # Save
            DH_SAVE = open(p + '\\' + SC + '\\' + str(yr) + '_GKFX.txt', 'w')
            DH_SAVE.write(fDH)
            DH_SAVE.close()
            
            # fh.to_csv(p + '\\' + SC + '\\' + str(yr) + '_DH_VAR_T.csv', index=False)
            
            
            # %% Create GAMS file for electricity
            fDE = 'DE commands'
            
            fe = pd.DataFrame([], columns=['S', 'T'], index=[])
            cols = fh.columns
            
            hours = 'T' + pd.Series(np.arange(1, 169, (168/week_length)).astype(int)).astype(str)
            hours[hours.str.len() == 2] = hours[hours.str.len()
                                                == 2].str.replace('T', 'T00')
            hours[hours.str.len() == 3] = hours[hours.str.len()
                                                == 3].str.replace('T', 'T0')
            # hours = hours[::int(168/week_length)]
            for i in range(len(S)):
                w = S[i]
                if len(w) == 1:
                    w = 'S0' + w
                else:
                    w = 'S' + w
                fe = fe.append(pd.DataFrame({cols[0]: ['%s . ' % w]*week_length,
                                             cols[1]: hours}),
                               ignore_index=True)
            
            
            node_idx = (f['Node'].str.find('PotWasteheat') == -1) & (f['Node'] != 'Node')
            for node in np.unique(f['Node'][node_idx]):
                # Specific node index
                idx1 = node == f['Node']
            
                # Selling Electricity
                idx3 = f['Unit'].str.find('SellingElectricity') != -1
            
                # Electricity Used
                idx4 = (f['Unit'].str.find('_El') != -1) 
                
                val1 = np.array(f['Value'][idx1 & idx3])
                val2 = np.array(f['Value'][idx1 & idx4])
                # Make sum, positive = demand, negative = production
                val = np.sum(val2 - val1)*168/week_length
                if val != 0:
                    sign = val / abs(val)
                else:
                    sign = 1
                string = "DE('%s','%s','PtX') = %d ;"%(yr, node, val)
                fDE = fDE + '\n' + string
                
                for i in range(len(S)):
                    date_start = times[week_length*i]
                    date_stop = times[week_length*i + week_length-1]
            
                    # Time index
                    idx2 = (f['Time'] >= date_start) & (f['Time'] <= date_stop)
            
                    
                    # Create data - variation should be positive, sum decides if it's demand or production
                    if sign > 0:
                        # If more electricity is used than sold
                        val1 = -np.array(f['Value'][idx1 & idx2 & idx3]) # Negative sell convention
                        val2 = np.array(f['Value'][idx1 & idx2 & idx4])  # Positive demand
                    else:
                        # If more electricity is sold than used
                        val1 = np.array(f['Value'][idx1 & idx2 & idx3])  # Positive sell convention
                        val2 = -np.array(f['Value'][idx1 & idx2 & idx4]) # Negative demand
                        
                        
                    val = (val1 + val2)*168/week_length
                    fe.loc[len(hours)*i:(i+1)*len(hours)-1, node] = val # Must be positive
                    # print(fe)
                
            
            
            # Save
            DE_SAVE = open(p + '\\' + SC + '\\' + str(yr) + '_DE.txt', 'w')
            DE_SAVE.write(fDE)
            DE_SAVE.close()
            
            fe.to_csv(p + '\\' + SC + '\\' + str(yr) + '_DE_VAR_T.csv', index=False)


#%% HOW TO SAVE IT DIRECTLY AS INC:
f = open('directsave.inc', 'w')
with open('directsave.inc', 'a') as f:
    dfAsString = fe.to_string(header=True, index=False)
    f.write(dfAsString)