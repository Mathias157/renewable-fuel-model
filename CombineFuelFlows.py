# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 11:05:54 2022

@author: mathi
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from formplot import *
from scipy.optimize import curve_fit
import os

style = 'ppt'

if style == 'report':
    plt.style.use('default')
    fc = 'white'
elif style == 'ppt':
    plt.style.use('dark_background')
    fc = 'none'

#%% ----------------------------- ###
###             Auto              ###
### ----------------------------- ###

SC = 'HUB_FullBio'
p = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\.spinetoolbox\items\fuelflows\output'
if SC != '':
    p = p + "\\" + SC + "\\"
    

i = 0
for file in os.listdir(p):
    if (file != 'FuelFlows.csv') & (file != 'FuelFlowsYearly.csv') & (file != 'TransportYearly.csv'):
        if i == 0:
            f = pd.read_csv(p + file)
            print(f)
        else:
            f = f.append(pd.read_csv(p + file))
        i += 1
        
f.to_csv(p + 'FuelFlows.csv', index=False)


#%% ----------------------------- ###
###             Manual            ###
### ----------------------------- ###

i = 0
for SC in ['Base', 'Base_SecondRun']:
    p = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\.spinetoolbox\items\fuelflows\output'
    if SC != '':
        p = p + "\\" + SC + "\\"
      
    if i == 0:
        f = pd.read_csv(p + 'FuelFlows.csv')
    
    else:
        f = f.append(pd.read_csv(p + 'FuelFlows.csv'))
        
    i += 1        

## To last
p = r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\.spinetoolbox\items\fuelflows\output'

# f = f.append(pd.read_csv(p + '\FuelFlows_Sensitivities.csv'))
# f = f[f['Scenario'] != 'HUB_IT__SpineOpt@2022-06-01T15:09:16']
# f = f[f['Scenario'] != 'H2B_IT__SpineOpt@2022-06-01T15:09:16']
# f = f[f['Scenario'] != 'LBP_IT__SpineOpt@2022-06-01T16:12:13']
# f = f[f['Scenario'] != 'Base_IT__SpineOpt@2022-06-01T20:21:57']



f.to_csv(p + '\FuelFlows_Iterations.csv', index=False)