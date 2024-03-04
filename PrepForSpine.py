# -*- coding: utf-8 -*-
"""
Created on Sun Jan 23 11:29:21 2022

@author: mathi
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from formplot import *
from scipy.optimize import curve_fit
import sqlite3

style = 'ppt'

if style == 'report':
    plt.style.use('default')
    fc = 'white'
elif style == 'ppt':
    plt.style.use('dark_background')
    fc = 'none'

#%% ----------- Prepare Heat Price ----------- ###

f = pd.read_excel('EI_Data.xlsx', sheet_name='TimeseriesPrices', 
                  header=2)

T = pd.Series(f['Time'])
#T = T.str.replace(' ', 'T')

vals = pd.Series(f['HeatSellPrice (€20/MWh)'])

new_file = pd.DataFrame({'Time' : T, 'HeatPrice (€20/MWh)' : vals})
new_file.to_csv('HeatPrice.csv', index=False)
# 2030-01-01T00:00:00

### ----------- Python SQLITE Manipulation ----------- ###

# con = sqlite3.connect(r'C:\Users\mathi\Danmarks Tekniske Universitet\Master Thesis at Sustainability Division - Documents\Codes\Electrolysis Infrastructure\.spinetoolbox\items\datainput\DataInput.SQLITE')

# cur = con.cursor()

# sql_query = "SELECT name FROM sqlite_master WHERE type='table';"
# cur.execute(sql_query)

# ### View all tables
# print('\nAll tables in database')
# print(cur.fetchall())

# ### View table TABLE
# TABLE = 'parameter_definition'
# cur.execute('SELECT * FROM %s'%TABLE)
# # View columns
# print('\nAll columns in table "%s"'%TABLE)
# print(cur.description)

# # Print rows
# sort = 'description'
# print('\nRows in %s sorted by %s'%(TABLE, sort))
# for row in cur.execute('SELECT * FROM %s ORDER BY %s'%(TABLE, sort)):
#         print(row)
# # {"data": {"2030-01-01T00:00:00": -14.45217777007603, "2030-01-01T01:00:00": -14.45217777007603, "2030-01-01T02:00:00": -14.45217777007603, "2030-01-01T03:00:00": -14.45217777007603, "2030-01-01T04:00:00": -22.22682657111151, "2030-01-01T05:00:00": -22.680461169433208, "2030-01-01T06:00:00": -22.680461169433208, "2030-01-01T07:00:00": -22.680461169433208, "2030-01-01T08:00:00": -22.22682657111151, "2030-01-01T09:00:00": -22.22682657111151, "2030-01-01T10:00:00": -22.396607941243605, "2030-01-01T11:00:00": -22.680461169433208, "2030-01-01T12:00:00": -22.680461169433208, "2030-01-01T13:00:00": -23.134095767754907, "2030-01-01T14:00:00": -21.381495279197004, "2030-01-01T15:00:00": -19.23650605341804, "2030-01-01T16:00:00": -22.22682657111151, "2030-01-01T17:00:00": -22.22682657111151, "2030-01-01T18:00:00": -22.22682657111151, "2030-01-01T19:00:00": -22.22682657111151, "2030-01-01T20:00:00": -22.22682657111151, "2030-01-01T21:00:00": -20.617825135416687, "2030-01-01T22:00:00": -21.095681260727343, "2030-01-01T23:00:00": -21.095681260727343, "2030-01-02T00:00:00": -21.421057106612295, "2030-01-02T01:00:00": -21.512983568565613, "2030-01-02T02:00:00": -21.421057106612295, "2030-01-02T03:00:00": -22.22682657111151, "2030-01-02T04:00:00": -22.22682657111151, "2030-01-02T05:00:00": -22.680461169433208, "2030-01-02T06:00:00": -22.680461169433208, "2030-01-02T07:00:00": -22.680461169433208, "2030-01-02T08:00:00": -22.680461169433208, "2030-01-02T09:00:00": -22.22682657111151, "2030-01-02T10:00:00": -22.67204130528943, "2030-01-02T11:00:00": -29.709663642232435, "2030-01-02T12:00:00": -27.47736157868853, "2030-01-02T13:00:00": -22.680461169433208, "2030-01-02T14:00:00": -21.750585214266785, "2030-01-02T15:00:00": -10.804911770139771, "2030-01-02T16:00:00": -22.22682657111151, "2030-01-02T17:00:00": -22.22682657111151, "2030-01-02T18:00:00": -22.680461169433208, "2030-01-02T19:00:00": -22.680461169433208, "2030-01-02T20:00:00": -22.22682657111151, "2030-01-02T21:00:00": -22.22682657111151, "2030-01-02T22:00:00": -20.752312280507738, "2030-01-02T23:00:00": -21.509177328610203, "2030-01-03T00:00:00": -21.82255775160539, "2030-01-03T01:00:00": -21.82255775160539, "2030-01-03T02:00:00": -21.82255775160539, "2030-01-03T03:00:00": -21.82255775160539, "2030-01-03T04:00:00": -21.82255775160539, "2030-01-03T05:00:00": -22.22682657111151, "2030-01-03T06:00:00": -22.22682657111151, "2030-01-03T07:00:00": -22.22682657111151, "2030-01-03T08:00:00": -22.22682657111151, "2030-01-03T09:00:00": -23.134095767754907, "2030-01-03T10:00:00": -22.680461169433208, "2030-01-03T11:00:00": -22.680461169433208, "2030-01-03T12:00:00": -22.680461169433208, "2030-01-03T13:00:00": -22.680461169433208, "2030-01-03T14:00:00": -22.14470406055848, "2030-01-03T15:00:00": -19.305479735034204, "2030-01-03T16:00:00": -19.382757940189443, "2030-01-03T17:00:00": -22.22682657111151, "2030-01-03T18:00:00": -22.22682657111151, "2030-01-03T19:00:00": -22.22682657111151, "2030-01-03T20:00:00": -21.15242883824432, "2030-01-03T21:00:00": -21.63639801560457, "2030-01-03T22:00:00": -21.82255775160539, "2030-01-03T23:00:00": -22.223597034179647, "2030-01-04T00:00:00": -22.22682657111151, "2030-01-04T01:00:00": -22.22682657111151, "2030-01-04T02:00:00": -22.22682657111151, "2030-01-04T03:00:00": -22.194992564211734, "2030-01-04T04:00:00": -22.22682657111151, "2030-01-04T05:00:00": -22.680461169433208, "2030-01-04T06:00:00": -22.680461169433208, "2030-01-04T07:00:00": -22.680461169433208, "2030-01-04T08:00:00": -22.680461169433208, "2030-01-04T09:00:00": -22.22682657111151, "2030-01-04T10:00:00": -22.22682657111151, "2030-01-04T11:00:00": -22.55358650425297, "2030-01-04T12:00:00": -23.035479550728446, "2030-01-04T13:00:00": -22.49522415827006, "2030-01-04T14:00:00": -21.95427672218342, "2030-01-04T15:00:00": -19.374107394836244, "2030-01-04T16:00:00": -22.22682657111151, "2030-01-04T17:00:00": -22.22682657111151, "2030-01-04T18:00:00": -22.22682657111151, "2030-01-04T19:00:00": -22.22682657111151, "2030-01-04T20:00:00": -22.22682657111151, "2030-01-04T21:00:00": -22.22682657111151, "2030-01-04T22:00:00": -21.793376578613934, "2030-01-04T23:00:00": -22.22417373720319, "2030-01-05T00:00:00": -22.504682087856224, "2030-01-05T01:00:00": -22.504682087856224, "2030-01-05T02:00:00": -22.504682087856224, "2030-01-05T03:00:00": -22.504682087856224, "2030-01-05T04:00:00": -22.680461169433208, "2030-01-05T05:00:00": -25.900770852917027, "2030-01-05T06:00:00": -26.7259175390074, "2030-01-05T07:00:00": -30.176562410095706, "2030-01-05T08:00:00": -27.953948957347382, "2030-01-05T09:00:00": -22.680461169433208, "2030-01-05T10:00:00": -23.43075180306724, "2030-01-05T11:00:00": -23.03663295677554, "2030-01-05T12:00:00": -22.504682087856224, "2030-01-05T13:00:00": -22.22682657111151, "2030-01-05T14:00:00": -22.22682657111151, "2030-01-05T15:00:00": -22.22682657111151, "2030-01-05T16:00:00": -22.680461169433208, "2030-01-05T17:00:00": -22.680461169433208, "2030-01-05T18:00:00": -22.680461169433208, "2030-01-05T19:00:00": -22.680461169433208, "2030-01-05T20:00:00": -22.680461169433208, "2030-01-05T21:00:00": -22.22682657111151, "2030-01-05T22:00:00": -22.22682657111151, "2030-01-05T23:00:00": -22.22682657111151, "2030-01-06T00:00:00": -22.22682657111151, "2030-01-06T01:00:00": -22.22682657111151, "2030-01-06T02:00:00": -22.22682657111151, "2030-01-06T03:00:00": -22.22682657111151, "2030-01-06T04:00:00": -22.22682657111151, "2030-01-06T05:00:00": -22.680461169433208, "2030-01-06T06:00:00": -22.680461169433208, "2030-01-06T07:00:00": -22.680461169433208, "2030-01-06T08:00:00": -22.680461169433208, "2030-01-06T09:00:00": -32.3492333810049, "2030-01-06T10:00:00": -22.68565149664512, "2030-01-06T11:00:00": -24.018758205875283, "2030-01-06T12:00:00": -24.92314388740095, "2030-01-06T13:00:00": -22.680461169433208, "2030-01-06T14:00:00": -22.22682657111151, "2030-01-06T15:00:00": -22.22682657111151, "2030-01-06T16:00:00": -22.680461169433208, "2030-01-06T17:00:00": -22.680461169433208, "2030-01-06T18:00:00": -22.680461169433208, "2030-01-06T19:00:00": -22.680461169433208, "2030-01-06T20:00:00": -22.680461169433208, "2030-01-06T21:00:00": -22.680461169433208, "2030-01-06T22:00:00": -22.680461169433208, "2030-01-06T23:00:00": -22.22682657111151, "2030-01-07T00:00:00": -22.823714200482165, "2030-01-07T01:00:00": -22.680461169433208, "2030-01-07T02:00:00": -22.680461169433208, "2030-01-07T03:00:00": -22.680461169433208, "2030-01-07T04:00:00": -29.753377731417263, "2030-01-07T05:00:00": -29.763643045236392, "2030-01-07T06:00:00": -29.763643045236392, "2030-01-07T07:00:00": -29.763643045236392, "2030-01-07T08:00:00": -29.763643045236392, "2030-01-07T09:00:00": -30.18913453600902, "2030-01-07T10:00:00": -30.203898133411812, "2030-01-07T11:00:00": -30.204820858249487, "2030-01-07T12:00:00": -30.204820858249487, "2030-01-07T13:00:00": -30.204820858249487, "2030-01-07T14:00:00": -30.18302148395943, "2030-01-07T15:00:00": -29.75695329016325, "2030-01-07T16:00:00": -24.010684363545632, "2030-01-07T17:00:00": -23.43075180306724, "2030-01-07T18:00:00": -23.43075180306724, "2030-01-07T19:00:00": -23.46904488383073, "2030-01-07T20:00:00": -23.960280519287664, "2030-01-07T21:00:00": -24.23156162156395, "2030-01-07T22:00:00": -24.4515161547446, "2030-01-07T23:00:00": -22.680461169433208, "2030-01-08T00:00:00": -14.45217777007603, "2030-01-08T01:00:00": -14.45217777007603, "2030-01-08T02:00:00": -14.45217777007603, "2030-01-08T03:00:00": -14.45217777007603, "2030-01-08T04:00:00": -22.22682657111151, "2030-01-08T05:00:00": -22.680461169433208, "2030-01-08T06:00:00": -22.680461169433208, "2030-01-08T07:00:00": -22.680461169433208, "2030-01-08T08:00:00": -22.22682657111151, "2030-01-08T09:00:00": -22.22682657111151, "2030-01-08T10:00:00": -22.396607941243605, "2030-01-08T11:00:00": -22.680461169433208, "2030-01-08T12:00:00": -22.680461169433208, "2030-01-08T13:00:00": -23.134095767754907, "2030-01-08T14:00:00": -21.381495279197004, "2030-01-08T15:00:00": -19.23650605341804, "2030-01-08T16:00:00": -22.22682657111151, "2030-01-08T17:00:00": -22.22682657111151, "2030-01-08T18:00:00": -22.22682657111151, "2030-01-08T19:00:00": -22.22682657111151, "2030-01-08T20:00:00": -22.22682657111151, "2030-01-08T21:00:00": -20.617825135416687, "2030-01-08T22:00:00": -21.095681260727343, "2030-01-08T23:00:00": -21.095681260727343, "2030-01-09T00:00:00": -21.421057106612295, "2030-01-09T01:00:00": -21.512983568565613, "2030-01-09T02:00:00": -21.421057106612295, "2030-01-09T03:00:00": -22.22682657111151, "2030-01-09T04:00:00": -22.22682657111151, "2030-01-09T05:00:00": -22.680461169433208, "2030-01-09T06:00:00": -22.680461169433208, "2030-01-09T07:00:00": -22.680461169433208, "2030-01-09T08:00:00": -22.680461169433208, "2030-01-09T09:00:00": -22.22682657111151, "2030-01-09T10:00:00": -22.67204130528943, "2030-01-09T11:00:00": -29.709663642232435, "2030-01-09T12:00:00": -27.47736157868853, "2030-01-09T13:00:00": -22.680461169433208, "2030-01-09T14:00:00": -21.750585214266785, "2030-01-09T15:00:00": -10.804911770139771, "2030-01-09T16:00:00": -22.22682657111151, "2030-01-09T17:00:00": -22.22682657111151, "2030-01-09T18:00:00": -22.680461169433208, "2030-01-09T19:00:00": -22.680461169433208, "2030-01-09T20:00:00": -22.22682657111151, "2030-01-09T21:00:00": -22.22682657111151, "2030-01-09T22:00:00": -20.752312280507738, "2030-01-09T23:00:00": -21.509177328610203, "2030-01-10T00:00:00": -21.82255775160539, "2030-01-10T01:00:00": -21.82255775160539, "2030-01-10T02:00:00": -21.82255775160539, "2030-01-10T03:00:00": -21.82255775160539, "2030-01-10T04:00:00": -21.82255775160539, "2030-01-10T05:00:00": -22.22682657111151, "2030-01-10T06:00:00": -22.22682657111151, "2030-01-10T07:00:00": -22.22682657111151, "2030-01-10T08:00:00": -22.22682657111151, "2030-01-10T09:00:00": -23.134095767754907, "2030-01-10T10:00:00": -22.680461169433208, "2030-01-10T11:00:00": -22.680461169433208, "2030-01-10T12:00:00": -22.680461169433208, "2030-01-10T13:00:00": -22.680461169433208, "2030-01-10T14:00:00": -22.14470406055848, "2030-01-10T15:00:00": -19.305479735034204, "2030-01-10T16:00:00": -19.382757940189443, "2030-01-10T17:00:00": -22.22682657111151, "2030-01-10T18:00:00": -22.22682657111151, "2030-01-10T19:00:00": -22.22682657111151, "2030-01-10T20:00:00": -21.15242883824432, "2030-01-10T21:00:00": -21.63639801560457, "2030-01-10T22:00:00": -21.82255775160539, "2030-01-10T23:00:00": -22.223597034179647, "2030-01-11T00:00:00": -22.22682657111151, "2030-01-11T01:00:00": -22.22682657111151, "2030-01-11T02:00:00": -22.22682657111151, "2030-01-11T03:00:00": -22.194992564211734, "2030-01-11T04:00:00": -22.22682657111151, "2030-01-11T05:00:00": -22.680461169433208, "2030-01-11T06:00:00": -22.680461169433208, "2030-01-11T07:00:00": -22.680461169433208, "2030-01-11T08:00:00": -22.680461169433208, "2030-01-11T09:00:00": -22.22682657111151, "2030-01-11T10:00:00": -22.22682657111151, "2030-01-11T11:00:00": -22.55358650425297, "2030-01-11T12:00:00": -23.035479550728446, "2030-01-11T13:00:00": -22.49522415827006, "2030-01-11T14:00:00": -21.95427672218342, "2030-01-11T15:00:00": -19.374107394836244, "2030-01-11T16:00:00": -22.22682657111151, "2030-01-11T17:00:00": -22.22682657111151, "2030-01-11T18:00:00": -22.22682657111151, "2030-01-11T19:00:00": -22.22682657111151, "2030-01-11T20:00:00": -22.22682657111151, "2030-01-11T21:00:00": -22.22682657111151, "2030-01-11T22:00:00": -21.793376578613934, "2030-01-11T23:00:00": -22.22417373720319, "2030-01-12T00:00:00": -22.504682087856224, "2030-01-12T01:00:00": -22.504682087856224, "2030-01-12T02:00:00": -22.504682087856224, "2030-01-12T03:00:00": -22.504682087856224, "2030-01-12T04:00:00": -22.680461169433208, "2030-01-12T05:00:00": -25.900770852917027, "2030-01-12T06:00:00": -26.7259175390074, "2030-01-12T07:00:00": -30.176562410095706, "2030-01-12T08:00:00": -27.953948957347382, "2030-01-12T09:00:00": -22.680461169433208, "2030-01-12T10:00:00": -23.43075180306724, "2030-01-12T11:00:00": -23.03663295677554, "2030-01-12T12:00:00": -22.504682087856224, "2030-01-12T13:00:00": -22.22682657111151, "2030-01-12T14:00:00": -22.22682657111151, "2030-01-12T15:00:00": -22.22682657111151, "2030-01-12T16:00:00": -22.680461169433208, "2030-01-12T17:00:00": -22.680461169433208, "2030-01-12T18:00:00": -22.680461169433208, "2030-01-12T19:00:00": -22.680461169433208, "2030-01-12T20:00:00": -22.680461169433208, "2030-01-12T21:00:00": -22.22682657111151, "2030-01-12T22:00:00": -22.22682657111151, "2030-01-12T23:00:00": -22.22682657111151, "2030-01-13T00:00:00": -22.22682657111151, "2030-01-13T01:00:00": -22.22682657111151, "2030-01-13T02:00:00": -22.22682657111151, "2030-01-13T03:00:00": -22.22682657111151, "2030-01-13T04:00:00": -22.22682657111151, "2030-01-13T05:00:00": -22.680461169433208, "2030-01-13T06:00:00": -22.680461169433208, "2030-01-13T07:00:00": -22.680461169433208, "2030-01-13T08:00:00": -22.680461169433208, "2030-01-13T09:00:00": -32.3492333810049, "2030-01-13T10:00:00": -22.68565149664512, "2030-01-13T11:00:00": -24.018758205875283, "2030-01-13T12:00:00": -24.92314388740095, "2030-01-13T13:00:00": -22.680461169433208, "2030-01-13T14:00:00": -22.22682657111151, "2030-01-13T15:00:00": -22.22682657111151, "2030-01-13T16:00:00": -22.680461169433208, "2030-01-13T17:00:00": -22.680461169433208, "2030-01-13T18:00:00": -22.680461169433208, "2030-01-13T19:00:00": -22.680461169433208, "2030-01-13T20:00:00": -22.680461169433208, "2030-01-13T21:00:00": -22.680461169433208, "2030-01-13T22:00:00": -22.680461169433208, "2030-01-13T23:00:00": -22.22682657111151, "2030-01-14T00:00:00": -22.823714200482165, "2030-01-14T01:00:00": -22.680461169433208, "2030-01-14T02:00:00": -22.680461169433208, "2030-01-14T03:00:00": -22.680461169433208, "2030-01-14T04:00:00": -29.753377731417263, "2030-01-14T05:00:00": -29.763643045236392, "2030-01-14T06:00:00": -29.763643045236392, "2030-01-14T07:00:00": -29.763643045236392, "2030-01-14T08:00:00": -29.763643045236392, "2030-01-14T09:00:00": -30.18913453600902, "2030-01-14T10:00:00": -30.203898133411812, "2030-01-14T11:00:00": -30.204820858249487, "2030-01-14T12:00:00": -30.204820858249487, "2030-01-14T13:00:00": -30.204820858249487, "2030-01-14T14:00:00": -30.18302148395943, "2030-01-14T15:00:00": -29.75695329016325, "2030-01-14T16:00:00": -24.010684363545632, "2030-01-14T17:00:00": -23.43075180306724, "2030-01-14T18:00:00": -23.43075180306724, "2030-01-14T19:00:00": -23.46904488383073, "2030-01-14T20:00:00": -23.960280519287664, "2030-01-14T21:00:00": -24.23156162156395, "2030-01-14T22:00:00": -24.4515161547446, "2030-01-14T23:00:00": -22.680461169433208}}