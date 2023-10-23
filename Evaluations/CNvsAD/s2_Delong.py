
import numpy as np
import pandas as pd
from Utility.Training_Utilities import *
from Utility.DelongTest import delong_roc_test
import time
pd.options.mode.chained_assignment = None  # default='warn'

dpath = '/Volumes/JasonWork/Projects/ADNI_CSF/'
mydf = pd.read_csv(dpath + 'Results/FinalModel/CNvsAD/s0_predictions.csv')
cols_lst = mydf.columns.tolist()[4:]
nb_preds = len(cols_lst)

'''
outfile = dpath + 'Results/FinalModel/CNvsAD/s2_DelongTest_CNvsAD.csv'
mydf = mydf.loc[mydf.target_y != 'other']
mydf.reset_index(inplace = True, drop = True)
mydf['target_y'].replace(['CN A-T-', 'A+T+'], [0, 1], inplace = True)
'''

y_neg, y_pos = 'CN A-T-', 'CN A+T+'
y_neg, y_pos = 'MCI A-T-', 'MCI A+T+'
y_neg, y_pos = 'Dementia A-T-', 'Dementia A+T+'

outfile = dpath + 'Results/FinalModel/CNvsAD/s1_DelongTest_'+y_neg+'vs'+y_pos+'.csv'
mydf['target_y'] = mydf['DX_AT_status'].copy()
mydf = mydf.loc[(mydf.target_y == y_neg) | (mydf.target_y == y_pos)]
mydf.reset_index(inplace = True, drop = True)
mydf['target_y'].replace([y_neg, y_pos], [0, 1], inplace = True)
mydf['target_y'].value_counts()


delong_df = pd.DataFrame(np.zeros((nb_preds, nb_preds)))
delong_df.columns = cols_lst
delong_df.index = cols_lst

for i in range(nb_preds):
    for j in range(nb_preds):
        try:
            tmpdf = mydf[['target_y', cols_lst[i], cols_lst[j]]]
            tmpdf.dropna(how='any', inplace=True)
            tmpdf.reset_index(inplace=True, drop=True)
            log10_p = delong_roc_test(tmpdf.target_y, tmpdf.iloc[:, 1], tmpdf.iloc[:, 2])
            delong_df.iloc[i, j] = 10 ** log10_p[0][0]
            print(str(i) + ' ' + str(j))
        except:
            delong_df.iloc[i, j] = np.nan

delong_df.to_csv(outfile, index = True)


