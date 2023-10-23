

import numpy as np
import pandas as pd
from Utility.Training_Utilities import *
from Utility.DelongTest import delong_roc_test
from lightgbm import LGBMClassifier
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
pd.options.mode.chained_assignment = None  # default='warn'

def get_nb_f(mydf):
    p_lst = mydf.Delong2.tolist()
    auc_lst = mydf.AUC_mean.tolist()
    i = 0
    while((p_lst[i]<0.05)&(auc_lst[i+1]>auc_lst[i])):
        i+=1
    return i+1

def get_nb_f(mydf):
    p_lst = mydf.Delong2.tolist()
    i = 0
    while((p_lst[i]<0.05)|(p_lst[i+1]<0.05)):
        i+=1
    return i


dpath = '/Volumes/JasonWork/Projects/ADNI_CSF/'
outimg = dpath + 'Results/CNvsDementia/sd_SFS_plot.png'
pro_imp_df = pd.read_csv(dpath + 'Results/CNvsDementia/sc_ReImportance.csv', usecols = ['Analytes', 'TotalGain_cv'])
pro_imp_df.rename(columns = {'TotalGain_cv': 'Pro_imp'}, inplace = True)
pro_auc_df = pd.read_csv(dpath + 'Results/CNvsDementia/sd_SFS.csv')
mydf = pd.merge(pro_imp_df, pro_auc_df, how = 'left', on = ['Analytes'])
mydf['AUC_lower'] = mydf['AUC_mean'] - mydf['AUC_std']
mydf['AUC_upper'] = mydf['AUC_mean'] + mydf['AUC_std']
mydf['AUC_upper'].iloc[mydf['AUC_upper']>=1] = 1
mydf['pro_idx'] = [i for i in range(1, len(mydf)+1)]
nb_f = get_nb_f(mydf)


fig, ax = plt.subplots(figsize = (18, 6.5))
palette = sns.color_palette("Blues",n_colors=len(mydf))
palette.reverse()
sns.barplot(ax=ax, x = "EntrezGeneSymbol", y = "Pro_imp", palette=palette, data=mydf.sort_values(by="Pro_imp", ascending=False))
y_imp_up_lim = round(mydf['Pro_imp'].max() + 0.01, 2)
ax.set_ylim([0, y_imp_up_lim])
ax.tick_params(axis='y', labelsize=14)
ax.set_xticklabels(mydf['EntrezGeneSymbol'], rotation=45, fontsize=10, horizontalalignment='right')
my_col = ['r']*nb_f + ['k']*(len(mydf)-nb_f)
for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), my_col):
    ticklabel.set_color(tickcolor)

ax.set_ylabel('Protein importance', weight='bold', fontsize=18)
#ax.set_title(my_title, y=1.0, pad=-25, weight='bold', fontsize=24)
ax.set_xlabel('')
ax.grid(which='minor', alpha=0.2, linestyle=':')
ax.grid(which='major', alpha=0.5,  linestyle='--')
ax.set_axisbelow(True)

ax2 = ax.twinx()
ax2.plot(np.arange(nb_f+1), mydf['AUC_mean'][:nb_f+1], 'red', alpha = 0.8, marker='o')
ax2.plot(np.arange(nb_f+1, len(mydf)), mydf['AUC_mean'][nb_f+1:], 'black', alpha = 0.8, marker='o')
ax2.plot([nb_f, nb_f+1], mydf['AUC_mean'][nb_f:nb_f+2], 'black', alpha = 0.8, marker='o')
plt.fill_between(mydf['pro_idx']-1, mydf['AUC_lower'], mydf['AUC_upper'], color = 'tomato', alpha = 0.2)
ax2.set_ylabel('Cumulative AUC', weight='bold', fontsize=18)
ax2.tick_params(axis='y', labelsize=14)
y_auc_up_lim = round(mydf['AUC_upper'].max() + 0.01, 2)
y_auc_low_lim = round(mydf['AUC_lower'].min() - 0.01, 2)
ax2.set_ylim([y_auc_low_lim, y_auc_up_lim])

fig.tight_layout()
plt.xlim([-.6, len(mydf)-.2])
plt.savefig(outimg, dpi = 200)
