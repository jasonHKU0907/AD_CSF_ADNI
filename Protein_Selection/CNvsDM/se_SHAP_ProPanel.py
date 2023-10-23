
import numpy as np
import pandas as pd
from Utility.Training_Utilities import *
from lightgbm import LGBMClassifier
from sklearn.calibration import CalibratedClassifierCV
import shap
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
pd.options.mode.chained_assignment = None  # default='warn'

def get_nb_f(mydf):
    p_lst = mydf.Delong2.tolist()
    i = 0
    while((p_lst[i]<0.05)|(p_lst[i+1]<0.05)):
        i+=1
    return i

nb_folds = 10

dpath = '/Volumes/JasonWork/Projects/ADNI_CSF/'
outimg = dpath + 'Results/CNvsDementia/sh_SHAP_ProPanel.png'
pro_auc_df = pd.read_csv(dpath + 'Results/CNvsDementia/sd_SFS.csv')
nb_f = get_nb_f(pro_auc_df)
pro_f_lst = pro_auc_df.Analytes.tolist()[:nb_f]
pro_f_name_lst = pro_auc_df.EntrezGeneSymbol.tolist()[:nb_f]
f_name_dict = dict(zip(pro_f_lst, pro_f_name_lst))
pro_df = pd.read_csv(dpath + 'Data/ADNI_CSF_DATA.csv', usecols=['RID'] + pro_f_lst)
pro_df.rename(columns = f_name_dict, inplace=True)


mydf = pd.read_csv(dpath + 'Data/ADNI_CSF_DATA.csv', usecols=['RID', 'DX3'])
mydf['target_y'] = mydf['DX3'].copy()
mydf = mydf.loc[mydf.target_y != 'MCI']
mydf.reset_index(inplace = True, drop = True)
mydf['target_y'].replace(['CN', 'Dementia'], [0, 1], inplace = True)

mydf = pd.merge(mydf, pro_df, how = 'left', on = ['RID'])
tmp_f_lst = pro_f_name_lst

mykf = StratifiedKFold(n_splits = nb_folds, random_state = 2023, shuffle = True)

my_params = {'n_estimators': 500,
             'max_depth': 15,
             'num_leaves': 10,
             'subsample': 0.7,
             'learning_rate': 0.01,
             'colsample_bytree': 0.7}


shap_cv_lst, test_idx_lst = [], []

for train_idx, test_idx in mykf.split(mydf, mydf.target_y):
    X_train, X_test = mydf[tmp_f_lst].iloc[train_idx, :], mydf[tmp_f_lst].iloc[test_idx, :]
    y_train, y_test = mydf.target_y.iloc[train_idx], mydf.target_y.iloc[test_idx]
    my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, n_jobs=4, verbosity=-1, seed=2023)
    my_lgb.set_params(**my_params)
    my_lgb.fit(X_train, y_train)
    explainer = shap.Explainer(my_lgb)
    shap_val = explainer(X_test)
    shap_cv_lst.append(shap_val[:,:,1])
    test_idx_lst.append(test_idx)


shap_values = shap_cv_lst[0]

for i in range(1, nb_folds):
    shap_values.values = np.concatenate((shap_values.values, shap_cv_lst[i].values), axis = 0)
    shap_values.base_values = np.concatenate((shap_values.base_values, shap_cv_lst[i].base_values), axis = 0)
    shap_values.data = np.concatenate((shap_values.data, shap_cv_lst[i].data), axis = 0)


nb2plot = len(tmp_f_lst)
shap.plots.beeswarm(shap_values, max_display=nb2plot, order=list(np.linspace(0, nb2plot, nb2plot+1).astype('uint8')))
plt.gcf().set_size_inches(12, 6)
ax = plt.gca()
ax.set_ylabel('Selected proteins', fontsize = 18, weight = 'bold')
ax.set_xlabel('SHAP values', fontsize = 14, weight = 'bold')
ax.tick_params(axis='x', labelsize=14)
ylabels = [tick.get_text() for tick in ax.get_yticklabels()]
ax.set_yticklabels(ylabels, fontsize = 16, color = 'black')
plt.tight_layout()

plt.savefig(outimg, dpi=300)


