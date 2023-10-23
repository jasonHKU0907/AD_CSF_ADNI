
import numpy as np
import pandas as pd
from Utility.Training_Utilities import *
from lightgbm import LGBMClassifier
import warnings
import re
import shap
from sklearn.model_selection import StratifiedKFold
pd.options.mode.chained_assignment = None  # default='warn'

nb_folds = 10
dpath = '/Volumes/JasonWork/Projects/ADNI_CSF/'
outfile = dpath + 'Results/CNvsDementia/sa_Importance.csv'
mydf = pd.read_csv(dpath + 'Data/ADNI_CSF_DATA.csv')
auc_df = pd.read_csv(dpath + 'Results/CNvsDementia/s0.1_AUC_direct.csv', usecols=['Analytes', 'AUC_BT_median', 'AUC_BT_out'])
dict_df = pd.read_csv(dpath + 'Data/ADNI_CSF_DICT.csv', usecols = ['Analytes', 'Target', 'EntrezGeneSymbol', 'TargetFullName'])
#my_f_df = pd.read_csv(dpath + 'Data/CNvsAD_protein_remove_replication.csv')
#my_f_lst = my_f_df.Analytes.tolist()
#my_f_lst = mydf.columns.tolist()[13:]
my_f_df = pd.read_csv(dpath + 'Data/CNvsDementia_Association.csv')
my_f_df = my_f_df.loc[my_f_df.P_bonferroni < 0.05]
my_f_lst = my_f_df.Analytes.tolist()
mydf['target_y'] = mydf['DX3'].copy()
mydf = mydf.loc[mydf.target_y != 'MCI']
mydf.reset_index(inplace = True, drop = True)
mydf['target_y'].replace(['CN', 'Dementia'], [0, 1], inplace = True)

my_X = mydf[my_f_lst]
y = mydf.target_y
mykf = StratifiedKFold(n_splits = nb_folds, random_state = 2023, shuffle = True)

my_params = {'n_estimators': 500,
             'max_depth': 15,
             'num_leaves': 10,
             'subsample': 0.7,
             'learning_rate': 0.01,
             'colsample_bytree': 0.7}

def normal_imp(mydict):
    mysum = sum(mydict.values())
    mykeys = mydict.keys()
    for key in mykeys:
        mydict[key] = mydict[key]/mysum
    return mydict

tg_imp_cv = Counter()
tc_imp_cv = Counter()
shap_imp_cv = np.zeros(len(my_f_lst))

for train_idx, test_idx in mykf.split(my_X, y):
    X_train, X_test = my_X.iloc[train_idx,:], my_X.iloc[test_idx,:]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
    my_lgb = LGBMClassifier(objective = 'binary', metric = 'auc', is_unbalance = True, verbosity = 1, seed = 2023)
    my_lgb.set_params(**my_params)
    my_lgb.fit(X_train, y_train)
    totalgain_imp = my_lgb.booster_.feature_importance(importance_type='gain')
    totalgain_imp = dict(zip(my_lgb.booster_.feature_name(), totalgain_imp.tolist()))
    totalcover_imp = my_lgb.booster_.feature_importance(importance_type='split')
    totalcover_imp = dict(zip(my_lgb.booster_.feature_name(), totalcover_imp.tolist()))
    tg_imp_cv += Counter(normal_imp(totalgain_imp))
    tc_imp_cv += Counter(normal_imp(totalcover_imp))
    explainer = shap.TreeExplainer(my_lgb)
    shap_values = explainer.shap_values(X_test)
    shap_values = np.abs(np.average(shap_values[0], axis=0))
    shap_imp_cv += shap_values / np.sum(shap_values)


shap_imp_df = pd.DataFrame({'Analytes': my_f_lst,
                            'ShapValues_cv': shap_imp_cv/nb_folds})
shap_imp_df.sort_values(by = 'ShapValues_cv', ascending = False, inplace = True)

tg_imp_cv = normal_imp(tg_imp_cv)
tg_imp_df = pd.DataFrame({'Analytes': list(tg_imp_cv.keys()),
                          'TotalGain_cv': list(tg_imp_cv.values())})

tc_imp_cv = normal_imp(tc_imp_cv)
tc_imp_df = pd.DataFrame({'Analytes': list(tc_imp_cv.keys()),
                          'TotalCover_cv': list(tc_imp_cv.values())})

my_imp_df = pd.merge(left = shap_imp_df, right = tg_imp_df, how = 'left', on = ['Analytes'])
my_imp_df = pd.merge(left = my_imp_df, right = tc_imp_df, how = 'left', on = ['Analytes'])
my_imp_df['Ensemble_cv'] = (my_imp_df['ShapValues_cv'] + my_imp_df['TotalGain_cv'] + my_imp_df['TotalCover_cv'])/3
my_imp_df.sort_values(by = 'TotalGain_cv', ascending = False, inplace = True)

my_imp_df = pd.merge(my_imp_df, auc_df, how = 'left', on = 'Analytes')
my_imp_df = pd.merge(my_imp_df, dict_df, how = 'left', on = 'Analytes')
my_imp_df.to_csv(outfile, index = False)

print('finished')

