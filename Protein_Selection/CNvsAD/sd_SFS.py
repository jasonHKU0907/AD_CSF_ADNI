
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from Utility.Training_Utilities import *
from Utility.DelongTest import delong_roc_test
from lightgbm import LGBMClassifier
import scipy.stats as stats
from sklearn.model_selection import StratifiedKFold
pd.options.mode.chained_assignment = None  # default='warn'

ImpMethod = 'TotalGain'
nb_folds = 10
dpath = '/Volumes/JasonWork/Projects/ADNI_CSF/'
outfile = dpath + 'Results/CNvsAD/sd_SFS.csv'
mydf = pd.read_csv(dpath + 'Data/ADNI_CSF_DATA.csv')
imp_f_df = pd.read_csv(dpath + 'Results/CNvsAD/sc_ReImportance.csv')
imp_f_df.sort_values(by = ImpMethod + '_cv', ascending=False, inplace = True)
imp_f_lst = imp_f_df.Analytes.tolist()

mydf['target_y'] = mydf['DX1'].copy()
mydf = mydf.loc[mydf.target_y != 'other']
mydf.reset_index(inplace = True, drop = True)
mydf['target_y'].replace(['CN A-T-', 'A+T+'], [0, 1], inplace = True)

y = mydf.target_y
mykf = StratifiedKFold(n_splits = nb_folds, random_state = 2023, shuffle = True)


my_params = {'n_estimators': 500,
             'max_depth': 15,
             'num_leaves': 10,
             'subsample': 0.7,
             'learning_rate': 0.01,
             'colsample_bytree': 0.7}

y_pred_lst_prev1, y_pred_lst_prev2, y_pred_lst_prev3 = np.zeros_like(y).tolist(), np.zeros_like(y).tolist(), np.zeros_like(y).tolist()

tmp_f, AUC_cv_lst= [], []
for f in imp_f_lst:
    tmp_f.append(f)
    my_X = mydf[tmp_f]
    AUC_cv, y_pred_lst, y_true_lst = [], [], []
    for train_idx, test_idx in mykf.split(my_X, y):
        X_train, X_test = my_X.iloc[train_idx, :], my_X.iloc[test_idx, :]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, n_jobs=4, verbosity=-1, seed=2023)
        my_lgb.set_params(**my_params)
        my_lgb.fit(X_train, y_train)
        y_pred_prob = my_lgb.predict_proba(X_test)[:, 1]
        AUC_cv.append(roc_auc_score(y_test, y_pred_prob))
        y_pred_lst += y_pred_prob.tolist()
        y_true_lst += y_test.tolist()
    auc_full = roc_auc_score(y_true_lst, y_pred_lst)
    log10_p1 = delong_roc_test(np.array(y_true_lst), np.array(y_pred_lst_prev1), np.array(y_pred_lst))
    log10_p2 = delong_roc_test(np.array(y_true_lst), np.array(y_pred_lst_prev2), np.array(y_pred_lst))
    log10_p3 = delong_roc_test(np.array(y_true_lst), np.array(y_pred_lst_prev3), np.array(y_pred_lst))
    y_pred_lst_prev3 = y_pred_lst_prev2
    y_pred_lst_prev2 = y_pred_lst_prev1
    y_pred_lst_prev1 = y_pred_lst
    tmp_out = np.array([np.mean(AUC_cv), np.std(AUC_cv), 10**log10_p1[0][0], 10**log10_p2[0][0], 10**log10_p3[0][0], auc_full])
    AUC_cv_lst.append(tmp_out)
    print((f, tmp_out))

AUC_df = pd.DataFrame(AUC_cv_lst, columns = ['AUC_mean', 'AUC_std', 'Delong1','Delong2','Delong3', 'AUC_all'])
AUC_df[['AUC_mean', 'AUC_std', 'AUC_all']] = np.round(AUC_df[['AUC_mean', 'AUC_std', 'AUC_all']], 3)
AUC_df = pd.concat((pd.DataFrame({'Analytes':tmp_f}), AUC_df), axis = 1)
AUC_df = pd.merge(AUC_df, imp_f_df, how = 'left', on ='Analytes')

AUC_df = AUC_df[['Analytes', 'AUC_mean', 'AUC_std', 'Delong1','Delong2','Delong3', 'AUC_all', 'Target', 'EntrezGeneSymbol', 'TargetFullName']]
AUC_df.to_csv(outfile, index = False)
print('finished')


