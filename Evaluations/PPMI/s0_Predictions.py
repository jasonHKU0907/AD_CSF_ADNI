
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from Utility.Training_Utilities import *
from Utility.DelongTest import delong_roc_test
from lightgbm import LGBMClassifier
import scipy.stats as stats
from sklearn.model_selection import StratifiedKFold
pd.options.mode.chained_assignment = None  # default='warn'

dpath = '/Volumes/JasonWork/Projects/ADNI_CSF/'
outfile = dpath + 'Results/FinalModel/PPMI/s0_Predictions.csv'

cov_f_lst = ['Age', 'Sex', 'Education', 'ApoE4']
bio_f_lst = ['abeta', 'asyn', 'ttau', 'ptau', 'ptauBabeta42', 'ttauBabeta42', 'abetaBasyn', 'ptauBasyn', 'ttauBasyn']
pro_f_lst1 = ['YWHAG', 'SMOC1', 'PIGR', 'TMOD2']
pro_f_lst2 = ['ACHE', 'YWHAG', 'MMP10', 'PCSK1', 'IRF1']
pro_f_lst3 = ['YWHAZ', 'YWHAB']

pro_f_lst = list(set(pro_f_lst1+pro_f_lst2))
read_f_lst = cov_f_lst+bio_f_lst+pro_f_lst+pro_f_lst3

mydf = pd.read_csv(dpath + 'Data/Data0921/PPMI_CSF_Data.csv', usecols=['PATNO', 'diagnosis1'] + read_f_lst)

mydf['target_y'] = mydf['diagnosis1'].copy()
mydf = mydf.loc[(mydf.target_y == 'Control_A-T-')|(mydf.target_y == 'Control_A+T+')]
mydf.reset_index(inplace = True, drop = True)
mydf['target_y'].replace(['Control_A-T-', 'Control_A+T+'], [0, 1], inplace = True)

my_params = {'n_estimators': 1000,
             'max_depth': 30,
             'num_leaves': 30,
             'subsample': 0.7,
             'learning_rate': 0.01,
             'colsample_bytree': 0.7}


def get_pred_probs(mydf, tmp_f_lst, my_params, colname):
    y_pred_lst, PATNO_lst = [], []
    for i in range(len(mydf)):
        train_idx = [ele for ele in range(len(mydf))]
        train_idx.remove(i)
        test_idx = [i]
        X_train, X_test = mydf[tmp_f_lst].iloc[train_idx, :], mydf[tmp_f_lst].iloc[test_idx, :]
        y_train, y_test = mydf.target_y.iloc[train_idx], mydf.target_y.iloc[test_idx]
        my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, n_jobs=4, verbosity=-1, seed=2023)
        #my_lgb.set_params(**my_params)
        my_lgb.fit(X_train, y_train)
        y_pred_prob = np.round(my_lgb.predict_proba(X_test)[:, 1], 10)
        y_pred_lst.append(y_pred_prob[0])
        PATNO_lst.append(mydf.PATNO.iloc[i])
    pred_df = pd.DataFrame({'PATNO': PATNO_lst, colname: y_pred_lst})
    return pred_df

pred_pro_CNvsAD = get_pred_probs(mydf, pro_f_lst1, my_params, 'ProPanel_CNvsAD')
pred_pro_CNvsAD_Cov = get_pred_probs(mydf, pro_f_lst1+cov_f_lst, my_params, 'ProPanel_CNvsAD_Cov')
pred_pro_CNvsDementia = get_pred_probs(mydf, pro_f_lst2, my_params, 'ProPanel_CNvsDementia')
pred_pro_CNvsDementia_Cov = get_pred_probs(mydf, pro_f_lst2+cov_f_lst, my_params, 'ProPanel_CNvsDementia_Cov')
pred_cov = get_pred_probs(mydf, cov_f_lst, my_params, 'Cov')
pred_bio_comb = get_pred_probs(mydf, ['abeta', 'ptau', 'ttau'], my_params, 'Comb_AB42_PT181_TTAU')

myout_df = pd.merge(mydf[['PATNO', 'diagnosis1', 'target_y']], pred_pro_CNvsAD, how = 'left', on = 'PATNO')
myout_df = pd.merge(myout_df, pred_pro_CNvsAD_Cov, how = 'left', on = 'PATNO')
myout_df = pd.merge(myout_df, pred_pro_CNvsDementia, how = 'left', on = 'PATNO')
myout_df = pd.merge(myout_df, pred_pro_CNvsDementia_Cov, how = 'left', on = 'PATNO')
myout_df = pd.merge(myout_df, pred_cov, how = 'left', on = 'PATNO')
myout_df = pd.merge(myout_df, pred_bio_comb, how = 'left', on = 'PATNO')

for i in range(len(pro_f_lst)):
    tmp_df = mydf[['target_y', pro_f_lst[i]]]
    tmp_df.dropna(how='any', inplace=True)
    tmp_df.reset_index(inplace=True, drop=True)
    if roc_auc_score(tmp_df.target_y, tmp_df[pro_f_lst[i]])<0.5:
        mydf[pro_f_lst[i]] = -mydf[pro_f_lst[i]]
    else:
        pass
    raw_pro_tmp = mydf[['PATNO'] + [pro_f_lst[i]]]
    myout_df = pd.merge(myout_df, raw_pro_tmp, how='left', on='PATNO')
    pred_pro_cov_tmp = get_pred_probs(mydf, [pro_f_lst[i]]+cov_f_lst, my_params, pro_f_lst[i]+'_Cov')
    myout_df = pd.merge(myout_df, pred_pro_cov_tmp, how='left', on='PATNO')

alone_f_lst = pro_f_lst3+bio_f_lst

for i in range(len(alone_f_lst)):
    tmp_df = mydf[['target_y', alone_f_lst[i]]]
    tmp_df.dropna(how='any', inplace=True)
    tmp_df.reset_index(inplace=True, drop=True)
    if roc_auc_score(tmp_df.target_y, tmp_df[alone_f_lst[i]])<0.5:
        mydf[alone_f_lst[i]] = -mydf[alone_f_lst[i]]
    else:
        pass
    raw_bio_tmp = mydf[['PATNO'] + [alone_f_lst[i]]]
    myout_df = pd.merge(myout_df, raw_bio_tmp, how='left', on='PATNO')

myout_df.to_csv(outfile, index = False)


