
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from Utility.Training_Utilities import *
from Utility.DelongTest import delong_roc_test
from lightgbm import LGBMClassifier
import scipy.stats as stats
from sklearn.model_selection import StratifiedKFold
pd.options.mode.chained_assignment = None  # default='warn'

nb_folds = 10

mykf = StratifiedKFold(n_splits = nb_folds, random_state = 2023, shuffle = True)

my_params = {'n_estimators': 500,
             'max_depth': 15,
             'num_leaves': 10,
             'subsample': 0.7,
             'learning_rate': 0.01,
             'colsample_bytree': 0.7}


def get_pred_probs_loocv(mydf, tmp_f_lst, my_params, colname):
    y_pred_lst, RID_lst = [], []
    for i in range(len(mydf)):
        train_idx = [ele for ele in range(len(mydf))]
        train_idx.remove(i)
        test_idx = [i]
        X_train, X_test = mydf[tmp_f_lst].iloc[train_idx, :], mydf[tmp_f_lst].iloc[test_idx, :]
        y_train, y_test = mydf.target_y.iloc[train_idx], mydf.target_y.iloc[test_idx]
        my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, n_jobs=4, verbosity=-1, seed=2023)
        my_lgb.set_params(**my_params)
        my_lgb.fit(X_train, y_train)
        y_pred_prob = np.round(my_lgb.predict_proba(X_test)[:, 1], 10)
        y_pred_lst.append(y_pred_prob[0])
        RID_lst.append(mydf.RID.iloc[i])
    pred_df = pd.DataFrame({'RID': RID_lst, colname: y_pred_lst})
    return pred_df

def get_pred_probs_all(mydf0, mydf1, tmp_f_lst, mykf, my_params, colname):
    X_train, X_test = mydf0[tmp_f_lst], mydf1[tmp_f_lst]
    y_train = mydf0.target_y
    my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, n_jobs=4, verbosity=-1, seed=2023)
    my_lgb.set_params(**my_params)
    my_lgb.fit(X_train, y_train)
    y_pred_prob = np.round(my_lgb.predict_proba(X_test)[:, 1], 10)
    pred_df1 = pd.DataFrame({'RID':mydf1.RID, colname: y_pred_prob})
    #pred_df0 = get_pred_probs_cv(mydf0, tmp_f_lst, mykf, my_params, colname)
    pred_df0 = get_pred_probs_loocv(mydf0, tmp_f_lst, my_params, colname)
    pred_df = pd.concat([pred_df0, pred_df1], axis = 0)
    return pred_df


dpath = '/Volumes/JasonWork/Projects/ADNI_CSF/'
outfile = dpath + 'Results/FinalModel/nonADvsAD/s0_predictions.csv'

target_df = pd.read_csv(dpath + 'Data/Data0921/ADNI_Target.csv', usecols=['RID', 'autopsy_DX'])
target_df['target_y'] = target_df['autopsy_DX'].copy()
target_df = target_df.loc[(target_df.target_y == 'Non-AD')|(target_df.target_y == 'AD')]
target_df.reset_index(inplace = True, drop = True)
target_df['target_y'].replace(['Non-AD', 'AD'], [0, 1], inplace = True)

pro_f_lst1 = ['YWHAG', 'SMOC1', 'PIGR', 'TMOD2', 'ACHE', 'PCSK1', 'MMP10', 'IRF1']
pro_f_lst2 = ['YWHAZ', 'YWHAB']
pro_df = pd.read_csv(dpath + 'Data/Data0921/ADNI_CSF_Data.csv', usecols=['RID'] + pro_f_lst1+pro_f_lst2)

bio_f_lst = ['CSF_ABETA42', 'CSF_PTAU181', 'CSF_TTAU', 'CSF_ABETA42B40',
             'CSF_PTAUBABETA42', 'CSF_TTAUBABETA42', 'CSFNFL', 'CSF_sTREM2',
             'CSF_PGRN', 'PLASMAPTAU181', 'PLASMA_NFL']

bio_df = pd.read_csv(dpath + 'Data/Data0921/ADNI_Biomarker.csv', usecols=['RID'] + bio_f_lst)

cov_f_lst = ['Age', 'Sex', 'Education', 'ApoE4']
cov_df = pd.read_csv(dpath + 'Data/Data0921/ADNI_Covariates.csv', usecols=['RID'] + cov_f_lst)

mydf = pd.merge(target_df, pro_df, how = 'left', on = ['RID'])
mydf = pd.merge(mydf, bio_df, how = 'left', on = ['RID'])
mydf = pd.merge(mydf, cov_df, how = 'left', on = ['RID'])

pred_df_CNvsAD = pd.read_csv(dpath + 'Results/FinalModel/CNvsAD/s0_predictions.csv', usecols = ['RID', 'ProPanel', 'ProPanel_Cov'])
pred_df_CNvsAD.rename(columns={'ProPanel':'ProPanel_CNvsAD', 'ProPanel_Cov':'ProPanel_Cov_CNvsAD'}, inplace = True)
pred_df_CNvsDementia = pd.read_csv(dpath + 'Results/FinalModel/CNvsDementia/s0_predictions.csv', usecols = ['RID', 'ProPanel', 'ProPanel_Cov'])
pred_df_CNvsDementia.rename(columns={'ProPanel':'ProPanel_CNvsDementia', 'ProPanel_Cov':'ProPanel_Cov_CNvsDementia'}, inplace = True)


pred_cov = get_pred_probs_loocv(mydf, cov_f_lst, my_params, 'Cov')
pred_bio_comb = get_pred_probs_loocv(mydf, ['CSF_ABETA42', 'CSF_PTAU181', 'CSF_TTAU'], my_params, 'Comb_AB42_PT181_TTAU')

myout_df = pd.merge(target_df, pred_cov, how = 'left', on = 'RID')
myout_df = pd.merge(myout_df, pred_bio_comb, how = 'left', on = 'RID')
myout_df = pd.merge(myout_df, pred_df_CNvsAD, how = 'left', on = 'RID')
myout_df = pd.merge(myout_df, pred_df_CNvsDementia, how = 'left', on = 'RID')


for i in range(len(pro_f_lst1)):
    tmp_f = pro_f_lst1[i]
    tmp_df = mydf[['target_y', tmp_f] + cov_f_lst]
    tmp_df.dropna(how='any', inplace=True)
    tmp_df.reset_index(inplace=True, drop=True)
    if roc_auc_score(tmp_df.target_y, tmp_df[tmp_f])<0.5:
        mydf[tmp_f] = -mydf[tmp_f]
    else:
        pass
    myout_df = pd.merge(myout_df, mydf[['RID', tmp_f]], how='left', on='RID')
    pred_pro_cov = get_pred_probs_loocv(mydf, [tmp_f]+cov_f_lst, my_params, tmp_f+'_Cov')
    myout_df = pd.merge(myout_df, pred_pro_cov, how='left', on='RID')

single_pro_lst = bio_f_lst+pro_f_lst2

for i in range(len(single_pro_lst)):
    tmp_f = single_pro_lst[i]
    tmp_df = mydf[['target_y', tmp_f]]
    tmp_df.dropna(how='any', inplace=True)
    tmp_df.reset_index(inplace=True, drop=True)
    if roc_auc_score(tmp_df.target_y, tmp_df[tmp_f])<0.5:
        mydf[tmp_f] = -mydf[tmp_f]
    else:
        pass
    myout_df = pd.merge(myout_df, mydf[['RID', tmp_f]], how='left', on='RID')

myout_df.to_csv(outfile, index = False)


