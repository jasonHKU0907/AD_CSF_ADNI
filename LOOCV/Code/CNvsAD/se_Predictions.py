
import numpy as np
import pandas as pd
import random
from tqdm import tqdm
from sklearn.metrics import  roc_auc_score
from lightgbm import LGBMClassifier
pd.options.mode.chained_assignment = None  # default='warn'

def get_bt_output(mydf, y_true_col, y_pred_col, nb_iters):
    auc_lst = []
    idx_lst = [ele for ele in range(len(mydf))]
    for i in range(nb_iters):
        random.seed(i)
        bt_idx = [random.choice(idx_lst) for _ in range(len(idx_lst))]
        mydf_bt = mydf.copy()
        mydf_bt = mydf_bt.iloc[bt_idx, :]
        mydf_bt.reset_index(inplace = True, drop = True)
        auc_lst.append(roc_auc_score(mydf_bt[y_true_col], mydf_bt[y_pred_col]))
    median_auc = np.round(np.quantile(np.array(auc_lst),0.5), 5)
    lbd_auc = np.round(np.quantile(auc_lst, 0.025), 5)
    ubd_auc = np.round(np.quantile(auc_lst, 0.975), 5)
    out_auc = f'{median_auc:.3f}' + ' [' + f'{lbd_auc:.3f}' + '-' + f'{ubd_auc:.3f}' + ']'
    return (median_auc, lbd_auc, ubd_auc, out_auc)

dpath = '/Volumes/JasonWork/Projects/ADNI_CSF/StrictCV/'
mydf = pd.read_csv(dpath + 'Data/ADNI_CSF_DATA.csv')
covdf = pd.read_csv(dpath + 'Data/ADNI_Covariates.csv')
cov_f_lst = ['Age', 'Sex', 'Education', 'ApoE4']
mydf = pd.merge(mydf, covdf, how = 'left', on = ['RID'])
mydf['target_y'] = mydf['DX1'].copy()
mydf = mydf.loc[mydf.target_y != 'other']
mydf['target_y'].replace(['CN A-T-', 'A+T+'], [0, 1], inplace = True)
cv_info_df = pd.read_csv(dpath + 'Data/CNvsAD_loocv_id.csv', usecols = ['RID', 'loocv_id', 'inner_cv_id'])
mydf = pd.merge(mydf, cv_info_df, how = 'inner', on = ['RID'])
loocv_id_lst = list(set(mydf.loocv_id))

cv_rid_lst, cv_y_true_lst  = [], []
cv_y_pred_pro_lst, cv_y_pred_cov_lst, cv_y_pred_pro_cov_lst = [], [], []

for loocv_id in tqdm(loocv_id_lst):
    traindf = mydf.loc[mydf.loocv_id != loocv_id]
    traindf.reset_index(inplace = True, drop = True)
    testdf = mydf.loc[mydf.loocv_id == loocv_id]
    testdf.reset_index(inplace = True, drop = True)
    auc_imp_df = pd.read_csv(dpath + 'Results/CNvsAD/sd_sf_selection/Test_loocv' + str(loocv_id) + '.csv')
    nb_f = auc_imp_df.SelectedAnalytes.sum()
    pro_f_lst = auc_imp_df.Analytes.tolist()[:nb_f]
    cv_y_train = traindf.target_y
    cv_pro_train, cv_pro_test = traindf[pro_f_lst], testdf[pro_f_lst]
    my_lgb_pro = LGBMClassifier(objective='binary', metric='auc', is_unbalance=False, verbosity=-1, seed=2022)
    my_lgb_pro.fit(cv_pro_train, cv_y_train)
    cv_y_pred_pro_lst.append(my_lgb_pro.predict_proba(cv_pro_test)[:,1][0])
    cv_cov_train, cv_cov_test = traindf[cov_f_lst], testdf[cov_f_lst]
    my_lgb_cov = LGBMClassifier(objective='binary', metric='auc', is_unbalance=False, verbosity=-1, seed=2022)
    my_lgb_cov.fit(cv_cov_train, cv_y_train)
    cv_y_pred_cov_lst.append(my_lgb_cov.predict_proba(cv_cov_test)[:,1][0])
    cv_pro_cov_train, cv_pro_cov_test = traindf[pro_f_lst+cov_f_lst], testdf[pro_f_lst+cov_f_lst]
    my_lgb_pro_cov = LGBMClassifier(objective='binary', metric='auc', is_unbalance=False, verbosity=-1, seed=2022)
    my_lgb_pro_cov.fit(cv_pro_cov_train, cv_y_train)
    cv_y_pred_pro_cov_lst.append(my_lgb_pro_cov.predict_proba(cv_pro_cov_test)[:,1][0])
    cv_rid_lst.append(testdf.RID.iloc[0])
    cv_y_true_lst.append(testdf.target_y.iloc[0])

pred_df = pd.DataFrame({'RID':cv_rid_lst, 'target_y':cv_y_true_lst,
                        'y_pred_pro':cv_y_pred_pro_lst,
                        'y_pred_cov':cv_y_pred_cov_lst,
                        'y_pred_pro_cov':cv_y_pred_pro_cov_lst})

pred_df = pd.merge(pred_df, mydf[['RID', 'loocv_id', 'inner_cv_id']], how = 'left', on = 'RID')
pred_df.to_csv(dpath + 'Results/CNvsAD/se0_loocv_prediction.csv', index = False)

pro_eval = get_bt_output(pred_df, 'target_y', 'y_pred_pro', nb_iters = 1000)
cov_eval = get_bt_output(pred_df, 'target_y', 'y_pred_cov', nb_iters = 1000)
pro_cov_eval = get_bt_output(pred_df, 'target_y', 'y_pred_pro_cov', nb_iters = 1000)
out_df = pd.DataFrame([pro_eval, cov_eval, pro_cov_eval])
out_df.columns = ['AUC_bt_median', 'AUC_bt_lower', 'AUC_bt_upper', 'AUC_bt_out']
out_df['Analytes'] = ['ProteinPanel', 'Covariates', 'ProteinPanel+Covariates']
out_df.to_csv(dpath + 'Results/CNvsAD/se1_loocv_evaluation.csv', index = False)

