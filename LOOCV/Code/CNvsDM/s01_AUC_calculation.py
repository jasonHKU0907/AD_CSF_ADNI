
import numpy as np
import pandas as pd
import random
from sklearn.metrics import roc_auc_score
from joblib import Parallel, delayed
from tqdm import tqdm
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
    median_auc = np.round(np.quantile(np.array(auc_lst),0.5), 10)
    lbd_auc = np.round(np.quantile(auc_lst, 0.025), 10)
    ubd_auc = np.round(np.quantile(auc_lst, 0.975), 10)
    out_auc = f'{median_auc:.3f}' + ' [' + f'{lbd_auc:.3f}' + '-' + f'{ubd_auc:.3f}' + ']'
    return (median_auc, lbd_auc, ubd_auc, out_auc)

def process(mydf, f):
    tmp_df = pd.DataFrame({'X': mydf[f], 'Y': mydf.target_y})
    nb_all = len(tmp_df)
    tmp_df.dropna(how='any', inplace=True)
    tmp_df.reset_index(inplace=True, drop=True)
    nb_tmp = len(tmp_df)
    tmp_auc = roc_auc_score(tmp_df['Y'], tmp_df['X'])
    if tmp_auc < 0.5:
        tmp_df['X'] = 1 - (tmp_df['X'] - tmp_df['X'].min()) / (tmp_df['X'].max() - tmp_df['X'].min())
    elif tmp_auc >= 0.5:
        tmp_df['X'] = (tmp_df['X'] - tmp_df['X'].min()) / (tmp_df['X'].max() - tmp_df['X'].min())
    tmp_auc = roc_auc_score(tmp_df['Y'], tmp_df['X'])
    bt_auc = get_bt_output(tmp_df, 'Y', 'X', nb_iters=1000)
    print([f, tmp_auc] + list(bt_auc))
    na_prop = (nb_all - nb_tmp)/nb_all*100
    na_prop= f'{na_prop:.1f}' + '%'
    return ([f, na_prop, tmp_auc] + list(bt_auc))

dpath = '/Volumes/JasonWork/Projects/ADNI_CSF/StrictCV/'
outfile = dpath + 'Results/CNvsAD/s0.1_AUC_direct.csv'
mydf = pd.read_csv(dpath + 'Data/ADNI_CSF_DATA.csv')
mydf_dict = pd.read_csv(dpath + 'Data/ADNI_CSF_DICT.csv', usecols = ['Analytes', 'Target', 'EntrezGeneSymbol', 'TargetFullName'])
mydf['target_y'] = mydf['DX3'].copy()
mydf = mydf.loc[mydf.target_y != 'MCI']
mydf['target_y'].replace(['CN', 'Dementia'], [0, 1], inplace = True)
cv_info_df = pd.read_csv(dpath + 'Data/CNvsDM_loocv_id.csv', usecols = ['RID', 'loocv_id', 'inner_cv_id'])
mydf = pd.merge(mydf, cv_info_df, how = 'inner', on = ['RID'])
loocv_id_lst = list(set(mydf.loocv_id))

for loocv_id in tqdm(loocv_id_lst):
    traindf = mydf.loc[mydf.loocv_id != loocv_id]
    traindf.reset_index(inplace = True, drop = True)
    my_f_df = pd.read_csv(dpath + 'Results/CNvsDM/s00_Associations/Test_loocv'+str(loocv_id)+'.csv')
    my_f_df = my_f_df.loc[my_f_df.P_BH_adjust < 0.05]
    my_f_lst = my_f_df.Analytes
    results = Parallel(n_jobs=20)(delayed(process)(traindf, f) for f in my_f_lst)
    results = pd.DataFrame(results)
    results.columns = ['Analytes', 'NA_prop', 'AUC_AllData', 'AUC_BT_median', 'AUC_bt_lbd', 'AUC_bt_ubd', 'AUC_bt_out']
    results = pd.merge(results[['Analytes', 'AUC_bt_out']], mydf_dict, how = 'left', on = 'Analytes')
    results.to_csv(dpath + 'Results/CNvsDM/s01_AUC_direct/Test_loocv' + str(loocv_id) + '.csv', index = False)
