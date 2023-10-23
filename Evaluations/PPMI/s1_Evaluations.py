
import numpy as np
import pandas as pd
from Utility.Training_Utilities import *
from sklearn.metrics import roc_auc_score
from joblib import Parallel, delayed
import time
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
        try:
            auc_lst.append(roc_auc_score(mydf_bt[y_true_col], mydf_bt[y_pred_col]))
        except:
            auc_lst.append(np.nan)
    median_auc = np.round(np.nanquantile(np.array(auc_lst),0.5), 5)
    lbd_auc = np.round(np.nanquantile(auc_lst, 0.025), 5)
    ubd_auc = np.round(np.nanquantile(auc_lst, 0.975), 5)
    out_auc = f'{median_auc:.3f}' + ' [' + f'{lbd_auc:.3f}' + '-' + f'{ubd_auc:.3f}' + ']'
    return (median_auc, lbd_auc, ubd_auc, out_auc)

def process(f, mydf):
    tmp_df = pd.DataFrame({'X': mydf[f], 'Y': mydf.target_y})
    tmp_df.dropna(how='any', inplace=True)
    tmp_df.reset_index(inplace=True, drop=True)
    tmp_auc = np.round(roc_auc_score(tmp_df['Y'], tmp_df['X']), 5)
    bt_auc = get_bt_output(tmp_df, 'Y', 'X', nb_iters=1000)
    print([f, tmp_auc] + list(bt_auc))
    return ([f, tmp_auc] + list(bt_auc))

dpath = '/Volumes/JasonWork/Projects/ADNI_CSF/'

mydf = pd.read_csv(dpath + 'Results/FinalModel/PPMI/s0_Predictions.csv')
outfile = dpath + 'Results/FinalModel/PPMI/s1_Evaluations.csv'
my_f_lst = mydf.columns.tolist()[3:]

results = []
for f in my_f_lst:
    try:
        results.append(process(f, mydf))
    except:
        pass

myout = pd.DataFrame(results)
myout.columns = ['Analytes', 'AUC_AllData', 'AUC_BT_median', 'AUC_BT_lbd', 'AUC_BT_ubd', 'AUC_BT_out']

myout.to_csv(outfile, index = False)

