
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from Utility.DelongTest import delong_roc_test
from lightgbm import LGBMClassifier
from tqdm import tqdm
pd.options.mode.chained_assignment = None  # default='warn'

dpath = '/Volumes/JasonWork/Projects/ADNI_CSF/StrictCV/'
mydf = pd.read_csv(dpath + 'Data/ADNI_CSF_DATA.csv')
mydf_dict = pd.read_csv(dpath + 'Data/ADNI_CSF_DICT.csv', usecols = ['Analytes', 'Target', 'EntrezGeneSymbol', 'TargetFullName'])
mydf['target_y'] = mydf['DX1'].copy()
mydf = mydf.loc[mydf.target_y != 'other']
mydf['target_y'].replace(['CN A-T-', 'A+T+'], [0, 1], inplace = True)
cv_info_df = pd.read_csv(dpath + 'Data/CNvsAD_loocv_id.csv', usecols = ['RID', 'loocv_id', 'inner_cv_id'])
mydf = pd.merge(mydf, cv_info_df, how = 'inner', on = ['RID'])
loocv_id_lst = list(set(mydf.loocv_id))
inner_cv_id_lst = list(set(mydf.inner_cv_id))


for loocv_id in tqdm(loocv_id_lst):
    traindf = mydf.loc[mydf.loocv_id != loocv_id]
    traindf.reset_index(inplace = True, drop = True)
    my_f_df = pd.read_csv(dpath + 'Results/CNvsAD/sc_Importance/Test_loocv' + str(loocv_id) + '.csv')
    my_f_lst = my_f_df.Analytes.tolist()
    in_cv_y_pred_lst1 = np.zeros_like(traindf.target_y).tolist()
    in_cv_y_pred_lst2 = np.zeros_like(traindf.target_y).tolist()
    in_cv_y_pred_lst3 = np.zeros_like(traindf.target_y).tolist()
    tmp_f, cv_AUC_lst = [], []
    for f in my_f_lst:
        tmp_f.append(f)
        my_X = traindf[tmp_f]
        in_cv_AUC_lst, in_cv_y_true_lst, in_cv_y_pred_lst = [], [], []
        for inner_cv_id in inner_cv_id_lst:
            in_train_idx = traindf['inner_cv_id'].index[traindf['inner_cv_id'] != inner_cv_id]
            in_test_idx = traindf['inner_cv_id'].index[traindf['inner_cv_id'] == inner_cv_id]
            in_cv_X_train, in_cv_y_train = traindf.iloc[in_train_idx][tmp_f], traindf.iloc[in_train_idx].target_y
            in_cv_X_test, in_cv_y_test = traindf.iloc[in_test_idx][tmp_f], traindf.iloc[in_test_idx].target_y
            my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, verbosity=-1, seed=2023)
            my_lgb.fit(in_cv_X_train, in_cv_y_train)
            in_cv_y_pred = my_lgb.predict_proba(in_cv_X_test)[:, 1]
            in_cv_y_true_lst += in_cv_y_test.tolist()
            in_cv_y_pred_lst += in_cv_y_pred.tolist()
            in_cv_AUC_lst.append(roc_auc_score(in_cv_y_test, in_cv_y_pred))
        log10_p = delong_roc_test(np.array(in_cv_y_true_lst), np.array(in_cv_y_pred_lst2), np.array(in_cv_y_pred_lst))
        in_cv_y_pred_lst2 = in_cv_y_pred_lst1
        in_cv_y_pred_lst1 = in_cv_y_pred_lst
        cv_out = np.array([np.round(np.mean(in_cv_AUC_lst),3), np.round(np.std(in_cv_AUC_lst),3), 10 ** log10_p[0][0]])
        cv_AUC_lst.append(cv_out)
        print((f, cv_out))
    cv_AUC_df = pd.DataFrame(cv_AUC_lst, columns=['AUC_mean', 'AUC_std', 'Delong'])
    cv_AUC_df = pd.concat((pd.DataFrame({'Analytes': tmp_f}), cv_AUC_df), axis=1)
    cv_AUC_df = pd.merge(cv_AUC_df, my_f_df, how='left', on='Analytes')
    nb_f = 0
    while ((mydf.Delong[i] < 0.05) | (mydf.Delong[i + 1] < 0.05)):
        nb_f += 1
    cv_AUC_df['SelectedAnalytes'] = [1]*nb_f + [0]*(len(cv_AUC_df)-nb_f)
    results = cv_AUC_df[['Analytes', 'AUC_mean', 'AUC_std', 'Delong', 'SelectedAnalytes', 'Target', 'EntrezGeneSymbol', 'TargetFullName']]
    results.to_csv(dpath + 'Results/CNvsAD/sd_sf_selection/Test_loocv' + str(loocv_id) + '.csv', index=False)
    print('finished')


