
import pandas as pd
from lightgbm import LGBMClassifier
from tqdm import tqdm
from collections import Counter
pd.options.mode.chained_assignment = None  # default='warn'

def normal_imp(mydict):
    mysum = sum(mydict.values())
    mykeys = mydict.keys()
    for key in mykeys:
        mydict[key] = mydict[key]/mysum
    return mydict

dpath = '/Volumes/JasonWork/Projects/ADNI_CSF/StrictCV/'
mydf = pd.read_csv(dpath + 'Data/ADNI_CSF_DATA.csv')
mydf_dict = pd.read_csv(dpath + 'Data/ADNI_CSF_DICT.csv', usecols = ['Analytes', 'Target', 'EntrezGeneSymbol', 'TargetFullName'])
mydf['target_y'] = mydf['DX3'].copy()
mydf = mydf.loc[mydf.target_y != 'MCI']
mydf['target_y'].replace(['CN', 'Dementia'], [0, 1], inplace = True)
cv_info_df = pd.read_csv(dpath + 'Data/CNvsDM_loocv_id.csv', usecols = ['RID', 'loocv_id', 'inner_cv_id'])
mydf = pd.merge(mydf, cv_info_df, how = 'inner', on = ['RID'])
loocv_id_lst = list(set(mydf.loocv_id))
inner_cv_id_lst = list(set(mydf.inner_cv_id))

for loocv_id in tqdm(loocv_id_lst):
    traindf = mydf.loc[mydf.loocv_id != loocv_id]
    traindf.reset_index(inplace = True, drop = True)
    analy_auc_df = pd.read_csv(dpath + 'Results/CNvsAD/s01_AUC_direct/Test_loocv' + str(loocv_id) + '.csv')
    my_f_df = pd.read_csv(dpath + 'Results/CNvsDM/sb_rmMultiColinearity/Test_loocv' + str(loocv_id) + '.csv')
    my_f_df = my_f_df.loc[my_f_df.Selected == '*']
    my_f_lst = my_f_df.Analytes
    tg_imp_cv = Counter()
    for inner_cv_id in inner_cv_id_lst:
        in_train_idx = traindf['inner_cv_id'].index[traindf['inner_cv_id'] != inner_cv_id]
        in_cv_X_train, in_cv_y_train = traindf.iloc[in_train_idx][my_f_lst], traindf.iloc[in_train_idx].target_y
        my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, verbosity=-1, seed=2023)
        my_lgb.fit(in_cv_X_train, in_cv_y_train)
        totalgain_imp = my_lgb.booster_.feature_importance(importance_type='gain')
        totalgain_imp = dict(zip(my_lgb.booster_.feature_name(), totalgain_imp.tolist()))
        tg_imp_cv += Counter(normal_imp(totalgain_imp))
    tg_imp_cv = normal_imp(tg_imp_cv)
    tg_imp_df = pd.DataFrame({'Analytes': list(tg_imp_cv.keys()), 'Importance': list(tg_imp_cv.values())})
    tg_imp_df.sort_values(by='Importance', ascending=False, inplace=True)
    tg_imp_df = pd.merge(tg_imp_df, analy_auc_df[['Analytes', 'AUC_bt_out']], how='left', on=['Analytes'])
    tg_imp_df = pd.merge(tg_imp_df, mydf_dict, how='left', on=['Analytes'])
    tg_imp_df.to_csv(dpath + 'Results/CNvsDM/sc_Importance/Test_loocv' + str(loocv_id) + '.csv', index = False)


