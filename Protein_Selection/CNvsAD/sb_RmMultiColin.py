

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Utility.Training_Utilities import *
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from collections import defaultdict

def best_analysts(my_X, f_lst, df_auc):
    f_df = pd.DataFrame({'Analytes':  my_X.columns[f_lst]})
    merged_df = pd.merge(f_df, df_auc, how='inner', on=['Analytes'])
    merged_df.sort_values(by = 'AUC_BT_median', ascending=False)
    return merged_df.Analytes[0]


def get_imp_analy(Imp_df, top_prop):
    imp_score, iter = 0, 0
    while imp_score < top_prop:
        imp_score += Imp_df.TotalGain_cv[iter]
        iter+=1
    return iter+1

top_prop=0.9

dpath = '/Volumes/JasonWork/Projects/ADNI_CSF/'
outfile = dpath + 'Results/CNvsAD/sb_rmMultiColinearity.csv'
outimg = dpath + 'Results/CNvsAD/sb_rmMultiColinearity.png'
Imp_df = pd.read_csv(dpath + 'Results/CNvsAD/sa_Importance.csv')
top_nb = get_imp_analy(Imp_df, top_prop)
Imp_df = Imp_df.iloc[:top_nb,:]
my_f_lst = Imp_df.Analytes.tolist()
mydf = pd.read_csv(dpath + 'Data/ADNI_CSF_DATA.csv', usecols = ['id', 'RID', 'DX1'] + my_f_lst)
mydf['target_y'] = mydf['DX1'].copy()
mydf = mydf.loc[mydf.target_y != 'other']
mydf.reset_index(inplace = True, drop = True)
mydf['target_y'].replace(['CN A-T-', 'A+T+'], [0, 1], inplace = True)

my_X = mydf[my_f_lst]
y = mydf.target_y
my_label = Imp_df.EntrezGeneSymbol

corr = np.array(my_X.corr(method='spearman'))
#corr = np.array(my_X.corr(method='pearson'))
corr = np.nan_to_num(corr)
corr = (corr + corr.T) / 2
np.fill_diagonal(corr, 1)
distance_matrix = 1 - np.abs(corr)
dist_linkage = hierarchy.ward(squareform(distance_matrix))


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 10))
dendro = hierarchy.dendrogram(dist_linkage, labels=my_label.tolist(), ax=ax2)
ax2.set_xticklabels(dendro["ivl"], rotation=60, fontsize=12, horizontalalignment='right')
ax2.axhline(y = 0.5, color = 'r', linewidth=2, linestyle = '--')
dendro_idx = np.arange(0, len(dendro["ivl"]))
ax1.imshow(corr[dendro["leaves"], :][:, dendro["leaves"]])
ax1.set_xticks(dendro_idx)
ax1.set_yticks(dendro_idx)
ax1.set_xticklabels(dendro["ivl"], rotation=60, fontsize=12, horizontalalignment='right')
ax1.set_yticklabels(dendro["ivl"], fontsize=12)
fig.tight_layout()
plt.show()
plt.savefig(outimg)

cluster_ids = hierarchy.fcluster(dist_linkage, .5, criterion="distance")
cluster_id_to_feature_ids = defaultdict(list)
for idx, cluster_id in enumerate(cluster_ids):
    cluster_id_to_feature_ids[cluster_id].append(idx)

Imp_df['Cluster_ids'] = cluster_ids
selected_f = [best_analysts(my_X, v, Imp_df) for v in cluster_id_to_feature_ids.values()]
select = ['*' if Analytes in selected_f else '' for Analytes in Imp_df.Analytes]
Imp_df['Selected'] = select

myout_df = Imp_df[['Analytes', 'TotalGain_cv', 'AUC_BT_median', 'AUC_BT_out', 'Cluster_ids', 'Selected', 'Target', 'EntrezGeneSymbol', 'TargetFullName']]
myout_df.to_csv(outfile, index = False)

print('Finished')
