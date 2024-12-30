
ReadMe

Code:

  To run the code, please customize the working directory within each code file and ensure you have establish a Results folder named the same as ours. Please run the code following the sequence of files listed below. Notably, to ensure a strict leave-one-fold cross-validation, we recorded all statistics that each sub-folder contains 354 (CNvsAD) or 304 (CNvsDM) .csv files represents the leave-one-out results at each separate steps.
  
  Code of modeling contains a couple a sequential steps that strictly follows the descriptions from the manuscript as follows:
  
  "Those dysregulated proteins in biologically or clinically defined AD compared with controls were fed into a preliminarily trained LGBM classifier61 and ranked on the basis of their contributions to the cor- responding discriminative task (biologically defined AD versus CN A−T− and clinically diagnosed AD versus CN). Such contributions were quantified as information gains, which could be interpreted as the protein’s ability to distinguish AD (biologically/clinically defined) from controls. Top-ranked proteins that donated over 90% of overall information gains were retained. Afterwards, hierarchical clustering on the Spearman rank-order correlations was conducted to further remove the redundant proteins with multicollinearity; the proteins with the largest AUC were chosen from each cluster where highly cor- related ones were grouped together. A sequential forward selection strategy was then employed to rerank the preselected proteins by feeding them into newly developed LGBM classifiers two at a time in succession, according to the updated importance ranking orders. The stopping point was determined by the optimal performance of the AUC, artificially defined as when no incremental performance was noticed in two consecutive DeLong tests. The model performance did not improve significantly in this case, even with the addition of proteins.
  Next, receiver operating characteristic analyses were performed within ADNI and PPMI cohorts to assess the discriminative accuracy of the above-selected important CSF proteins alone or combined with demographic measures, including age, sex, education and APOE ε4 alleles. Performance was evaluated through a nested leave-one-out cross-validation (LOOCV) procedure where the model was trained using all individuals except one left as a test case, and the procedure was repeated until all participants had been used as test cases. On the basis of the predicted risk probabilities of each individual obtained through LOOCV, we then conducted a bootstrap strategy with 1,000 iterations by drawing samples with replacement and evaluating their perfor- mance; then, AUC statistics of medians and 95% CIs were reported. To evaluate the integrated ability of the proteins, we further leveraged the protein panels, a combination of the selected proteins described above, as an alternative to single proteins for modelling. Differences in the AUCs between various models were estimated using DeLong tests."
  
  - CNvsAD: code of CN A−T− versus biological defined AD. Besides brief descriptions listed below, please also see notes within the code files for step-to-step illustration.

    - s00_association_analysis.txt: R code that perform association analysis of each protein using regressions, corresponding to Results/CNvsAD/s00_Associations.

    - s01_AUC_calculation.py: Caculated AUC of each protein through 1000 iterated bootstrap, corresponding to Results/CNvsAD/s01_AUC_direct.

    - sa_ImpRanking.py: Calculate the information gain represents protein importance, corresponding to Results/CNvsAD/sa_Importance.

    - sb_RmMultiColin.py: Removal multi-collinear proteins using hiearical clustering based on Spearman correlations, in a subset of top 90% important proteins from sa_Importance, corresponding to Results/CNvsAD/sb_rmMultiColinearity.

    - sc_ImpRanking.py: Recalculated CSF protein importance, corresponding to Results/CNvsAD/sc_Importance.

    - sd_SFS.py: Sequential forward selection of CSF proteins to determine the protein panel, corresponding to Results/CNvsAD/sd_sf_selection.

    - se_Predictions.py: Calculate the predicted probabilities based on proteins panel and perform AUC evaluation using bootstap with 1,000 iterations, corresponding to Results/CNvsAD/ se0_loocv_prediction.csv and se1_loocv_evaluation.csv.


  - CNvsDM: code of target information of CN from clinically diagnosed AD dementia. Besides brief descriptions listed below, please also see notes within the code files for step-to-step illustration.

    - s00_association_analysis.txt: R code that perform association analysis of each protein using regressions, corresponding to Results/CNvsDM/s00_Associations.

    - s01_AUC_calculation.py: Caculated AUC of each protein through 1000 iterated bootstrap, corresponding to Results/CNvsDM/s01_AUC_direct.

    - sa_ImpRanking.py: Calculate the information gain represents protein importance, corresponding to Results/CNvsDM/sa_Importance.

    - sb_RmMultiColin.py: Removal multi-collinear proteins using hiearical clustering based on Spearman correlations, in a subset of top 90% important proteins from sa_Importance, corresponding to Results/CNvsDM/sb_rmMultiColinearity.

    - sc_ImpRanking.py: Recalculated CSF protein importance, corresponding to Results/CNvsDM/sc_Importance.

    - sd_SFS.py: Sequential forward selection of CSF proteins to determine the protein panel, corresponding to Results/CNvsDM/sd_sf_selection.

    - se_Predictions.py: Calculate the predicted probabilities based on proteins panel and perform AUC evaluation using bootstap with 1,000 iterations, corresponding to Results/CNvsDM/ se0_loocv_prediction.csv and se1_loocv_evaluation.csv.





Results:

  - CNvsAD: result of CN A−T− versus biological defined AD, each sub-folders contain 354 .csv files represents the leave-one-out results at each separate steps.

    - s00_Associations: Association results from regression analysis, corresponding to Code/CNvsAD/s00_association_analysis.txt.

    - s01_AUC_direct: AUCs of each individual CSF protein caculated through 1000 iterated Bootstrap, corresponding to Code/CNvsAD/s01_AUC_calculation.py.

    - sa_Importance: CSF protein importance based on adjusted significant proteins in s00_Associations, corresponding to Code/CNvsAD/sa_ImpRanking.py.

    - sb_rmMultiColinearity: Removal multi-collinear proteins based on top 90% important proteins from sa_Importance, corresponding to Code/CNvsAD/sb_RmMultiColin.py.

    - sc_Importance: Recalculated CSF protein importance based on proteins selected from sb_rmMultiColinearity, corresponding to Code/CNvsAD/sc_ImpRanking.py.

    - sd_sf_selection: Sequential forward selection procedure of CSF protein based on sc_Importance, corresponding to Code/CNvsAD/sd_SFS.py.

    - se0_loocv_prediction.csv: Predicted probabilities based on selected proteins from sd_sf_selection, corresponding to Code/CNvsAD/se_Predictions.py.

    - se1_loocv_evaluation.csv: Evaluation of predicted probabilities calculated in se0_loocv_prediction.csv, corresponding to Code/CNvsAD/se_Predictions.py.


  - CNvsDM: result of CN from clinically diagnosed AD dementia, each sub-folders contain 304 .csv files represents the leave-one-out results at each separate steps.

    - s00_Associations: Association results from regression analysis, corresponding to Code/CNvsDM/s00_association_analysis.txt.

    - s01_AUC_direct: AUCs of each individual CSF protein caculated through 1000 iterated Bootstrap, corresponding to Code/CNvsDM/s01_AUC_calculation.py.

    - sa_Importance: CSF protein importance based on adjusted significant proteins in s00_Associations, corresponding to Code/CNvsDM/sa_ImpRanking.py.

    - sb_rmMultiColinearity: Removal multi-collinear proteins based on top 90% important proteins from sa_Importance, corresponding to Code/CNvsDM/sb_RmMultiColin.py.

    - sc_Importance: Recalculated CSF protein importance based on proteins selected from sb_rmMultiColinearity, corresponding to Code/CNvsDM/sc_ImpRanking.py.

    - sd_sf_selection: Sequential forward selection procedure of CSF protein based on sc_Importance, corresponding to Code/CNvsDM/sd_SFS.py.

    - se0_loocv_prediction.csv: Predicted probabilities based on selected proteins from sd_sf_selection, corresponding to Code/CNvsDM/se_Predictions.py.

    - se1_loocv_evaluation.csv: Evaluation of predicted probabilities calculated in se0_loocv_prediction.csv, corresponding to Code/CNvsDM/se_Predictions.py.

