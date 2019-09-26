from itertools import combinations
import numpy as np
import glob
import pandas as pd
import matplotlib.pyplot as plt

bio_conc_input = glob.glob(r"Biomarker_Concentration_Input/*.xlsx")
bio_conc_data = pd.read_excel(bio_conc_input[0]) 
data_matrix = np.ndarray.tolist(bio_conc_data.values)


max_biomarkers = #####
eczema_samples = #####
number_of_samples = ####
biomarker_eczema_STD = #[##, ## , ##,##, ##, ##, ##, ##, ##, ##, ##, ##, ##]#
eczema_mean = #[##, ## , ##,##, ##, ##, ##, ##, ##, ##, ##, ##, ##]#
exclusions_allowed = #
true_condition = np.zeros(number_of_samples)
for a in range(0,eczema_samples):
    true_condition[a] = 0.5

total_combos_iter = 0
combo_matrix = [0,0,0,0,0,0,0,0,0,0,0,0,0]
combo_list = []
combo_matrix_depths = [0,0,0,0,0,0,0,0,0,0,0,0,0]
sensi_by_combo = []
speci_by_combo = []

#Generate list of all combos for max_i Choose 1, max_i Choose 2, ..., max_i Choose max_i
for i in range(1,max_biomarkers):
    iter_combos = combinations(list(range(0,max_biomarkers)),i)
    total_combos_iter = total_combos_iter + len(list(iter_combos))
    combo_matrix[i] = [x for x in combinations(list(range(0,max_biomarkers),i)]
    combo_matrix_depths[i] = len(list(combinations(list(range(0,max_biomarkers),i)))
    for b in range(0,combo_matrix_depths[i]):
        combo_list.append(combo_matrix[i][b])

#Create Scoring Matrix to count true positive, false positive, true negative, and false negative for each combo in list
#Row in this matrix are sample number, columns are biomarker type, and field is 1 or 0 for each based on >2 standard deviation criterion
#Calculate key clinical sensitivity/specifcity performance indicators for every combo available.
for i in range(1,max_biomarkers):
    for j in range (0, len(combo_matrix[i])):
        posi_calc_matrix = np.ndarray.tolist(np.zeros((number_of_samples,max_cytokines)))
        global_positive_calc = np.zeros(number_of_samples)
        
        for l in combo_matrix[i][j]:
            tp_count = 0
            fp_count = 0
            tn_count = 0
            fn_count = 0
            for k in range(0, number_of_samples):
                if abs(data_matrix[k][l] - eczema_mean[l]) > 2*biomarker_eczema_STD[l]:
                    posi_calc_matrix[k][l] = 1.0
                if (sum(posi_calc_matrix[k]) >= i - exclusions_allowed) & (i - exclusions_allowed > 0):
                    global_positive_calc[k] = 1
                if true_condition[k] + global_positive_calc[k] == 1.5:
                    fp_count= fp_count + 1
                if true_condition[k] + global_positive_calc[k] == 1.0:
                    tp_count= tp_count + 1
                if true_condition[k] + global_positive_calc[k] == 0.5:
                    tn_count= tn_count + 1
                if true_condition[k] + global_positive_calc[k] == 0:
                    fn_count= fn_count + 1
        sensi_by_combo.append(tp_count/(tp_count+fn_count))
        speci_by_combo.append(tn_count/(tn_count+fp_count))
                
#ROC type plot compares sensitivity to specificity for every combo and places as scatter.
#Only points with high performance in both values are considered strong candidates for a potential clinical test.
plt.plot(speci_by_combo, sensi_by_combo, 'ro', markersize = 4)
plt.axis([0,1.05,0,1.05])
plt.xlabel('Specificity')
plt.ylabel('Sensitivity')
plt.title('Blank ROC Scatter for Biomarker Panel Combos')


plt.show()

