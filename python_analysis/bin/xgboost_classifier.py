import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import date
from pathlib import Path
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_curve, roc_auc_score

np.set_printoptions(linewidth=np.inf)


''' For each clinical variable: cancer/normal, cancer specific death, MSI/MSS, KRAS, BRAF:
    - Make a XGBoost classifier, train/test sets are divided 80/20.
    - Plot a ROC curve and the top most informative (feature importance) ASVs.'''


def clin_info(meta_file, clin_var):
    ''' Output clin data for samples included '''
    metadata = pd.DataFrame(pd.read_csv(meta_file, sep=',', header=0))
    clin_var_info = metadata[['Novogene_ID', clin_var]]
    return clin_var_info

def make_predictions(infile, clin_var_info):
    ''' Random forest classifier with 1000 trees, train/test divided 80/20 '''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    indata = indata.T
    indata.columns = indata.iloc[0]
    indata = indata.drop(indata.index[0])
    indata = indata.drop(columns=['Age', 'Sex', 'Stage', 'Cancer'])
    indata.reset_index(inplace=True)
    indata = indata.rename(columns={'index': 'Novogene_ID'})
    indata = pd.merge(indata, clin_var_info, on='Novogene_ID')

    # drop samples with Nan in clinical variable
    filt = indata.columns[-1]
    indata = indata[indata[filt] != 'Nan']
   
    # Transform labels to numbers
    labelencoder = LabelEncoder()
    labels = labelencoder.fit_transform(indata.iloc[:, -1])
    pred_data = indata.drop(indata.columns[-1], axis=1).drop(indata.columns[0], axis=1)
    data_train, data_test, labels_train, labels_test = train_test_split(pred_data, labels, test_size=0.2, random_state=44)
    data_train = data_train.astype('float')
    data_test = data_test.astype('float')
    data_for_mean = data_test.copy()
    data_for_mean['clin_var_info'] = labels_test
    data_for_mean.to_csv('data_for_mean_MSI_'  + date.today().strftime('%y%m%d') + '.csv')
    means_df = data_for_mean.astype(float).groupby('clin_var_info').mean()
    feat_colors = []
    for _ , values in means_df.items():
        if values[0] > values[1]:
            feat_colors.append('cornflowerblue')
        elif values[1] > values[0]:
            feat_colors.append('indianred')
        elif values[0] == values[1]:
            feat_colors.append('gray')

    weights = [1, 10, 25, 50, 100, 200, 300, 500, 1000, 2000, 3000]
    roc_auc_list = []               
    for w in weights:
        model = XGBClassifier(scale_pos_weight=w)
        model.fit(data_train, labels_train)
        pred_prob = model.predict_proba(data_test)[:, 1]
        roc_auc = roc_auc_score(labels_test, pred_prob)
        roc_auc_list.append(roc_auc)

    max_index = roc_auc_list.index(max(roc_auc_list))
    best_w = weights[max_index]
    print("best weight for " + str(filt) + " is " + str(best_w))

    model = XGBClassifier(scale_pos_weight=best_w)
    model.fit(data_train, labels_train)
    pred_prob = model.predict_proba(data_test)[:, 1]
    roc_auc = roc_auc_score(labels_test, pred_prob)
    feat_importance = pd.DataFrame(model.feature_importances_)
    feat_importance.index = pred_data.columns
    fpr, tpr, thresholds = roc_curve(labels_test, pred_prob, pos_label=1)    
    return fpr, tpr, roc_auc, feat_importance, feat_colors

    
   
def plot_importances(feat_importance, out_prefix, clin_var, feat_colors):
    ''' Plot feature importances of top 20 most informative ASVs '''
    feat_importance['col'] = feat_colors
    feat_out = feat_importance.sort_values(by=0, ascending=False)
    #feat_out.to_csv(str(out_prefix) + '_feat_importance_' + clin_var + '_' + date.today().strftime('%y%m%d') + '.csv')
    top = feat_importance.sort_values(by=0, ascending=True).tail(20)
    top.reset_index(inplace=True)
    #names = top['index'].str.split('_').str[3:].apply('_'.join).str.split('.').str[0]
    names = top['index'].str.split('_').str[3:]
    top = top.drop(columns=['index'])
    top.index = names
    top[0].plot(kind='barh', legend=False,  figsize=(12, 5), color=top['col'])
    plt.xlabel('Feature importance')
    plt.ylabel(None)
    plt.subplots_adjust(left=0.6)
    #plt.savefig(str(out_prefix) + '_feat_importances_' + clin_var + '_' + date.today().strftime('%y%m%d') + '.png', dpi=200)
    plt.show()
    plt.close()
      
    
def plot_roc(fpr, tpr, roc_auc, out_prefix, clin_var):
    ''' Plot ROC curve '''
    plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], 'k--', label='Random classifier')
    plt.xlabel('False Positive Rate (1 - Specificity)')
    plt.ylabel('True Positive Rate (Sensitivity)')
    plt.legend(loc="lower right")
    plt.savefig(str(out_prefix) + '_roc_' + clin_var + '_'  + date.today().strftime('%y%m%d') + '.png', dpi=200)
    plt.show()
    plt.close()
    


def main():
    infile = Path(sys.argv[1])
    meta_file = Path(sys.argv[2])
    out_prefix = Path(sys.argv[3])
 
    # Cancer vs normal
    cancer_info = clin_info(meta_file, 'Cancer')
    cancer_fpr, cancer_tpr, cancer_auc, cancer_feat_importances, cancer_colors = make_predictions(infile, cancer_info)
    plot_importances(cancer_feat_importances, out_prefix, 'Cancer', cancer_colors)
    plot_roc(cancer_fpr, cancer_tpr, cancer_auc, out_prefix, 'Cancer')
    
    # #MSI vs MSS
    msi_info = clin_info(meta_file, 'MSI')
    msi_fpr, msi_tpr, msi_auc, msi_feat_importances, msi_colors = make_predictions(infile, msi_info)
    plot_importances(msi_feat_importances, out_prefix, 'MSI', msi_colors)
    plot_roc(msi_fpr, msi_tpr, msi_auc, out_prefix, 'MSI') 
       
    # # Cancer specific death
    # death_info = clin_info(meta_file, 'Can_spec_death')
    # death_fpr, death_tpr, death_auc, death_feat_importances, death_colors = make_predictions(infile, death_info)
    # plot_importances(death_feat_importances, out_prefix, 'Can_spec_death', death_colors)
    # plot_roc(death_fpr, death_tpr, death_auc, out_prefix, 'Can_spec_death')   
    # # KRAS
    # kras_info = clin_info(meta_file, 'KRAS')
    # kras_fpr, kras_tpr, kras_auc, kras_feat_importances, kras_colors = make_predictions(infile, kras_info)
    # plot_importances(kras_feat_importances, out_prefix, 'KRAS', kras_colors)
    # plot_roc(kras_fpr, kras_tpr, kras_auc, out_prefix, 'KRAS')    
    # # BRAF
    # braf_info = clin_info(meta_file, 'BRAF')
    # braf_fpr, braf_tpr, braf_auc, braf_feat_importances, braf_colors = make_predictions(infile, braf_info)
    # plot_importances(braf_feat_importances, out_prefix, 'BRAF', braf_colors)
    # plot_roc(braf_fpr, braf_tpr, braf_auc, out_prefix, 'BRAF')


if __name__ == "__main__":
    main()      
