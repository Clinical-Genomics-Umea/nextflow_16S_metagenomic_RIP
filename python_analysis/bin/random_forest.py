import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import date
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_curve, roc_auc_score

'''
For each clinical variable: cancer/normal, cancer specific death, MSI/MSS, KRAS, BRAF:
    - Make a random forest classifier with 1000 trees, train/test sets are divided 80/20.
    - Plot a ROC curve and the top most informative (feature importance) ASVs.
'''

def clin_info(meta_file, clin_var):
    ''' Output clin data for samples included '''
    metadata = pd.DataFrame(pd.read_csv(meta_file, sep=',', header=0))
    clin_var = metadata[['Novogene_ID', clin_var]]
    return clin_var    


def make_predictions(infile, clin_var):
    ''' Random forest classifier with 1000 trees, train/test divided 80/20 '''
    indata = pd.DataFrame(pd.read_csv(infile, sep=',', header=0))
    indata = indata.T
    indata.columns = indata.iloc[0]
    indata = indata.drop(indata.index[0])
    indata = indata.drop(columns=['Age', 'Sex', 'Stage', 'Cancer'])
    indata.reset_index(inplace=True)
    indata = indata.rename(columns={'index': 'Novogene_ID'})
    indata = pd.merge(indata, clin_var, on='Novogene_ID')
    
    # drop samples with Nan in clinical variable
    filt = indata.columns[-1]
    indata = indata[indata[filt] != 'Nan']
      
    # Transform labels to numbers
    labelencoder = LabelEncoder()
    labels = labelencoder.fit_transform(indata.iloc[:, -1])
   
    pred_data = indata.drop(indata.columns[-1], axis=1).drop(indata.columns[0], axis=1)
    data_train, data_test, labels_train, labels_test = train_test_split(pred_data, labels, test_size=0.2, random_state=44)
    model = RandomForestClassifier(n_estimators=1000, max_features="sqrt", random_state=44, bootstrap=True)
    model.fit(data_train, labels_train)

    feat_importance = pd.DataFrame(model.feature_importances_)
    feat_importance.index = pred_data.columns
    pred_prob = model.predict_proba(data_test)[:, 1]
    fpr, tpr, thresholds = roc_curve(labels_test, pred_prob, pos_label=1)
    roc_auc = roc_auc_score(labels_test, pred_prob)  
    return fpr, tpr, roc_auc, feat_importance
    

   
def plot_importances(feat_importance, out_prefix, clin_var):
    ''' Plot feature importances of top 20 most informative ASVs'''
    feat_out = feat_importance.sort_values(by=0, ascending=False)
    feat_out.to_csv(str(out_prefix) + '_feat_importance_' + clin_var + '_' + date.today().strftime('%y%m%d') + '.csv')
    top = feat_importance.sort_values(by=0, ascending=True).tail(20)
    top.reset_index(inplace=True)
    names = top['index'].str.split('_').str[4:].apply('_'.join).str.split('.').str[0]
    top = top.drop(columns=['index'])
    top.index = names  
    top.plot(kind='barh', legend=False,  figsize=(10, 5))
    plt.xlabel('Feature importance')
    plt.ylabel(None)
    plt.subplots_adjust(left=0.4)
    plt.savefig(str(out_prefix) + '_feat_importances_' + clin_var + '_' + date.today().strftime('%y%m%d') + '.png', dpi=200)
    plt.close()
      
    
def plot_roc(fpr, tpr, roc_auc, out_prefix, clin_var):
    ''' Plot ROC curve '''
    plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], 'k--', label='Random classifier')
    plt.xlabel('False Positive Rate (1 - Specificity)')
    plt.ylabel('True Positive Rate (Sensitivity)')
    plt.legend(loc="lower right")
    plt.savefig(str(out_prefix) + '_roc_' + clin_var + '_'  + date.today().strftime('%y%m%d') + '.png', dpi=200)
    plt.close()
    


def main():
    infile = Path(sys.argv[1]) 
    meta_file = Path(sys.argv[2])
    out_prefix = Path(sys.argv[3])
    # Cancer vs normal
    cancer_info = clin_info(meta_file, 'Cancer')
    cancer_fpr, cancer_tpr, cancer_auc, cancer_feat_importances = make_predictions(infile, cancer_info)
    plot_importances(cancer_feat_importances, out_prefix, 'Cancer')
    plot_roc(cancer_fpr, cancer_tpr, cancer_auc, out_prefix, 'Cancer')
    # MSI vs MSS
    msi_info = clin_info(meta_file, 'MSI')
    msi_fpr, msi_tpr, msi_auc, msi_feat_importances = make_predictions(infile, msi_info)
    plot_importances(msi_feat_importances, out_prefix, 'MSI')
    plot_roc(msi_fpr, msi_tpr, msi_auc, out_prefix, 'MSI')    
    # Cancer specific death
    death_info = clin_info(meta_file, 'Can_spec_death')
    death_fpr, death_tpr, death_auc, death_feat_importances = make_predictions(infile, death_info)
    plot_importances(death_feat_importances, out_prefix, 'Can_spec_death')
    plot_roc(death_fpr, death_tpr, death_auc, out_prefix, 'Can_spec_death')   
    # KRAS
    kras_info = clin_info(meta_file, 'KRAS')
    kras_fpr, kras_tpr, kras_auc, kras_feat_importances = make_predictions(infile, kras_info)
    plot_importances(kras_feat_importances, out_prefix, 'KRAS')
    plot_roc(kras_fpr, kras_tpr, kras_auc, out_prefix, 'KRAS')    
    # BRAF
    braf_info = clin_info(meta_file, 'BRAF')
    braf_fpr, braf_tpr, braf_auc, braf_feat_importances = make_predictions(infile, braf_info)
    plot_importances(braf_feat_importances, out_prefix, 'BRAF')
    plot_roc(braf_fpr, braf_tpr, braf_auc, out_prefix, 'BRAF')


if __name__ == "__main__":
    main()      
