from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
import pandas as pd
import os
import numpy as np


def get_datasets(data_directory='./'):
    path_by_dataset = {'organoids': os.path.join(data_directory, 'organoids_treutlein_dataset/RNA_ATAC_metacells_sce_peaks_obs2000_var4000.h5ad'),
                       'gbm': os.path.join(data_directory, 'gbm_multiome_scdori/23_01_23_atac_compressed_n1000.h5ad'),
                       'noack_2022': os.path.join(data_directory, 'noack_2022/*/merged_scATAC_integrated_cicero_faye_chong_obs*_var*.h5ad'),
                       'pancreatic_endocrinogenesis': os.path.join(data_directory, 'pancreatic_endocrinogenesis/*/pancreas_multiome_2022_processed_atac_obs*_var*.h5ad'),
                       'pbmc': os.path.join(data_directory, 'pbmc/*/processed_laura_2023_obs*_var*.h5ad')}
    return path_by_dataset

def iterative_metric_calc(pred, dataloader, function, lightweight_memory = True):
    if not lightweight_memory:    
        y = dataloader.dataset.rounds.copy()
        y_flatten = y.flatten()
        pred_flatten = pred.values.flatten().round()

        label_binarizer = LabelBinarizer().fit(np.concatenate((y_flatten, pred_flatten)))
        y_onehot_test = label_binarizer.transform(y_flatten)
        pred_onehot = label_binarizer.transform(pred_flatten)
        # y_onehot_test.shape  # (n_samples, n_classes)
        metric = function(y_onehot_test, pred_onehot)
        return metric
    else:
        print('calculate metric with iterative approach...')
        metric_sample = []
        pred_values = pred.values
        step = 1500
        for yi in range(0, dataloader.dataset.rounds.shape[0], step):
            print(yi, 'out of', pred_values.shape[0])            
            y = dataloader.dataset.rounds[yi: min(yi + step, dataloader.dataset.rounds.shape[0]),:].copy()                

            print(len(set(y.flatten())))
            
            y_flatten = y.flatten()

            pred_flatten = pred_values[yi: min(yi + step, pred_values.shape[0]),:].flatten().round()

            label_binarizer = LabelBinarizer().fit(np.concatenate((y_flatten, pred_flatten)))
            y_onehot_test = label_binarizer.transform(y_flatten)
            pred_onehot = label_binarizer.transform(pred_flatten)
            # y_onehot_test.shape  # (n_samples, n_classes)
            metric_sample.append(function(y_onehot_test, pred_onehot))
        return np.mean(metric_sample)
            

def get_auroc(model, dataloader):
    import mubind as mb
    pred = mb.tl.kmer_enrichment(model, dataloader)
    pred = pred[[c for c in pred if c.startswith('p')]]            
    
    y = dataloader.dataset.rounds.copy()
    y[y > 0] = 1
    return roc_auc_score(y.flatten(), pred.values.flatten())

def get_auprc(model, dataloader):
    import mubind as mb
    pred = mb.tl.kmer_enrichment(model, dataloader)
    pred = pred[[c for c in pred if c.startswith('p')]]            
    y = dataloader.dataset.rounds.copy()
    y[y > 0] = 1
    return average_precision_score(y.flatten(), pred.values.flatten())
    

from sklearn.preprocessing import LabelBinarizer
import sklearn
def get_auprc_multiclass(model, dataloader):
    # multiclass
    import mubind as mb
    pred = mb.tl.kmer_enrichment(model, dataloader)
    pred = pred[[c for c in pred if c.startswith('p')]]            
    pr_auc_multi= iterative_metric_calc(pred, dataloader, average_precision_score)    
    return pr_auc_multi

def get_auroc_multiclass(model, dataloader):
    # multiclass
    import mubind as mb
    pred = mb.tl.kmer_enrichment(model, dataloader)
    pred = pred[[c for c in pred if c.startswith('p')]]            
    roc_auc_multi= iterative_metric_calc(pred, dataloader, roc_auc_score)    
    return roc_auc_multi


def prepare_df(seqs, adata):
    # remove Ns
    for s in seqs:
        if 'N' in s:
            assert False
        # seqs = [[s[0], s[1].replace('N', '')] for s in seqs]
    counts = adata.X.A.T

    next_data = pd.DataFrame(counts)
    next_data['var'] = next_data.var(axis=1)
    next_data.index = [str(i) + '-' + s[1] for i, s in enumerate(seqs)]
    next_data.index.name = 'seq'
    print(next_data.shape)
    n_cells = adata.shape[0]
    n_peaks = adata.shape[1]
    top_var = next_data[['var']].sort_values('var', ascending=False).index[:n_peaks]
    # next_data = next_data.head(10000)
    next_data_sel = next_data.reindex(top_var) # .reset_index(drop=True)
    del next_data_sel['var']
    df = next_data_sel.copy() # sample
    zero_counts = df.sum(axis=1) == 0
    df = df[~zero_counts] # remove zeroes
    df.shape

    print('# cells', n_cells)
    print('# peaks', n_peaks)

    print('selected', df.shape)

    df.index = [v.split('-')[1] for v in df.index]
    df.index.name = 'seq'

    return df
