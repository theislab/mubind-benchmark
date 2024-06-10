from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
import pandas as pd

def get_datasets():
    path_by_dataset = {'organoids': '/mnt/f/workspace/theislab/mubind/data/organoids_treutlein_dataset/RNA_ATAC_metacells_sce_peaks_obs2000_var4000.h5ad',
                       'gbm': '/mnt/f/workspace/theislab/mubind/data/gbm_multiome_scdori/23_01_23_atac_compressed_n1000.h5ad',
                       'noack_2022': '/mnt/f/workspace/theislab/mubind/data/noack_2022/*/merged_scATAC_integrated_cicero_faye_chong_obs*_var*.h5ad',
                       'pancreatic_endocrinogenesis': '/mnt/f/workspace/theislab/mubind/data/pancreatic_endocrinogenesis/*/pancreas_multiome_2022_processed_atac_obs*_var*.h5ad',
                       'pbmc': '/mnt/f/workspace/theislab/mubind/data/pbmc/*/processed_laura_2023_obs*_var*.h5ad'}
    return path_by_dataset

def get_auroc(model, dataloader):
    import mubind as mb
    pred = mb.tl.kmer_enrichment(model, dataloader)
    pred = pred[[c for c in pred if c.startswith('p')]]            
    y = dataloader.dataset.rounds
    y[y > 0] = 1 
    return roc_auc_score(y.flatten(), pred.values.flatten())

def get_auprc(model, dataloader):
    import mubind as mb
    pred = mb.tl.kmer_enrichment(model, dataloader)
    pred = pred[[c for c in pred if c.startswith('p')]]            
    y = dataloader.dataset.rounds
    y[y > 0] = 1 
    return average_precision_score(y.flatten(), pred.values.flatten())

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
