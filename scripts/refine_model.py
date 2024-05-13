import mubind as mb
import numpy as np
import pandas as pd
import torch
import glob
import itertools
from os.path import exists
import bindome as bd
import torch.optim as topti
import torch.utils.data as tdata
import matplotlib.pyplot as plt
import pickle
import sys

# Use a GPU if available, as it should be faster.
# device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
# print("Using device: " + str(device))
if __name__ == '__main__':
    """
    read fastq files, prepare input files for modeling
    """

    import argparse
    import os

    parser = argparse.ArgumentParser(
        description='Precompute diffusion connectivities for knn data integration methods.')

    parser.add_argument('-i', '--queries', help='path to queries.tsv with entries')
    parser.add_argument('--out_model', required=True, help='directory to save refined model')
    parser.add_argument('--out_tsv', required=True, help='output path for metrics')
    parser.add_argument('--early_stopping', default=100, help='# epochs for early stopping', type=int)
    parser.add_argument('--experiment', help='the experiment type we are going to model', required=True)
    parser.add_argument('--is_count_data', default=True)
    parser.add_argument('--outdense', default=False)
    parser.add_argument('--use_mono', help='using mononucleotide features', default=True)
    parser.add_argument('--use_dinuc', help='using dinucleotide features', default=False)
    parser.add_argument('--dinuc_mode', help='using dinucleotide features', default='local')
    parser.add_argument('--n_epochs', nargs='+', default=[50], help='# of epochs for training', type=int)
    parser.add_argument('--batch_sizes', nargs='+', default=[32], help='batch sizes for training', type=int)
    parser.add_argument('--learning_rates', nargs='+', default=[0.001], help='learning rates for training', type=float)
    parser.add_argument('--n_kernels', nargs='+', default=[1], help='# of kernels for training', type=int)
    # a specific flag to indicate that in this script we can also load and reuse binding modes from an input file.
    # a specific flag to indicate that in this script we can also load and reuse binding modes from an input file.
    parser.add_argument('--filters', help='a path to prior filters, if existing', default=False)

    # Use a GPU if available, as it should be faster.
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print("Using device: " + str(device))

    args = parser.parse_args()

    # make output diretory
    outdir = os.path.abspath(os.path.join(args.out_tsv, os.pardir))
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        os.makedirs(outdir + '/models')

    criterion = mb.tl.PoissonLoss()

    import os
    queries = pd.read_csv(args.queries, sep='\t', index_col=0)

    input_paths = list(queries['counts_path'])

    if len(input_paths) == 0:
        print('# no count paths provided. Write empty dataframe and finish')
        metrics = pd.DataFrame([], columns=list(queries.columns[:-1]) + ['batch_size', 'learning_rate', 'n_epochs',
                                                                              'n_kernels', 'best_loss',
                                                                              'r2_counts', 'r2_foldchange', 'r2_enr',
                                                                              'r2_fc', 'pearson_foldchange',
                                                                              'running_time'])
        metrics.to_csv(args.out_tsv)
        sys.exit()

    import scipy
    df = []
    batch = 0
    batch_by_name = {}
    for p in input_paths:
        # print(p)

        counts_path = p
        sparse_counts_path = counts_path.replace(".tsv.gz", '_sparse.npz')
        rownames_path = counts_path.replace(".tsv.gz", '_sparse_rownames_int.npz')
        colnames_path = counts_path.replace(".tsv.gz", '_colnames.npz')
        X = scipy.sparse.load_npz(sparse_counts_path)
        rownames_int = np.load(rownames_path, allow_pickle=True)
        rownames_int = rownames_int['arr_0']
        print(type(rownames_int))
        row_names = pd.Series(rownames_int).apply(mb.tl.bin2string)
        col_names = np.load(colnames_path, allow_pickle=True)
        col_names = col_names['arr_0']
        df2 = pd.DataFrame(X.toarray(), index=row_names, columns=col_names)
        # df2 = pd.read_csv(p, sep='\t', index_col=0)  # .head(100)

        assert 'batch' in df2
        # print(df2.columns)
        # df2 = df2.sample(100000)
        n_rounds = len(df2.columns) - 2
        df2.columns = list(range(n_rounds)) + ['batch', 'is_count_data']
        df2['batch'] = batch
        batch_by_name[batch] = os.path.basename(p)
        df2['n_rounds'] = n_rounds
        # df2 = mb.pp.sample_rounds(df2, n_rounds, n_sample_per_round)
        # print(df2.shape)
        # print(p, df2.shape, n_rounds)
        batch += 1
        df.append(df2)
    
    df = pd.concat(df)
    df = df[[c for c in df.columns if not c in ['batch', 'is_count_data', 'n_rounds']] + ['batch', 'is_count_data',
                                                                                          'n_rounds']]
    dataset = mb.datasets.SelexDataset(df, n_rounds=df['n_rounds'], labels=list(df.columns[:-3]), store_rev=False)
    train = tdata.DataLoader(dataset=dataset,
                             # batch_size=256,
                             batch_size=512,
                             shuffle=False)

    n_rounds = train.dataset.n_rounds
    n_batches = train.dataset.n_batches
    enr_series = train.dataset.enr_series

    # load the prior filters
    filters = pickle.load(open(args.filters, 'rb'))

    model = mb.models.Mubind(
        datatype="selex",
        kernels=[0] + [m.shape[-1] for m in filters],
        n_rounds=n_rounds,
        init_random=False,
        n_batches=n_batches,
        enr_series=enr_series,
        use_dinuc_full=True,
    ).to(device)
    for i, mono_best in enumerate(filters):
        # print(mono_best.shape, model.binding_modes.conv_mono[i + 1].weight.shape)
        # print(model.binding_modes.conv_mono[i + 1].weight.device)
        new_w = mono_best.reshape([1, 1] + list(mono_best.shape))
        model.binding_modes.conv_mono[i + 1].weight = torch.nn.Parameter(torch.tensor(new_w, dtype=torch.float))
        # print(model.binding_modes.conv_mono[i + 1].weight.device)
    # move the model a final time to the GPU
    model = model.to(device)
    # mb.pl.conv_mono(model, n_rows=len(reduced_groups) + 1, n_cols=2, title=False, xticks=False)


    # fit the model using the data
    # metrics by batch, using combined model
    metrics = []
    metrics = pd.DataFrame(metrics, columns=list(queries.columns[:-1]) + ['batch_size', 'learning_rate', 'n_epochs',
                                                                              'n_kernels', 'best_loss',
                                                                              'r2_counts', 'r2_foldchange', 'r2_enr', 'r2_fc', 'pearson_foldchange',
                                                                              'running_time'])
    
    refined_model_outdir = args.out_model.replace('models', 'refine_model')
    # print(os.path.exists(refined_model_outdir), refined_model_outdir)

    refined_model_path_h5 = os.path.join(refined_model_outdir, 'model.h5')
    refined_model_path_pkl = os.path.join(refined_model_outdir, 'model.pkl')
    motif_img_path = refined_model_path_h5.replace('.h5', '_filters.png')

    if not exists(refined_model_path_h5):
        torch.save(model.state_dict(), refined_model_path_h5)
        pickle.dump(model, open(refined_model_path_pkl, 'wb'))

    metrics.to_csv(args.out_tsv)
