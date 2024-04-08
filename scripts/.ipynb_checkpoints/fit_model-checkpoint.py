import mubind as mb
import numpy as np
import pandas as pd
import torch
import itertools
from os.path import exists
import bindome as bd
import torch.optim as topti
import torch.utils.data as tdata
import matplotlib.pyplot as plt
import pickle

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
    parser.add_argument('--out_model', required=True, help='directory to save learned model')
    parser.add_argument('--out_tsv', required=True, help='output path for metrics')
    parser.add_argument('--early_stopping', default=100, help='# epochs for early stopping', type=int)
    parser.add_argument('--is_count_data', default=True)

    parser.add_argument('--outdense', default=False)

    parser.add_argument('--n_epochs', nargs='+', default=[50], help='# of epochs for training', type=int)
    parser.add_argument('--batch_sizes', nargs='+', default=[32], help='batch sizes for training', type=int)
    parser.add_argument('--learning_rates', nargs='+', default=[0.001], help='learning rates for training', type=float)
    parser.add_argument('--n_kernels', nargs='+', default=[32], help='# of kernels for training', type=int)

    # Use a GPU if available, as it should be faster.
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print("Using device: " + str(device))

    args = parser.parse_args()

    queries = pd.read_csv(args.queries, sep='\t', index_col=0)

    out_model = args.out_model
    if not os.path.exists(out_model):
        os.mkdir(out_model)

    metrics = []
    for (n_epochs, lr, batch_size, n_kernels) in itertools.product(args.n_epochs, args.learning_rates, args.batch_sizes,
                                                                   args.n_kernels):
        hyperparam_id = f'N_{n_epochs}_L_{lr}_B_{batch_size}_K_{n_kernels}'
        print(f'Training with n_epochs={n_epochs}, lr={lr}, batch_size={batch_size}')

        for ri, r in queries.iterrows():
            counts_path = r['counts_path']

            model_path = f'{out_model}/{hyperparam_id}_{os.path.basename(counts_path).replace(".tsv.gz", ".h5")}'
            pkl_path = model_path.replace('.h5', '.pkl')

            data = pd.read_csv(counts_path, sep=',' if 'simulated' in counts_path else '\t', index_col=0)

            if 'simulated' not in counts_path:
                n_rounds = len(data.columns) - 3
                n_batches = len(set(data.batch))
            else:
                n_rounds = 1
                n_batches = 1
            print('# rounds', n_rounds)
            print('# batches', n_batches)

            print(data.shape)
            # data = data.head(1000)
            # print(data.shape)

            labels = list(data.columns[:n_rounds + 1])
            print('labels', labels)
            dataset = mb.datasets.SelexDataset(data, n_rounds=n_rounds, labels=labels)
            train = tdata.DataLoader(dataset=dataset, batch_size=batch_size, shuffle=True)
            # train_test = tdata.DataLoader(dataset=dataset, batch_size=1, shuffle=False)
            # train_test = tdata.DataLoader(dataset=dataset, batch_size=1, shuffle=False)
            # train_test = tdata.DataLoader(dataset=dataset, batch_size=1, shuffle=False)
            ### steps to train model

            print(exists(model_path), model_path)
            if not exists(model_path):
                print('training starts...')
                model, best_loss = mb.tl.train_iterative(train, device, show_logo=False, log_each=50,
                                                         num_epochs=n_epochs, n_kernels=n_kernels,
                                                         early_stopping=args.early_stopping, lr=lr)
                torch.save(model.state_dict(), model_path)
                pickle.dump(model, open(pkl_path, 'wb'))
            else:
                model = pickle.load(open(pkl_path, 'rb'))
                model.load_state_dict(torch.load(model_path))

            r2 = mb.pl.kmer_enrichment(model, train, k=8, show=False)
            print("R^2:", r2)

            if args.outdense:
                # print r2 for all epochs
                for idx, val in enumerate(model.r2_history):
                    metrics.append(list(r.values[:-1]) + [batch_size, lr, idx, n_kernels, -1, val, -1])

            best_r2 = max(model.r2_history)
            metrics.append(list(r.values[:-1]) + [batch_size, lr, n_epochs, n_kernels, model.best_loss, best_r2, model.total_time])

        metrics = pd.DataFrame(metrics, columns=list(queries.columns[:-1]) + ['batch_size', 'learning_rate', 'n_epochs',
                                                                              'n_kernels', 'best_loss', 'best_r2',
                                                                              'running_time'])
    metrics.to_csv(args.out_tsv)
