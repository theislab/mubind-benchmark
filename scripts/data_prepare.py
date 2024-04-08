#!/usr/bin/env python
# coding: utf-8
import mubind as mb
import numpy as np
import pandas as pd
import bindome as bd
import sys
import os
from pathlib import Path
import scipy

if __name__ == '__main__':
    """
    read fastq files, prepare input files for modeling
    """

    import argparse
    import os

    parser = argparse.ArgumentParser(
        description='Prepare data related to counts/intensities for modeling with mubind.')

    parser.add_argument('--annot', help='annotations directory')
    parser.add_argument('--experiment', help='the code of the experiment type e.g. SELEX/PBM', required=True)
    parser.add_argument('--tf_name', required=True)
    parser.add_argument('--n_sample', nargs='+', default=[None], type=int)
    parser.add_argument('-o', '--output', required=True, help='output directory for counts and queries metadata file')

    args = parser.parse_args()
    queries = []
    total_counts = []
    queries_tsv_outpath = args.output
    queries_directory = Path(queries_tsv_outpath).parent.absolute()

    if args.tf_name == 'SELEX_SIM':
        counts_path = f'{args.annot}/selex/simulated'
        for i in range(20):
            queries.append(['SELEX_SIM', 'k_r0', 'library',
                            os.path.abspath(f'{counts_path}/{i}_counts.csv.gz'),
                            args.n_sample])
    else:
        bd.constants.ANNOTATIONS_DIRECTORY = args.annot
        if not os.path.exists(queries_directory):
            print(queries_directory)
            os.makedirs(queries_directory)

        if args.experiment == 'PBM':
            print('loading PBM...')
            ad = bd.datasets.PBM.uniprobe()
            if ad is not None: # check that h5ad is available in the annotations directory
                ad_sel = ad[ad.obs.index.str.lower().str.contains(args.tf_name.lower()), :]
                for next_pbm in ad_sel.obs_names:
                    next_data = ad_sel.to_df()
                    next_data = next_data[next_data.index == next_pbm]
                    next_data = next_data.T
                    next_data = next_data[~np.isnan(next_data).any(axis=1)]
                    next_outpath = str(queries_directory) + '/' + next_pbm + '.tsv.gz'
                    next_data.index.name = 'seq'
                    # next_data = next_data.head(10000)
                    if args.n_sample is not None and args.n_sample != -1:
                        next_data = next_data.sample(n=args.n_sample)

                    # assign batch and data type
                    next_data['batch'] = 1
                    next_data['is_count_data'] = 0

                    next_data.to_csv(next_outpath, sep='\t')
                    queries.append([args.tf_name, np.nan, next_pbm, next_outpath, next_data.shape[0]])

        if args.experiment == 'SELEX':

            data = bd.bindome.datasets.SELEX.get_data()
            tf_query = args.tf_name
            tf_queries = {tf_query}
            # assert False
            model_by_k = {}
            for tf in tf_queries:  # set(data['tf_name']):
                if 'ZERO' in tf:
                    continue
                print(tf)

                for library, grp in data.groupby('library'):
                    # print(grp)

                    try:
                        data_sel_tf = grp[(grp['tf_name'] == tf)]  # & (grp['cycle'] == '1')]
                    except:
                        data_sel_tf = grp[(grp['tf.name'] == tf)]
                    if data_sel_tf.shape[0] == 0:
                        continue

                    print('loading', tf, ':', library)

                    # print(grp)
                    # assert False

                    reads_tf = mb.bindome.datasets.SELEX.load_read_counts(tf, data=data_sel_tf)
                    try:
                        data_sel_zero = grp[(grp['cycle'] == 0) & grp['library'].isin(set(grp[grp['tf_name'] == tf][
                                                                                          'library']))]  # & grp['accession'].isin(set(grp[grp['tf_name'] == tf]['accession']))]
                    except:
                        data_sel_zero = grp[(grp['cycle'] == 0) & grp['library'].isin(set(grp[grp['tf.name'] == tf][
                                                                                          'library']))]
                    reads_zero = mb.bindome.datasets.SELEX.load_read_counts(data=data_sel_zero, library=library)

                    print('# zero files found', data_sel_zero.shape)
                    print(reads_tf.keys())
                    print(reads_zero.keys())

                    for k_r0 in reads_zero:
                        print(k_r0)

                        k_model = tf + '-' + k_r0 + '-' + library

                        n_rounds = 1

                        next_r0 = reads_zero[k_r0].copy()
                        new_cols = ['seq', k_r0]
                        next_r0.columns = ['seq', k_r0]
                        print(next_r0.head())
                        next_data = None
                        print('')
                        for k_tf in reads_tf:
                            print(k_tf)
                            # print(next_r0.head())
                            # print(reads_tf[k_tf].head())
                            next_data = (next_r0 if next_data is None else next_data)

                            next_data_tf = reads_tf[k_tf].copy()
                            next_data_tf.columns = ['seq', k_tf]
                            print(next_data.columns)
                            next_data = next_data.merge(next_data_tf,
                                                        on='seq',
                                                        how='outer').fillna(0)
                            
                            print(next_data.head())
                            print('')
                            new_cols.append(k_tf)
                        next_data.columns = new_cols

                        for i, k in enumerate(new_cols[1:]): # convert all columns into int
                            next_data[k] = next_data[k].astype(int)

                        # print('next data')
                        # print(next_data.head())
                        for n_sample in args.n_sample:
                            # next_data = next_data.head(10000)
                            if args.n_sample is not None and args.n_sample != -1:
                                select_random = True
                                if not select_random and next_data.shape[1] > 2: # sort by descending enrichment (seq + [0, 1, 2, ...])
                                    fg = (next_data[next_data.columns[2:]].max(axis=1) + 1)
                                    # print(next_data[next_data.columns[1]])
                                    bg = next_data[next_data.columns[1]] + 1

                                    enrichment_descending = (fg / bg).sort_values(ascending=False)
                                    next_data_sample = next_data.reindex(enrichment_descending.index[:n_sample]).copy()
                                else:
                                    print('sampling randomly...')
                                    next_data_sample = next_data.sample(n=min(n_sample, next_data.shape[0]))

                                print('# total/sampled reads:', next_data.shape[0], next_data_sample.shape[0])
                            else:
                                next_data_sample = next_data.copy()

                            # print(next_data.shape)
                            next_data_sample = next_data_sample.set_index('seq')
                            # print(next_data.head())

                            # remove indexes with N as those are not encoded with mb.tl.string2bin
                            print(next_data_sample.head())
                            if next_data_sample.index.astype(str).str.contains('N').any():
                                print('removing Ns...before/after')
                                print(next_data_sample.shape[0])
                                next_data_sample = next_data_sample[~next_data_sample.index.str.contains('N')]
                                print(next_data_sample.shape[0])

                            # assign batch and data type
                            next_data_sample['batch'] = 1
                            next_data_sample['is_count_data'] = 1

                            next_outpath = str(queries_directory) + '/' + k_model + '_%s.tsv.gz' % n_sample

                            next_outpath_sparse = next_outpath.replace(".tsv.gz", '_sparse.npz')
                            next_outpath_rownames = next_outpath.replace(".tsv.gz", '_sparse_rownames_int.npz')
                            next_outpath_colnames = next_outpath.replace(".tsv.gz", '_colnames.npz')

                            X = scipy.sparse.csr_matrix(next_data_sample.to_numpy())
                            colnames = np.array(next_data_sample.columns)
                            rownames = np.array(next_data_sample.index)
                            rownames_int = pd.Series(rownames).apply(mb.tl.string2bin)
                            scipy.sparse.save_npz(next_outpath_sparse, X)
                            np.savez(next_outpath_rownames, rownames_int)
                            np.savez(next_outpath_colnames, colnames)

                            print('next counts saved at')
                            print(next_outpath_sparse)
                            print(next_outpath_rownames)
                            print(next_outpath_colnames)
                            # next_data_sample.to_csv(next_outpath, sep='\t')

                            # repeats for each query (e.g. two in ALX1), want to save sep file for each query so should save file here
                            queries.append([tf_query, k_r0, library, next_outpath, n_sample, next_data_sample.shape[0]])
                            total_counts = next_data_sample.sum(axis=0)
                            # total_counts = pd.DataFrame(total_counts, columns=['query', 'total_counts'])
                            total_counts.to_csv('/'.join(queries_tsv_outpath.split('/')[:-1]) +
                                                f'/{k_r0}_%s_counts.tsv' % n_sample, sep='\t')

    queries = pd.DataFrame(queries, columns=['tf_name', 'r0', 'library', 'counts_path', 'n_sample_parm', 'n_sample_obs'])
    queries.to_csv(queries_tsv_outpath, sep='\t')

    print('queries saved at')
    print(queries_tsv_outpath)

