#!/usr/bin/env python
# coding: utf-8
import mubind as mb
import numpy as np
import pandas as pd
import bindome as bd
import sys
import os
from pathlib import Path

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
    parser.add_argument('--n_sample', default=None, type=int)
    parser.add_argument('-o', '--output', required=True, help='output directory for counts and queries metadata file')

    args = parser.parse_args()
    queries = []
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
            model_by_k = {}
            for tf in tf_queries:  # set(data['tf_name']):
                if 'ZERO' in tf:
                    continue
                print(tf)

                for library, grp in data.groupby('library'):
                    try: 
                        data_sel_tf = grp[(grp['tf_name'] == tf)]  # & (grp['cycle'] == '1')]
                    except:
                        data_sel_tf = grp[(grp['tf.name'] == tf)]
                    if data_sel_tf.shape[0] == 0:
                        continue

                    print('loading', tf, ':', library)

                    reads_tf = mb.bindome.datasets.SELEX.load_read_counts(tf, data=data_sel_tf)
                    try:
                        data_sel_zero = grp[(grp['cycle'] == 0) & grp['library'].isin(set(grp[grp['tf_name'] == tf][
                                                                                          'library']))]  # & grp['accession'].isin(set(grp[grp['tf_name'] == tf]['accession']))]
                    except:
                        data_sel_zero = grp[(grp['cycle'] == 0) & grp['library'].isin(set(grp[grp['tf.name'] == tf][
                                                                                          'library']))]  # & 
                    reads_zero = mb.bindome.datasets.SELEX.load_read_counts(data=data_sel_zero, library=library)

                    print('# zero files found', data_sel_zero.shape)
                    print(reads_tf.keys())
                    print(reads_zero.keys())

                    for k_r0 in reads_zero:
                        print(k_r0)

                        k_model = tf + '-' + k_r0 + '-' + library

                        n_rounds = 1

                        next_data = reads_zero[k_r0].copy()
                        new_cols = ['seq', k_r0]
                        for k_tf in reads_tf:
                            next_data = next_data.merge(reads_tf[k_tf], on='seq', how='outer').fillna(
                                0)  # .astype(int)
                            new_cols.append(k_tf)
                        # next_data = reads_zero[k_r0].merge(reads_tf[k_tf], on='seq', how='outer').fillna(0) # .astype(int)
                        # new_cols = ['seq', k_r0, k_tf]
                        next_data.columns = new_cols
                        for i, k in enumerate(new_cols[1:]):
                            next_data[k] = next_data[k].astype(int)
                            # next_data[i] = next_data[k].astype(int)

                        # next_data = next_data.head(10000)
                        if args.n_sample is not None and args.n_sample != -1:
                            next_data = next_data.sample(n=args.n_sample)

                        # print(next_data.shape)
                        next_data = next_data.set_index('seq')
                        # print(next_data.head())

                        # not needed for the current model, because the enrichment is not predicted
                        # next_data = mb.tl.calculate_enrichment(next_data, cols=next_data.columns[1:])

                        # assign batch and data type
                        next_data['batch'] = 1
                        next_data['is_count_data'] = 1

                        next_outpath = str(queries_directory) + '/' + k_model + '.tsv.gz'
                        next_data.to_csv(next_outpath, sep='\t')

                        queries.append([tf_query, k_r0, library, next_outpath, next_data.shape[0]])

    queries = pd.DataFrame(queries, columns=['tf_name', 'r0', 'library', 'counts_path', 'n_sample'])
    queries.to_csv(queries_tsv_outpath, sep='\t')
    sys.exit()

