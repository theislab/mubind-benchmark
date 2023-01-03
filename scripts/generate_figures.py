import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

if __name__ == '__main__':
    """
    prepares plots by collecting the results
    """
    parser = argparse.ArgumentParser(description='Prepare plots.')
    parser.add_argument('-o', '--output', required=True, help='output directory for plots')
    parser.add_argument('-i', '--input', required=True, help='input metrics files', nargs="*")
    args = parser.parse_args()

    # create output directory
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    print('Starting plotting...')
    tf_names = []
    results = dict()
    for f in args.input:
        # iterates over metrics file for genes
        tf_name = f.split('/')[-2]
        tf_names.append(tf_name)
        metrics = pd.read_csv(f)
        results[tf_name] = metrics

    r2_y_vars = ['r2_counts', 'r2_foldchange', 'r2_enr', 'r2_fc', 'pearson_foldchange']
    params_x = ['n_epochs', 'n_kernels', 'learning_rate', 'batch_size']

    # PLOT #1. one file for each tf. r2 values wrt all hyperparams
    for tf in tf_names:
        plot_name = f'{args.output}/{tf}_full.png'
        # print(tf, results[tf])
        if len(results[tf]) == 0:
            continue
        g = sns.PairGrid(results[tf],
                        y_vars=r2_y_vars, x_vars=params_x,
                        despine=True)
        g.map(sns.boxplot)
        g.savefig(plot_name)
        print(f'Saved {plot_name}.')


    # PLOT #2. single file. r2 values for libraries
    merged_results = pd.concat(results.values())
    if merged_results.shape[0] != 0:
        plot_name = f'{args.output}/library_r2.png'

        g = sns.PairGrid(merged_results,
                        y_vars=r2_y_vars, x_vars=['library'],
                        despine=True, aspect=len(merged_results.groupby('library'))/2)
        plt.xticks(rotation=90)
        g.map(sns.boxplot)
        g.savefig(plot_name)
        print(f'Saved {plot_name}.')



