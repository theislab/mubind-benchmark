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

    # 1. one file for each tf, plot r2 values wrt all hyperparams 
    for tf in tf_names:
        plot_name = f'{args.output}/{tf}_full.png'
        g = sns.PairGrid(results[tf], y_vars=['r2_counts', 'r2_foldchange', 'r2_enr', 'r2_fc', 'pearson_foldchange'],
                        x_vars=['n_epochs', 'n_kernels', 'learning_rate', 'batch_size'],
                        despine=True)
        g.map(sns.boxplot)
        g.savefig(plot_name)
        print(f'Saved {plot_name}.')

    
    

