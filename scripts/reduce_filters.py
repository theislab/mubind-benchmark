import mubind as mb
import numpy as np
import pandas as pd
import torch
import itertools
import argparse
import os
from os.path import exists
import bindome as bd
import torch.optim as topti
import torch.utils.data as tdata
import matplotlib.pyplot as plt
import pickle
import glob

# Use a GPU if available, as it should be faster.
# device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
# print("Using device: " + str(device))
if __name__ == '__main__':
    """
    read fastq files, prepare input files for modeling
    """
    parser = argparse.ArgumentParser(
        description='Precompute diffusion connectivities for knn data integration methods.')

    parser.add_argument('-i', '--input', required=True, help='the path to fit_model metricss=')
    parser.add_argument('-o', '--output', required=True, help='path to write models to output')
    args = parser.parse_args()

    models_directory = os.path.join(os.path.abspath(os.path.join(args.input, os.pardir)), 'models')
    # print('reading models from', models_directory)
    models_list = glob.glob(models_directory + '/*.pkl')

    # make output directory
    outdir = os.path.abspath(os.path.join(args.output, os.pardir))
    if not os.path.exists(outdir):
        os.mkdir(outdir)


    binding_modes = mb.tl.binding_modes(models_directory + '/*')
    if len(binding_modes) != 0:
        reduced_groups = mb.tl.reduce_filters(binding_modes)

    else:
        reduced_groups = []

    pickle.dump(reduced_groups, open(args.output, 'wb'))
