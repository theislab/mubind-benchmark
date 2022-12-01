import numpy as np
import pandas as pd
import torch
import torch.optim as topti
import torch.utils.data as tdata
import os

def test_simdata_train():

    # print('here...l')
    os.system('snakemake --configfile config.yaml --dry-run')

