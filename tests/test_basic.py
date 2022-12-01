import numpy as np
import pandas as pd
import torch
import torch.optim as topti
import torch.utils.data as tdata
import os

def test_configtest_exists():
    output_state = os.system('snakemake --configfile config_test.yaml -c1 --dry-run')
    assert output_state == 0

def test_configtest_executes():
    output_state = os.system('snakemake --configfile config_test.yaml -c1')
    assert output_state == 0

