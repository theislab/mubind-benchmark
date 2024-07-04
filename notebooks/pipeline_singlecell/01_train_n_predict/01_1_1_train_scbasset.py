#!/usr/bin/env python
# coding: utf-8

# ### List datasets that will be used for training with scBasset

# In[1]:


# %load_ext autoreload
# %autoreload 2


# In[2]:


# maybe anndata and others required
# !pip install anndata h5py


# In[3]:


import scbasset
import anndata


# In[4]:
# cd '~/workspace/theislab/mubind-pipeline/notebooks/pipeline/01_train_n_predict'


# In[5]:


import utils
import glob


# ### Load the datasets used in the core benchmark

# In[6]:

path_by_dataset = utils.get_datasets()
path_by_dataset


# In[7]:


# cd ~/workspace/scbasset/bin


# ## Run python using the scbasset environment. If not found, please set up manually

# In[8]:


pythonbin = '/home/ilibarra/.conda/envs/scbasset/bin/python' # ''


# In[9]:


overwrite_prepare = True
overwrite_train = True
n_epochs = 10


# In[10]:


# adata = anndata.read_h5ad(ad_path)
# print(adata.shape)
# adata.var['chr'] = adata.var_names.str.split('-').str[0]
# adata.var['chr'] = np.where(~adata.var['chr'].str.contains('chr'), 'chr', '') + adata.var['chr']
# adata.var['start'] = adata.var_names.str.split('-').str[1]
# adata.var['end'] = adata.var_names.str.split('-').str[2]
# adata.write(ad_path, compression='lzf')
# genome_fasta = '/mnt/c/Users/ignacio.ibarra/Dropbox/annotations/%s/genome/%s.fa' % (species, species)
# os.path.exists(genome_fasta)
# print(genome_fasta)


# In[11]:


for k in path_by_dataset:
    print(k, path_by_dataset[k])


# In[19]:


import os
import anndata
import os
import numpy as np

genomes_path = '/mnt/c/Users/IgnacioIbarra/Dropbox/annotations/%s/genome/%s.fa'

for k in path_by_dataset:
    try:
        print(k)
        # if not 'noack_2022' in k:
        #     continue
        for ad_path in glob.glob(path_by_dataset[k]):
            if not '_obs1000_' in ad_path and not '_obs500_' in ad_path: #  or not '/random/' in ad_path:
                continue
            print(os.path.exists(ad_path), ad_path)
            # print('\nnext path')
            # print(k, ad_path)

            # continue

            n_cells = ad_path.split('_')[-2]
            sampling_method = ad_path.split('/')[-2]

            print(k, sampling_method, ad_path)

            species = 'hg38' if (k != 'noack_2022' and k != 'pancreatic_endocrinogenesis') else 'mm10'

            filename = os.path.basename(ad_path)
            input_dir = os.path.join(ad_path.replace(filename, 'scbasset_input'), n_cells)

            for use_poisson in [1, 0]:
                out_dir = os.path.join(ad_path.replace(filename,
                                                    'scbasset_output'),
                                                    n_cells, 'poisson' if use_poisson else 'bce')
                
                print('in', input_dir)
                print('out', out_dir)

                if not os.path.exists(input_dir):
                    os.makedirs(input_dir)
                if not os.path.exists(out_dir):
                    os.makedirs(out_dir)
                
                print(species)
                train_path = os.path.join(out_dir, 'ad.h5ad')
                train_seqs_path = os.path.join(out_dir, 'train_seqs.h5')
                if not os.path.exists(train_path) or not os.path.exists(train_seqs_path) or overwrite_prepare:
                    adata = anndata.read_h5ad(ad_path)
                    print(adata.shape)

                    # assert False
                    adata.var['chr'] = adata.var_names.str.split('-').str[0]
                    adata.var['chr'] = np.where(~adata.var['chr'].str.contains('chr'), 'chr', '') + adata.var['chr']
                    adata.var['start'] = adata.var_names.str.split('-').str[1]
                    adata.var['end'] = adata.var_names.str.split('-').str[2]
                    adata.write(ad_path, compression='lzf')
                    genome_fasta = genomes_path % (species, species)
                    os.path.exists(genome_fasta)
                    print(genome_fasta)
                    # preprocess fasta
                    print('preprocess fasta...')
                    get_ipython().system('echo $pythonbin scbasset_preprocess.py --ad_file $ad_path --input_fasta $genome_fasta --out_path $input_dir --batch 50')
                    get_ipython().system('$pythonbin scbasset_preprocess.py --ad_file $ad_path --input_fasta $genome_fasta --out_path $input_dir --batch 50')
                    print('')
                else:
                    print('skip prepare (already done...)')
                
                best_model_path = os.path.join(out_dir, 'best_model.h5')
                if not os.path.exists(best_model_path) or overwrite_train:
                    print('')
                    print('run scbasset...')
                    get_ipython().system('echo $pythonbin scbasset_train.py --input_folder $input_dir --epochs $n_epochs --out_path $out_dir')
                    get_ipython().system('$pythonbin scbasset_train.py --input_folder $input_dir --epochs $n_epochs --out_path $out_dir --use_poisson $use_poisson')
                    print('')

                else:
                    print('skip train (already done...)')
                print('')

                assert False

    except Exception:
        pass

    print('')

