# mubind  (snakemake)

[![Tests][badge-tests]][link-tests]

[badge-tests]: https://img.shields.io/github/workflow/status/ilibarra/mubind_pipeline/Test/main
[link-tests]: https://github.com/theislab/mubind_pipeline/actions/workflows/test.yml

Snakemake pipeline related to mubind main repository.

### Test
To test that the pipeline is running properly, as a small test, please execute `config_test.yaml`
```bash
snakemake -c 1 --printshellcmds --configfile config_test.yaml --nolock
```
The parameter `-c` indicates that only one core is used.

### Full pipeline

#### Run until the data preparation step
```bash
snakemake -c 5 --printshellcmds --configfile config.yaml --until data_prepare --nolock
```

#### Fit model preparation step
```bash
snakemake -c 5 --printshellcmds --configfile config.yaml --until fit_model --nolock
```
This is useful once the datasets are prepared, to only define GPU nodes if using `--cluster`
