import pathlib
from pipeline_config import *
from os.path import join

# configfile: "config_test.yaml"
cfg = ParsedConfig(config, gene_names='gene_names_all.txt')
# cfg = ParsedConfig(config, gene_names='gene_names_cardiac_complexes.txt')
# cfg = ParsedConfig(config, gene_names='gene_names_test.txt')


rule all:
    input:
        files = expand(str(cfg.ROOT) + '/{gene_names}/{experiment}/refine_model/metrics.tsv',
                       gene_names=cfg.GENE_NAMES, experiment=cfg.EXPERIMENT),
        script = 'scripts/generate_figures.py',
    output:
        # results = 'figures/{gene_names}_full.png',
        results = directory(str(cfg.ROOT) + '/results'),
    params:
        cmd = f'conda run -n {cfg.py_env} python',
    shell: "{params.cmd} {input.script} -i {input.files} -o {output.results}"


rule merge_results:
    input:
        files = str(cfg.ROOT) + '/{gene_names}/{experiment}/refine_model/metrics.tsv',
        script = "scripts/merge_metrics.py",
    output:
        results = str(cfg.ROOT) + '/results.tsv.gz'
    params:
        cmd = f"conda run -n {cfg.py_env} python",
    shell: "{params.cmd} {input.script} -i {input.files} -o {output.results}"

rule reduce_filters:
    input:
        metrics = str(cfg.ROOT) + '/{gene_names}/{experiment}/fit_model/metrics.tsv',
        script = "scripts/reduce_filters.py",
    output:
        results = str(cfg.ROOT) + '/{gene_names}/{experiment}/reduce_filters/reduced.pkl',
    params:
        cmd = f"conda run -n {cfg.py_env} python",
    shell: "{params.cmd} {input.script} -i {input.metrics} -o {output.results}"

rule fit_model:
    input:
        script = "scripts/fit_model.py",
        queries = str(cfg.ROOT) + '/{gene_names}/{experiment}/prepare/queries.tsv'
    output:
        metrics = str(cfg.ROOT) + '/{gene_names}/{experiment}/fit_model/metrics.tsv',
    message: "Merge all metrics"
    # log:
    #     str(cfg.ROOT) + '/logs/{gene_names}/{experiment}/fit_model.out.txt'
    params:
        cmd = f"conda run -n {cfg.py_env} python",
        model = str(cfg.ROOT) + '/{gene_names}/{experiment}/fit_model/models',
        experiment = '{experiment}',
        use_mono = expand('{n_epochs}', n_epochs=cfg.HYPERPARAM['use_mono']),
        use_dinuc = expand('{n_epochs}', n_epochs=cfg.HYPERPARAM['use_dinuc']),
        dinuc_mode = expand('{dinuc_mode}', dinuc_mode=cfg.HYPERPARAM['dinuc_mode']),
        n_epochs = expand('{n_epochs}', n_epochs=cfg.HYPERPARAM['n_epochs']),
        batch_sizes = expand('{batch_sizes}', batch_sizes=cfg.HYPERPARAM['batch_sizes']),
        learning_rates = expand('{learning_rates}', learning_rates=cfg.HYPERPARAM['learning_rates']),
        n_kernels = expand('{n_kernels}', n_kernels=cfg.HYPERPARAM['n_kernels']),
        width= expand('{n_kernels}',n_kernels=cfg.HYPERPARAM['width']),
        seq_bias= expand('{seq_bias}',seq_bias=cfg.HYPERPARAM['seq_bias']),

    shell:
        """
        {params.cmd} {input.script} -i {input.queries} --experiment {params.experiment} --out_model {params.model} \
        --out_tsv {output.metrics} --use_mono {params.use_mono} --use_dinuc {params.use_dinuc} \
        --dinuc_mode {params.dinuc_mode} --n_epochs {params.n_epochs} --batch_sizes {params.batch_sizes} \
        --learning_rates {params.learning_rates} --n_kernels {params.n_kernels} --width {params.width} \
        --seq_bias {params.seq_bias} # 1> {log}
        """

rule refine_model:
    input:
        script = "scripts/refine_model.py",
        queries = str(cfg.ROOT) + '/{gene_names}/{experiment}/prepare/queries.tsv',
        reduced_filters = str(cfg.ROOT) + '/{gene_names}/{experiment}/reduce_filters/reduced.pkl'
    output:
        metrics = str(cfg.ROOT) + '/{gene_names}/{experiment}/refine_model/metrics.tsv',
    message: "Merge all metrics"
    # log:
    #     str(cfg.ROOT) + '/logs/{gene_names}/{experiment}/fit_model.out.txt'
    params:
        cmd = f"conda run -n {cfg.py_env} python",
        model = str(cfg.ROOT) + '/{gene_names}/{experiment}/models',
        experiment = '{experiment}',
        use_mono = expand('{n_epochs}', n_epochs=cfg.HYPERPARAM['use_mono']),
        use_dinuc = expand('{n_epochs}', n_epochs=cfg.HYPERPARAM['use_dinuc']),
        dinuc_mode = expand('{dinuc_mode}', dinuc_mode=cfg.HYPERPARAM['dinuc_mode']),
        n_epochs = expand('{n_epochs}', n_epochs=cfg.HYPERPARAM['n_epochs']),
        batch_sizes = expand('{batch_sizes}', batch_sizes=cfg.HYPERPARAM['batch_sizes']),
        learning_rates = expand('{learning_rates}', learning_rates=cfg.HYPERPARAM['learning_rates']),
        n_kernels = expand('{n_kernels}', n_kernels=cfg.HYPERPARAM['n_kernels']),
    shell:
        """
        {params.cmd} {input.script} -i {input.queries} --experiment {params.experiment} --out_model {params.model} \
        --out_tsv {output.metrics} --use_mono {params.use_mono} --use_dinuc {params.use_dinuc} \
        --dinuc_mode {params.dinuc_mode} --n_epochs {params.n_epochs} --batch_sizes {params.batch_sizes} \
        --learning_rates {params.learning_rates} --filters {input.reduced_filters} # 1> {log}
        """


rule data_prepare:
    input:
        script = "scripts/data_prepare.py",
    output:
        queries = str(cfg.ROOT) + '/{gene_names}/{experiment}/prepare/queries.tsv'
    message:
        """
        Preparing adata
        wildcards: {wildcards}
        parameters: {params}
        output: {output}
        """
    # log:
    #     str(cfg.ROOT) + '/logs/{gene_names}/{experiment}/data_prepare.out.txt'
    params:
        cmd       = f"conda run -n {cfg.py_env} python",
        # tf_name = expand('{gene_names}', gene_names=cfg.GENE_NAMES),
        # tf_name = lambda wildcards: 'wildcards.gene_names',
        experiment = '{experiment}',
        tf_name = '{gene_names}'.split('//')[0],
        annot = f"{cfg.ANNOTATIONS}",
        n_sample = expand('{n_sample}', n_sample=cfg.HYPERPARAM['n_sample']),
    shell:
        """
        {params.cmd} {input.script} --annot {params.annot} --experiment {params.experiment} --tf_name {params.tf_name} -o {output.queries} --n_sample {params.n_sample} # 1> {log}
        """

