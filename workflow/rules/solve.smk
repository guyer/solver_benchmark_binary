rule solve:
    output:
        "results/{id}/solver.log"
    input:
        config="results/{id}/config.json",
        benchmark="results/{id}/benchmark.py"
    params:
        config=lambda w, input: read_config(input.config),
        output=lambda w, output: os.path.dirname(output[0])
    conda:
        get_conda_environment_from_id
    log:
        "results/{id}/solver.stderr"
    shell:
        r"""
        FIPY_SOLVERS={params.config[suite]} \
            python {input.benchmark:q} \
            --solver={params.config[solver]} \
            --preconditioner={params.config[preconditioner]} \
            --numberOfElements={params.config[size]} \
            --output={params.output:q} \
            --restart resources/t=300.0.npz \
            --totaltime=301 \
            --checkpoint_interval=1. \
            2> {log:q} \
            || touch {output:q}
        """

rule copy_script:
    output:
        benchmark="results/{id}/benchmark.py"
    input:
        benchmark=get_benchmark
    conda:
        "../envs/snakemake.yml"
    log:
        "logs/copy_script_{id}.log"
    shell:
        "cp {input.benchmark:q} {output.benchmark:q} 2> {log:q}"

rule ipynb2py:
    input:
        "workflow/notebooks/{notebook}.py.ipynb"
    output:
        temp("workflow/scripts/{notebook}.py")
    conda:
        "../envs/snakemake.yml"
    log:
        stdout="workflow/scripts/{notebook}.stdout",
        stderr="workflow/scripts/{notebook}.stderr"
    shell:
        r"""
        jupyter nbconvert {input:q} --to python \
            --output-dir=workflow/scripts/ \
            --output {wildcards.notebook:q}.py \
            > {log.stdout:q} 2> {log.stderr:q}
        """

rule make_config:
    output:
        "results/{id}/config.json"
    input:
        "config/all_permutations.csv"
    log:
        "logs/make_config_{id}.log"
    run:
        extract_config_by_id(wildcards, output[0], log[0])
