rule solve:
    output:
        "results/fipy~{rev}/suite~{suite}/{id}/solver.log"
    input:
        env="results/fipy~{rev}/suite~{suite}/environment.yml",
        config=config["simulations"],
        benchmark=get_benchmark,
    params:
        config=get_config_by_id,
    conda:
        "../../results/fipy~{rev}/suite~{suite}/environment.yml"
    log:
        "logs/fipy~{rev}/suite~{suite}/{id}/notebooks/benchmark.log"
    benchmark:
        "benchmarks/fipy~{rev}/suite~{suite}/benchmark-{id}.tsv"
    shell:
        r"""
        FIPY_SOLVERS={wildcards.suite} \
            python {input.benchmark:q} \
            --solver={params.config[solver]} \
            --preconditioner={params.config[preconditioner]} \
            --numberOfElements={params.config[size]} \
            --restart={params.config[restart]} \
            --totaltime={params.config[totaltime]} \
            --checkpoint_interval={params.config[checkpoint_interval]} \
            --solve_log={output:q}
            2> {log:q} \
            || touch {output:q}
        """

rule ipynb2py:
    output:
        "workflow/scripts/{notebook}.py"
    input:
        "workflow/notebooks/{notebook}.py.ipynb"
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
