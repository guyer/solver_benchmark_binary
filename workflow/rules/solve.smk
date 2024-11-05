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
        get_conda_environment
    log:
        "results/{id}/solver.stderr"
    shell:
        "FIPY_SOLVERS={params.config[suite]}"
        " python {input.benchmark}"
        " --solver={params.config[solver]}"
        " --preconditioner={params.config[preconditioner]}"
        " --numberOfElements={params.config[size]}"
        " --output={params.output}"
        " --restart resources/t=300.0.npz"
        " --totaltime=301"
        " --checkpoint_interval=1."
        " 2> {log}"
        " || touch {output[0]}"

rule copy_script:
    output:
        benchmark="results/{id}/benchmark.py"
    input:
        benchmark=get_benchmark
    shell:
        "cp {input.benchmark} {output.benchmark}"

rule ipynb2py:
    input:
        "workflow/notebooks/{notebook}.py.ipynb"
    output:
        temp("workflow/scripts/{notebook}.py")
    conda:
        "snakemake"
    log:
        stdout="workflow/scripts/{notebook}.stdout",
        stderr="workflow/scripts/{notebook}.stderr"
    shell:
        "jupyter nbconvert {input} --to python --output-dir=workflow/scripts/ --output {wildcards.notebook}.py > {log.stdout} 2> {log.stderr}"

rule make_config:
    output:
        "results/{id}/config.json"
    input:
        "results/permutations.json"
    run:
        permutations.loc[wildcards.id].to_json(output[0])
