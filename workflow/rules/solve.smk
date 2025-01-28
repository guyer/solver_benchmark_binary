rule solve:
    output:
        "results/fipy~{rev}/suite~{suite}/{id}/solver.log"
    input:
        env="results/fipy~{rev}/suite~{suite}/environment.yml"
    params:
        config=get_config_by_id,
    conda:
        "../../results/fipy~{rev}/suite~{suite}/environment.yml"
    log:
        notebook="logs/fipy~{rev}/suite~{suite}/{id}/notebooks/benchmark.ipynb"
    benchmark:
        "benchmarks/fipy~{rev}/suite~{suite}/benchmark-{id}.tsv"
    notebook:
        "../notebooks/binary_phase_field.py.ipynb"
