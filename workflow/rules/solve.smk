rule solve:
    output:
        "results/fipy~{rev}/suite~{suite}/{id}/solver.log"
    input:
        config="results/fipy~{rev}/suite~{suite}/{id}/config.json",
    params:
        config=lambda w, input: read_config(input.config),
    conda:
        "../../results/fipy~{rev}/suite~{suite}/environment.yml"
    log:
        notebook="logs/fipy~{rev}/suite~{suite}/{id}/notebooks/benchmark.ipynb"
    notebook:
        "../notebooks/binary_phase_field.py.ipynb"

rule make_config:
    output:
        "results/fipy~{rev}/suite~{suite}/{id}/config.json"
    input:
        "results/fipy~{rev}/suite~{suite}/permutations.csv"
    log:
        "logs/fipy~{rev}/suite~{suite}/{id}/make_config.log"
    run:
        extract_config_by_id(wildcards, input[0], output[0], log[0])
