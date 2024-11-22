rule solve:
    output:
        "results/{id}/solver.log"
    input:
        config="results/{id}/config.json",
    params:
        config=lambda w, input: read_config(input.config),
    conda:
        get_conda_environment_from_id
    log:
        notebook="logs/notebooks/benchmark_{id}.ipynb"
    notebook:
        "../notebooks/binary_phase_field.py.ipynb"

rule make_config:
    output:
        "results/{id}/config.json"
    input:
        "config/all_permutations.csv"
    log:
        "logs/make_config_{id}.log"
    run:
        extract_config_by_id(wildcards, output[0], log[0])
