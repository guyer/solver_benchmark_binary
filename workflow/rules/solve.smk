rule solve:
    output:
        "results/{id}/solver.log"
    params:
        config=lambda w, input: read_config(input.config),
    conda:
        get_conda_environment_from_id
    log:
        notebook="logs/notebooks/benchmark_{id}.ipynb"
    notebook:
        "../notebooks/binary_phase_field.py.ipynb"
