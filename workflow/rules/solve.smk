rule solve:
    output:
        "results/{id}/solver.log"
    input:
        "config/all_permutations.csv"
    conda:
        get_conda_environment_from_id
    log:
        notebook="logs/notebooks/benchmark_{id}.ipynb"
    notebook:
        "../notebooks/binary_phase_field.py.ipynb"
