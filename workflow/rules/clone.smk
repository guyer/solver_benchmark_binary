checkpoint render_conda_template:
    localrule: True
    output:
        "workflow/envs/fipy_benchmark_petsc.yml",
    input:
        template="workflow/envs/benchmark_petsc.yml",
    log:
        "logs/render_conda_template_petsc.log"
    template_engine:
        "yte"
