checkpoint render_conda_template:
    localrule: True
    output:
        "workflow/envs/fipy~{rev}/benchmark_{suite}.yml",
    input:
        template="workflow/envs/benchmark_{suite}.yml",
    log:
        "logs/fipy~{rev}/render_conda_template_{suite}.log"
    template_engine:
        "yte"
