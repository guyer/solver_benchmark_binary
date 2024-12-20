rule list_solvers:
    localrule: True
    output:
        "results/fipy~{rev}/suite~{suite}/solvers.txt"
    input:
        "results/fipy~{rev}/suite~{suite}/environment.yml"
    conda:
        get_conda_environment
    log:
        "logs/fipy~{rev}/suite~{suite}/list_solvers.log"
    script:
        "../scripts/list_solvers.py"

rule list_preconditioners:
    localrule: True
    output:
        "results/fipy~{rev}/suite~{suite}/preconditioners.txt"
    input:
        "results/fipy~{rev}/suite~{suite}/environment.yml"
    conda:
        get_conda_environment
    log:
        "logs/fipy~{rev}/suite~{suite}/list_preconditioners.log"
    script:
        "../scripts/list_preconditioners.py"

checkpoint render_conda_template:
    localrule: True
    output:
        "results/fipy~{rev}/suite~{suite}/environment.yml"
    input:
        template="workflow/envs/benchmark_{suite}.yml",
    log:
        "logs/fipy~{rev}/suite~{suite}/render_conda_template.log"
    template_engine:
        "yte"
