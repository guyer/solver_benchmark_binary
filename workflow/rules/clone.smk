rule list_solvers:
    localrule: True
    output:
        "resources/fipy~{rev}/{suite}_solvers.txt"
    input:
        "workflow/envs/fipy~{rev}/benchmark_{suite}.yml"
    conda:
        get_conda_environment
    log:
        "logs/fipy~{rev}/suite~{suite}/list_solvers.log"
    script:
        "../scripts/list_solvers.py"

rule list_preconditioners:
    localrule: True
    output:
        "resources/fipy~{rev}/{suite}_preconditioners.txt"
    input:
        "workflow/envs/fipy~{rev}/benchmark_{suite}.yml"
    conda:
        get_conda_environment
    log:
        "logs/fipy~{rev}/suite~{suite}/list_preconditioners.log"
    script:
        "../scripts/list_preconditioners.py"

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
