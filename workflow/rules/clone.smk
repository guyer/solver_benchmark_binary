checkpoint list_solvers:
    localrule: True
    output:
        "resources/fipy~{rev}/{suite}_solvers.txt"
    input:
        "workflow/envs/fipy~{rev}/benchmark_{suite}.yml"
    conda:
        "../envs/fipy~{rev}/benchmark_{suite}.yml"
    log:
        "logs/fipy~{rev}/suite~{suite}/list_solvers.log"
    script:
        "../scripts/list_solvers.py"

checkpoint list_preconditioners:
    localrule: True
    output:
        "resources/fipy~{rev}/{suite}_preconditioners.txt"
    input:
        "workflow/envs/fipy~{rev}/benchmark_{suite}.yml"
    conda:
        "../envs/fipy~{rev}/benchmark_{suite}.yml"
    log:
        "logs/fipy~{rev}/suite~{suite}/list_preconditioners.log"
    script:
        "../scripts/list_preconditioners.py"

rule make_conda_env:
    localrule: True
    output:
        "workflow/envs/fipy~{rev}/benchmark_{suite}.yml",
    input:
        template="workflow/envs/benchmark_{suite}.yml",
    log:
        "logs/fipy~{rev}/make_conda_env_{suite}.log"
    template_engine:
        "jinja2"
