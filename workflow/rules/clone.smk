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
        env="workflow/envs/fipy~{rev}/benchmark_{suite}.yml",
        post="workflow/envs/fipy~{rev}/benchmark_{suite}.post-deploy.sh"
    input:
        env="workflow/envs/benchmark_{suite}.yml",
        repo=workflow.source_path("resources/fipy~{rev}/repo/")
    conda:
        "../envs/snakemake.yml"
    log:
        "logs/fipy~{rev}/make_conda_env_{suite}.log"
    script:
        "../scripts/make_conda_env.sh"

rule clone_repo:
    localrule: True
    output:
        repo=directory("resources/fipy~{rev}/repo/"),
    params:
        fipy_repo=lambda _: config["fipy_url"]
    conda:
        "../envs/snakemake.yml"
    log:
        "logs/fipy~{rev}/clone_repo.log"
    script:
        "../scripts/clone_repo.sh"
