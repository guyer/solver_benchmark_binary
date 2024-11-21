from pathlib import Path
from textwrap import dedent

checkpoint list_solvers:
    output:
        "resources/fipy~{rev}/{suite}_solvers.txt"
    input:
        conda="workflow/envs/fipy~{rev}/benchmark_{suite}.yml",
        script="workflow/scripts/solvers.py",
    conda:
        "../envs/fipy~{rev}/benchmark_{suite}.yml"
    log:
        "logs/fipy~{rev}/suite~{suite}/list_solvers.log"
    shell:
        r"""
        FIPY_SOLVERS={wildcards.suite} \
            python {input.script} \
            > {output:q} 2> {log:q}
        """

checkpoint list_preconditioners:
    output:
        "resources/fipy~{rev}/{suite}_preconditioners.txt"
    input:
        conda="workflow/envs/fipy~{rev}/benchmark_{suite}.yml",
        script="workflow/scripts/preconditioners.py",
    conda:
        "../envs/fipy~{rev}/benchmark_{suite}.yml"
    log:
        "logs/fipy~{rev}/suite~{suite}/list_preconditioners.log"
    shell:
        r"""
        FIPY_SOLVERS={wildcards.suite} \
            python {input.script} \
            > {output:q} 2> {log:q}
        """

rule make_conda_env:
    output:
        env="workflow/envs/fipy~{rev}/benchmark_{suite}.yml",
        post="workflow/envs/fipy~{rev}/benchmark_{suite}.post-deploy.sh"
    input:
        env="workflow/envs/benchmark_{suite}.yml",
        repo="resources/fipy~{rev}/repo/"
    conda:
        "../envs/snakemake.yml"
    log:
        "logs/fipy~{rev}/make_conda_env_{suite}.log"
    shell:
        dedent("""
        exec 2> "{log:q}"  # send all stderr from this script to the log file

        cp {input.env:q} {output.env:q}
        
        cat <<EOF >> {output.post:q}
        #!/usr/bin/env bash
        
        pip install --editable "resources/fipy~{wildcards.rev}/repo"
        EOF        
        """)

rule clone_repo:
    output:
        repo=directory("resources/fipy~{rev}/repo/"),
    params:
        fipy_repo=lambda _: config["fipy_url"]
    conda:
        "../envs/snakemake.yml"
    log:
        "logs/fipy~{rev}/clone_repo.log"
    shell:
        """
        exec 2> "{log:q}"  # send all stderr from this script to the log file

        git clone --filter=blob:none {params.fipy_repo:q} {output[repo]:q}
        pushd {output[repo]:q}
        git checkout {wildcards.rev:q}
        popd
        """
