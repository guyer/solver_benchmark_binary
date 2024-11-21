from pathlib import Path
from textwrap import dedent

checkpoint list_solvers:
    output:
        "resources/fipy~{rev}/{suite}_solvers.txt"
    input:
#         conda="../envs/fipy~{rev}/benchmark_{suite}.yml",
        conda=(Path(workflow.current_basedir)
               / "../envs/fipy~{rev}/benchmark_{suite}.yml").as_posix(),
        script="workflow/scripts/solvers.py",
#         clone="resources/fipy~{rev}/repo/"
    conda:
        get_conda_environment_from_rev_and_suite
#         "../envs/fipy~{rev}/benchmark_{suite}.ymlA"
#         (Path(workflow.current_basedir)
#              / "../envs/fipy~{rev}/benchmark_{suite}.yml").as_posix()
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
#         conda="../envs/fipy~{rev}/benchmark_{suite}.yml",
        conda=(Path(workflow.current_basedir)
               / "../envs/fipy~{rev}/benchmark_{suite}.yml").as_posix(),
        script="workflow/scripts/preconditioners.py",
#         clone="resources/fipy~{rev}/repo/"
    conda:
        get_conda_environment_from_rev_and_suite
#         "../envs/fipy~{rev}/benchmark_{suite}.ymlB"
    log:
        "logs/fipy~{rev}/suite~{suite}/list_preconditioners.log"
    shell:
        r"""
        FIPY_SOLVERS={wildcards.suite} \
            python {input.script} \
            > {output:q} 2> {log:q}
        """



checkpoint make_conda_env:
    output:
#         env="../envs/fipy~{rev}/benchmark_{suite}.yml",
#         post="../envs/fipy~{rev}/benchmark_{suite}.post-deploy.sh"
        env=(Path(workflow.current_basedir)
             / "../envs/fipy~{rev}/benchmark_{suite}.yml").as_posix(),
        post=(Path(workflow.current_basedir)
              / "../envs/fipy~{rev}"
              / "benchmark_{suite}.post-deploy.sh").as_posix()
    input:
        env="workflow/envs/benchmark_{suite}.yml",
        repo="resources/fipy~{rev}/repo/"
    log:
        "logs/fipy~{rev}/make_conda_env_{suite}.log"
    shell:
        dedent("""
        pwd > {log}

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
    shell:
        """
        git clone --filter=blob:none {params.fipy_repo:q} {output[repo]:q}
        pushd {output[repo]:q}
        git checkout {wildcards.rev:q}
        popd
        """
