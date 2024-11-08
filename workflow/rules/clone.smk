from textwrap import dedent

checkpoint list_solvers:
    output:
        "clones/fipy~{rev}/{suite}_solvers.txt"
    input:
        conda="workflow/envs/fipy~{rev}/benchmark_{suite}.yml",
        script="workflow/scripts/solvers.py"
    conda:
       "envs/fipy~{rev}/benchmark_{suite}.yml"
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
        "clones/fipy~{rev}/{suite}_preconditioners.txt"
    input:
        conda="workflow/envs/fipy~{rev}/benchmark_{suite}.yml",
        script="workflow/scripts/preconditioners.py"
    conda:
       "envs/fipy~{rev}/benchmark_{suite}.yml"
    log:
        "logs/fipy~{rev}/suite~{suite}/list_preconditioners.log"
    shell:
        r"""
        FIPY_SOLVERS={wildcards.suite} \
            python {input.script} \
            > {output:q} 2> {log:q}
        """

rule miney:
    output:
        repo=directory("clones/fipy~{rev}")
    input:
        env=expand("workflow/envs/fipy~{{rev}}/benchmark_{suite}.yml",
                   suite=config["suites"]),
        post=expand("workflow/envs/fipy~{{rev}}/benchmark_{suite}.post-deploy.yml",
                    suite=config["suites"])
    params:
        fipy_repo=lambda _: config["fipy_url"]
    shell:
        """
        git clone --filter=blob:none {params.fipy_repo:q} {output[repo]:q}
        pushd clones/fipy_{wildcards.rev:q}
        git checkout {wildcards.rev:q}
        popd
        """

rule meeny:
    output:
        env="workflow/envs/fipy~{rev}/benchmark_{suite}.yml",
        post="workflow/envs/fipy~{rev}/benchmark_{suite}.post-deploy.yml"
    input:
        env="workflow/envs/benchmark_{suite}.yml",
    shell:
        dedent("""
        cp {input.env:q} {output.env:q}
        
        cat <<EOF >> {output.post:q}
        #!/usr/bin/env bash
        
        pip install --editable "clones/fipy_{wildcards.rev}"
        EOF        
        """)
        