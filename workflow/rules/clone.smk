from textwrap import dedent

checkpoint list_solvers:
    output:
        "clones/fipy~{rev}/{suite}_solvers.txt"
    input:
        conda="workflow/rules/../envs/fipy~{rev}/benchmark_{suite}.yml",
        script="workflow/scripts/solvers.py",
        clone="clones/fipy~{rev}/repo/"
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
        "clones/fipy~{rev}/{suite}_preconditioners.txt"
    input:
        conda="workflow/rules/../envs/fipy~{rev}/benchmark_{suite}.yml",
        script="workflow/scripts/preconditioners.py",
        clone="clones/fipy~{rev}/repo/"
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

rule miney:
    output:
        repo=directory("clones/fipy~{rev}/repo/"),
    input:
        env=expand("workflow/rules/../envs/fipy~{{rev}}/benchmark_{suite}.yml",
                   suite=config["suites"]),
        post=expand("workflow/rules/../envs/fipy~{{rev}}/benchmark_{suite}.post-deploy.sh",
                    suite=config["suites"])
    params:
        fipy_repo=lambda _: config["fipy_url"]
    shell:
        """
        git clone --filter=blob:none {params.fipy_repo:q} {output[repo]:q}
        pushd {output[repo]:q}
        git checkout {wildcards.rev:q}
        popd
        """

rule meeny:
    output:
        env="workflow/rules/../envs/fipy~{rev}/benchmark_{suite}.yml",
        post="workflow/rules/../envs/fipy~{rev}/benchmark_{suite}.post-deploy.sh"
    input:
        env="workflow/envs/benchmark_{suite}.yml"
    log:
        "logs/fipy~{rev}/meeny_{suite}.log"
    shell:
        dedent("""
        cp {input.env:q} {output.env:q}
        
        cat <<EOF >> {output.post:q}
        #!/usr/bin/env bash
        
        pip install --editable "clones/fipy~{wildcards.rev}/repo"
        EOF        
        """)
