from itertools import product
import uuid
import platform

rule aggregate_param_sweeps:
    output:
        "results/permutations.json"
    input:
        expand("results/{solversuite}-permutations.json",
               solversuite=config["suites"])
    log:
        "results/permutations.log"
    run:
        concat_json(input, output[0], log[0])

checkpoint add_param_sweep:
    output:
        "results/{solversuite}-permutations.json"
    input:
        "results/{solversuite}-preconditioners.txt",
        "results/{solversuite}-solvers.txt",
    run:
        solvers = get_checkpoint_list(check=checkpoints.list_solvers,
                                      solversuite=wildcards.solversuite)
        preconditioners = get_checkpoint_list(check=checkpoints.list_preconditioners,
                                              solversuite=wildcards.solversuite)

        df = pd.DataFrame(data=list(product(config["benchmarks"],
                                            solvers,
                                            preconditioners,
                                            SIZES)),
                          columns=["benchmark", "solver", "preconditioner", "size"])

        df["uuid"] = [str(uuid.uuid4()) for item in df.iterrows()]
        df = df.set_index("uuid")

        df["suite"] = wildcards.solversuite
        df["hostname"] = platform.node()
        df["version"] = git_version(path=".")
        df["fipy_version"] = git_version(path=config["fipy_path"])

        df.to_json(output[0])

checkpoint list_solvers:
    output:
        temp("results/{solversuite}-solvers.txt")
    conda:
       "benchmark_{solversuite}"
    log:
        "results/{solversuite}-solvers.log"
    shell:
        "FIPY_SOLVERS={wildcards.solversuite}"
        " python workflow/scripts/solvers.py"
        " > {output} 2> {log}"

checkpoint list_preconditioners:
    output:
        temp("results/{solversuite}-preconditioners.txt")
    conda:
       "benchmark_{solversuite}"
    log:
        "results/{solversuite}-preconditioners.log"
    shell:
        "FIPY_SOLVERS={wildcards.solversuite}"
        " python workflow/scripts/preconditioners.py"
        " > {output} 2> {log}"
