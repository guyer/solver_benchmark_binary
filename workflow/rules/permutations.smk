from itertools import product
import uuid
import platform
import numpy as np

# calculate dimensions that produce steps in orders of magnitude
# in number of cells for square 2D grids
SIZES = (10**(np.arange(np.log10(config["size"]["min"]),
                        np.log10(config["size"]["max"])+1,
                        1.)
              /2)).round().astype(int)**2

checkpoint aggregate_param_sweeps2:
    output:
        "config/all_permutations.csv"
    input:
        expand("config/fipy~{rev}/all_permutations.csv",
               rev=config["fipy_revs"])
    log:
        "logs/aggregate_param_sweeps2"
    run:
        concat_csv(input, output[0], log[0])

rule aggregate_param_sweeps:
    output:
        "config/fipy~{rev}/all_permutations.csv"
    input:
        expand("config/fipy~{{rev}}/{suite}_permutations.csv",
               suite=config["suites"])
    log:
        "logs/fipy~{rev}/aggregate_param_sweeps.log"
    run:
        concat_csv(input, output[0], log[0])

checkpoint add_param_sweep:
    output:
        "config/fipy~{rev}/{suite}_permutations.csv"
    input:
        preconditioners="clones/fipy~{rev}/{suite}_preconditioners.txt",
        solvers="clones/fipy~{rev}/{suite}_solvers.txt",
        clone="clones/fipy~{rev}/repo/"
    run:
        solvers = get_checkpoint_list(check=checkpoints.list_solvers,
                                      rev=wildcards.rev,
                                      suite=wildcards.suite)
        preconditioners = \
                  get_checkpoint_list(check=checkpoints.list_preconditioners,
                                      rev=wildcards.rev,
                                      suite=wildcards.suite)

        df = pd.DataFrame(data=list(product(config["benchmarks"],
                                            solvers,
                                            preconditioners,
                                            SIZES)),
                          columns=["benchmark",
                                   "solver",
                                   "preconditioner",
                                   "size"])

        df["uuid"] = [str(uuid.uuid4()) for item in df.iterrows()]
        df = df.set_index("uuid")

        df["fipy_rev"] = wildcards.rev
        df["fipy_version"] = git_version(path=input.clone)
        df["suite"] = wildcards.suite
        df["hostname"] = platform.node()

        df.to_csv(output[0])
