from itertools import product
import uuid
import platform
import numpy as np
import logging

# calculate dimensions that produce steps in orders of magnitude
# in number of cells for square 2D grids
SIZES = (10**(np.arange(np.log10(config["size"]["min"]),
                        np.log10(config["size"]["max"])+1,
                        1.)
              /2)).round().astype(int)**2

checkpoint aggregate_permutations:
    output:
        "config/all_permutations.csv"
    input:
        expand("config/fipy~{rev}/{suite}_permutations.csv",
               rev=config["fipy_revs"],
               suite=config["suites"])
    log:
        "logs/aggregate_permutations.log"
    run:
        concat_csv(input, output[0], log[0])

rule rev_and_suite_permutations:
    output:
        "config/fipy~{rev}/{suite}_permutations.csv"
    input:
        preconditioners="resources/fipy~{rev}/{suite}_preconditioners.txt",
        solvers="resources/fipy~{rev}/{suite}_solvers.txt",
        clone="resources/fipy~{rev}/repo/"
    log:
        "logs/fipy~{rev}/{suite}_permutations.log"
    run:
        # https://stackoverflow.com/a/55849527/2019542
        logger = logging.getLogger('rev_and_suite_permutations')
        fh = logging.FileHandler(str(log))
        fh.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        logger.addHandler(fh)

        try:
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
        except Exception as e: 
            logger.error(e, exc_info=True)
