from itertools import product
import uuid
import platform
import numpy as np
import pandas as pd
import logging

def get_checkpoint_list(listf):
    with open(listf, 'r') as f:
        items = f.read().split()
    return items

# https://stackoverflow.com/a/55849527/2019542
logger = logging.getLogger('rev_and_suite_permutations')
fh = logging.FileHandler(str(snakemake.log[0]))
fh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)

try:
    global checkpoints

    solvers = get_checkpoint_list(snakemake.input.solvers)
    preconditioners = get_checkpoint_list(snakemake.input.preconditioners)

    # calculate dimensions that produce steps in orders of magnitude
    # in number of cells for square 2D grids
    SIZES = (10**(np.arange(np.log10(snakemake.config["size"]["min"]),
                            np.log10(snakemake.config["size"]["max"])+1,
                            1.)
                  /2)).round().astype(int)**2

    df = pd.DataFrame(data=list(product(snakemake.config["benchmarks"],
                                        solvers,
                                        preconditioners,
                                        SIZES)),
                      columns=["benchmark",
                               "solver",
                               "preconditioner",
                               "size"])

    df["uuid"] = [str(uuid.uuid4()) for item in df.iterrows()]
    df = df.set_index("uuid")

    df["fipy_rev"] = snakemake.wildcards.rev
    df["suite"] = snakemake.wildcards.suite
    df["hostname"] = platform.node()

    df.to_csv(snakemake.output[0])
except Exception as e:
    logger.error(e, exc_info=True)
    raise e
