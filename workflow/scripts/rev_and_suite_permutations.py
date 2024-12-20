from itertools import product
import platform
import numpy as np
import pandas as pd
import logging

def get_list_from_file(listf):
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
    solvers = get_list_from_file(snakemake.input.solvers)
    preconditioners = get_list_from_file(snakemake.input.preconditioners)

    # calculate dimensions that produce steps in orders of magnitude
    # in number of cells for square 2D grids
    min_size = snakemake.config["size"].get("min", 10)
    max_size = snakemake.config["size"].get("max", 1000000)
    size_steps = snakemake.config["size"].get("steps", 6)
    dimension = 2
    sizes = np.logspace(np.log10(min_size) / dimension,
                        np.log10(max_size) / dimension,
                        size_steps,
                        dtype=int)**dimension

    df = pd.DataFrame(data=list(product(snakemake.config["benchmarks"],
                                        solvers,
                                        preconditioners,
                                        sizes)),
                      columns=["benchmark",
                               "solver",
                               "preconditioner",
                               "size"])

    df["fipy_rev"] = snakemake.wildcards.rev
    df["suite"] = snakemake.wildcards.suite
    df["hostname"] = platform.node()

    df.to_csv(snakemake.output[0], index_label="index")
except Exception as e:
    logger.error(e, exc_info=True)
    raise e
