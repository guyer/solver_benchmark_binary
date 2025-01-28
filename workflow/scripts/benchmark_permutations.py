import platform
import numpy as np
import pandas as pd
import logging

def get_list_from_file(listf):
    with open(listf, 'r') as f:
        items = f.read().split()
    return items

# https://stackoverflow.com/a/55849527/2019542
logger = logging.getLogger('benchmark_permutations')
fh = logging.FileHandler(str(snakemake.log[0]))
fh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)

try:
    min_size = snakemake.config["size"].get("min", 10)
    max_size = snakemake.config["size"].get("max", 1000000)
    size_steps = snakemake.config["size"].get("steps", 6)

    benchmarks = []
    for name, benchmark in snakemake.config["benchmarks"].items():
        # calculate dimensions that produce steps in orders of magnitude
        # in number of cells for square 2D grids
        dimension = benchmark.get("dimension", 2)
        sizes = np.logspace(np.log10(min_size) / dimension,
                            np.log10(max_size) / dimension,
                            size_steps,
                            dtype=int)**dimension
        df = pd.DataFrame(data=sizes, columns=["numberOfElements"])
        df["benchmark"] = name
        df["hostname"] = platform.node()
        benchmarks.append(df)

    df = pd.concat(benchmarks)
    df.to_csv(snakemake.output[0], index=False)
except Exception as e:
    logger.error(e, exc_info=True)
    raise e
