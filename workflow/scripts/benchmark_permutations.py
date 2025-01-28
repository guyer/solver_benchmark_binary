import platform
import numpy as np
import pandas as pd
import logging

def get_list_from_file(listf):
    with open(listf, 'r') as f:
        items = f.read().split()
    return items

# https://stackoverflow.com/a/55849527
# https://stackoverflow.com/a/65732832
logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
)

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
    logging.error(e, exc_info=True)
    raise e
