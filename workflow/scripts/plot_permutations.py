import pandas as pd
import logging

# https://stackoverflow.com/a/55849527/2019542
logger = logging.getLogger('plot_permutations')
fh = logging.FileHandler(str(snakemake.log[0]))
fh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)

try:
    df = pd.read_json(snakemake.input[0])
    df = df.query(f"fipy_rev == '{snakemake.wildcards.rev}'")
    plot_all(df, snakemake.output.total, ymin=1e-4, ymax=1e2)
    plot_all(df, snakemake.output.prepare, data_set="prepare_seconds", ylabel="prepare time", ymin=1e-4, ymax=1e2)
    plot_all(df, snakemake.output.solve, data_set="solve_seconds", ylabel="solve time", ymin=1e-4, ymax=1e2)
except Exception as e:
    logger.error(e, exc_info=True)
    pass
