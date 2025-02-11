import numpy as np
import pandas as pd
from pathlib import Path

from plot_permutations import plot_all

p = pd.read_csv("config/all_permutations.csv")

def get_elapsed(row):
    benchmark = f"benchmarks/fipy~{row['fipy_rev']}/suite~{row['suite']}/benchmark-{row['id']}.tsv"
    benchmark = Path(benchmark)
    if benchmark.exists():
        return pd.read_table(benchmark).iloc[0]["s"]
    else:
        return np.nan

# merge elapsed time with permutations
p["elapsed_seconds"] = p.apply(get_elapsed, axis=1)

# make digestible by plot_all()
p["converged"] = True
p = p.rename(columns={"suite": "package.solver", "solver": "solver_class"})

plot_all(p, "benchmark.png")
