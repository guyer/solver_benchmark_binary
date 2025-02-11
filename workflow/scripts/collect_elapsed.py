import numpy as np
import pandas as pd
from pathlib import Path

p = pd.read_csv("config/all_permutations.csv")

def get_elapsed(row):
    benchmark = f"benchmarks/fipy~{row['fipy_rev']}/suite~{row['suite']}/benchmark-{row['id']}.tsv"
    benchmark = Path(benchmark)
    if benchmark.exists():
        return pd.read_table(benchmark).iloc[0]["s"]
    else:
        return np.nan

p["elapsed"] = p.apply(get_elapsed, axis=1)

p.to_csv("elapsed_permuations.csv")
