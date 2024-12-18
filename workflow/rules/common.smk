import os
import pandas as pd

def get_conda_environment_from_id(wildcards):
    return "../envs/fipy_benchmark_petsc.yml"

def get_all_permutation_ids(wildcards):
    df = get_all_permutations(wildcards)

    return df.index

def get_all_permutations(wildcards):
    path = checkpoints.aggregate_permutations.get().output[0]
    return pd.read_csv(path,
                       index_col="uuid")
