import os
import pandas as pd

def concat_csv(input, output, log):
    try:
        li = [pd.read_csv(fname, index_col=False) for fname in input]
        if li:
            df = pd.concat(li, ignore_index=True)
        else:
            df = pd.DataFrame()
        df.to_csv(output, index=False)
    except Exception as e:
        with open(log, 'w') as f:
            f.write(repr(e))
        raise e

def concat_json(input, output, log):
    try:
        li = [pd.read_json(fname) for fname in input]
        if li:
            df = pd.concat(li, ignore_index=True)
        else:
            df = pd.DataFrame()
        df.to_json(output)
    except Exception as e:
        with open(log, 'w') as f:
            f.write(repr(e))
        raise e

def get_conda_environment_from_id(wildcards):
    permutations = get_all_permutations(wildcards)
    rev = permutations.loc[wildcards.id, 'fipy_rev']
    suite = permutations.loc[wildcards.id, 'suite']
    return f"../envs/fipy~{rev}/benchmark_{suite}.yml"

def get_benchmark(wildcards):
    permutations = get_all_permutations(wildcards)
    benchmarks = permutations.loc[wildcards.id, 'benchmark']
    return f"workflow/scripts/{benchmarks}.py"

def get_all_permutation_ids(wildcards):
    df = get_all_permutations(wildcards)

    return df.index

def get_repo(wildcards):
    path = "../resources/fipy~{wildcards.rev}/repo/".format(wildcards=wildcards)
    return os.path.join(workflow.basedir, path)

def get_all_permutations(wildcards):
    path = checkpoints.aggregate_permutations.get().output[0]
    if exists(path):
        df = pd.read_csv(path,
                         index_col="uuid")
    else:
        df = pd.DataFrame(columns=["benchmark",
                                   "solver",
                                   "preconditioner",
                                   "size",
                                   "uuid",
                                   "fipy_rev",
                                   "fipy_version",
                                   "suite",
                                   "hostname"])
        df.set_index("uuid")

    return df

def extract_config_by_id(wildcards, output, log):
    import logging

    # https://stackoverflow.com/a/55849527/2019542
    logger = logging.getLogger('make_config')
    fh = logging.FileHandler(str(log))
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    try:
        permutations = get_all_permutations(wildcards)
        permutations.loc[wildcards.id].to_json(output)
    except Exception as e:
        logger.error(e, exc_info=True)
        raise e

def read_config(path):
    import json

    with open(path, 'r') as f:
        return json.load(f)

def get_preconditioners(wildcards):
    return checkpoints.list_preconditioners.get(rev=wildcards.rev,
                                                suite=wildcards.suite).output[0]

def get_solvers(wildcards):
    return checkpoints.list_solvers.get(rev=wildcards.rev,
                                        suite=wildcards.suite).output[0]
