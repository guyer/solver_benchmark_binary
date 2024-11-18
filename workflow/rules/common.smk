def get_checkpoint_list(check, rev, suite):
    with open(check.get(rev=rev, suite=suite).output[0], 'r') as f:
        items = f.read().split()
    return items

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
        raise

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
        raise

def git_version(path):
    import subprocess
    result = subprocess.run(["git", "-C", path, "describe", "--always", "--dirty"],
                            capture_output=True, text=True)
    return result.stdout.strip()

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

def get_all_permutations(wildcards):
    path = checkpoints.aggregate_param_sweeps2.get().output[0]
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

def read_config(path):
    import json

    with open(path, 'r') as f:
        return json.load(f)
