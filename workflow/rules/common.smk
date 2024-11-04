def get_checkpoint_list(check, solversuite):
    with open(check.get(solversuite=solversuite).output[0], 'r') as f:
        items = f.read().split()
    return items

def get_params(wildcards):
    p = checkpoints.params
    param_file = p.get(path=wildcards.path,
                       solversuite=wildcards.solversuite).output[0]
    paramspace = Paramspace(pd.read_csv(param_file))
    return expand(f"results//{{params}}/solver.log",
                  params=paramspace.instance_patterns)

def concat_csv(input, output, log):
    try:
        li = [pd.read_csv(fname, index=False) for fname in input]
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

def get_conda_environment(wildcards):
    return f"benchmark_{permutations.loc[wildcards.id, 'suite']}"

def get_benchmark(wildcards):
    return f"workflow/scripts/{permutations.loc[wildcards.id, 'benchmark']}.py"

def get_all_plots(wildcards):
    if exists(checkpoints.total_times.get().output[0]):
        df = pd.read_json(checkpoints.total_times.get().output[0])

        gb = df.groupby(by=["fipy_version"])

        plots = expand("results/plots/fipy~{key}/{plot}.png",
                       key=list(gb.groups.keys()),
                       plot=["total", "prepare", "solve"])
    else:
        plots = []

    return plots

def read_config(path):
    import json

    with open(path, 'r') as f:
        return json.load(f)
