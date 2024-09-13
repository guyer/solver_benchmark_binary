def get_solvers(wildcards):
    s = checkpoints.list_solvers
    with open(s.get(path=wildcards.path,
                    solversuite=wildcards.solversuite).output[0], 'r') as f:
        solvers = f.read().split()
    return expand(f"results/{wildcards.path}/suite~{wildcards.solversuite}/solver~{{solvers}}/all_preconditioners.csv",
                  solvers=solvers)

def get_preconditioners(wildcards):
    p = checkpoints.list_preconditioners
    with open(p.get(path=wildcards.path,
                    solversuite=wildcards.solversuite).output[0], 'r') as f:
        preconditioners = f.read().split()
    return expand(f"results/{wildcards.path}/suite~{wildcards.solversuite}/solver~{wildcards.solver}/preconditioner~{{preconditioners}}/all_sizes.csv",
                  preconditioners=preconditioners)

def get_sizes(wildcards):
    return expand(f"results/{wildcards.path}/size~{{sizes}}/solver.csv",
                  sizes=SIZES)

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
        df = pd.concat(li, ignore_index=True)
        df.to_csv(output, index=False)
    except Exception as e:
        with open(log, 'w') as f:
            f.write(e)
        raise

def git_version(path):
    import subprocess
    result = subprocess.run(["git", "-C", path, "describe", "--always", "--dirty"],
                            capture_output=True, text=True)
    return result.stdout.strip()

def tatanka(wildcards):
    self_version = git_version(path=".")
    fipy_version = git_version(path=FIPY_PATH)
    return expand(f"results/benchmark~{{benchmark}}/all_suites.csv",
                  benchmark=BENCHMARKS)
