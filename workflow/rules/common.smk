def get_suites(wildcards):
    suites, = glob_wildcards("results/benchmark~{wildcards.benchmark}/suite~{solversuite,[A-Za-z0-9]+}/all_solvers.csv")
    return expand(f"results/{wildcards.path}/suite~{{solversuite}}/all_solvers.csv",
                  solversuite=suites)

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
    return expand(f"results/{wildcards.path}/size~{{sizes}}/all_hostnames.csv",
                  sizes=SIZES)

def get_hostnames(wildcards):
    hostnames, = glob_wildcards("results/{wildcards.path}/hostname~{hostname,[A-Za-z0-9]+}/all_selfversions.csv")
    return expand(f"results/{wildcards.path}/hostname~{{hostname}}/all_selfversions.csv",
                  hostname=hostnames)

def get_selfversions(wildcards):
    selfversions, = glob_wildcards("results/{wildcards.path}/self~{selfversion,[A-Za-z0-9]+}/all_fipyversions.csv")
    return expand(f"results/{wildcards.path}/self~{{selfversion}}/all_fipyversions.csv",
                  selfversion=selfversions)

def get_fipyversions(wildcards):
    fipyversions, = glob_wildcards("results/{wildcards.path}/fipy~{fipyversion,[A-Za-z0-9]+}/solver.csv")
    return expand(f"results/{wildcards.path}/fipy~{{fipyversion}}/solver.csv",
                  fipyversion=fipyversions)

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
