# def versions():
#     return shell("git rev-parse --short HEAD")

def get_solvers(wildcards):
    s = checkpoints.list_solvers
    with open(s.get(path=wildcards.path,
                    solversuite=wildcards.solversuite).output[0], 'r') as f:
        solvers = f.read().split()
    return expand(f"results/{wildcards.path}/{wildcards.solversuite}/solver~{{solvers}}/all_preconditioners.csv",
                  solvers=solvers)

def get_preconditioners(wildcards):
    p = checkpoints.list_preconditioners
    with open(p.get(path=wildcards.path,
                    solversuite=wildcards.solversuite).output[0], 'r') as f:
        preconditioners = f.read().split()
    return expand(f"results/{wildcards.path}/{wildcards.solversuite}/solver~{wildcards.solver}/preconditioner~{{preconditioners}}/all_sizes.csv",
                  preconditioners=preconditioners)

def get_sizes(wildcards):
    return expand(f"results/{wildcards.path}/size~{{sizes}}/solver.csv",
                  sizes=SIZES)

def get_params(wildcards):
    p = checkpoints.params
    param_file = p.get(path=wildcards.path,
                       solversuite=wildcards.solversuite).output[0]
    paramspace = Paramspace(pd.read_csv(param_file))
    return expand(f"results/{wildcards.path}/{wildcards.solversuite}/{{params}}/solver.log",
                  params=paramspace.instance_patterns)
