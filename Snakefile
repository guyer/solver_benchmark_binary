from snakemake.utils import Paramspace
import numpy as np
import pandas as pd
from itertools import product

VERSIONS, FIPYVERSIONS = glob_wildcards("results/{version,[A-Za-z0-9]+}/{fipyversion,[^/]+}/")
# calculate dimensions that produce six orders of magnitude in number of cells
SIZES = (10**(np.arange(1, 6.5, 1.)/2)).round().astype(int)**2

# def versions():
#     return shell("git rev-parse --short HEAD")

# VERSIONS += [versions()]
# FIPYVERSIONS += ["3.4.5+7.g3a4f82063"]

rule versions:
    input:
    output:
    shell:
        "git rev-parse --short HEAD"
        #> {output}"

rule all:
    input:
        expand("results/{version}/{fipyversion}/nucleation.py",
               zip, version=VERSIONS, fipyversion=FIPYVERSIONS)

rule current_version:
    input:
    output:
        "current_version.txt"
    conda:
        "fipy_solver_benchmarking_311"
    shell:
        "echo (mkdir -p results/$(git rev-parse --short HEAD)/"
        "$(python -c \"import fipy; print(fipy.__version__)\")) > {output}"
        
rule ipynb2py:
    input:
        "codes/notebooks/{notebook}.ipynb"
    output:
        "{path}/{notebook}.py"
    conda:
        "fipy_solver_benchmarking_311"
    shell:
        "jupyter nbconvert {input} --to python --output-dir={wildcards.path}"

checkpoint list_solvers:
    output:
        "results/{path}/{solversuite}/solvers.txt"
    conda:
       "benchmark_{solversuite}"
    shell:
        "FIPY_SOLVERS={wildcards.solversuite} python codes/scripts/solvers.py > {output}"

checkpoint list_preconditioners:
    output:
        "results/{path}/{solversuite}/preconditioners.txt"
    conda:
       "benchmark_{solversuite}"
    shell:
        "FIPY_SOLVERS={wildcards.solversuite} python codes/scripts/preconditioners.py > {output}"

# checkpoint params:
#     output:
#         "results/{path}/{solversuite}/params.csv"
#     input:
#         "results/{path}/{solversuite}/preconditioners.txt",
#         "results/{path}/{solversuite}/solvers.txt",
#     run:
#         s = checkpoints.solvers
#         solve_file = s.get(path=wildcards.path,
#                            solversuite=wildcards.solversuite).output[0]
#         with open(solve_file, 'r') as f:
#             solvers = f.read().split()
# 
#         p = checkpoints.preconditioners
#         precon_file = p.get(path=wildcards.path,
#                             solversuite=wildcards.solversuite).output[0]
#         with open(p.get(path=wildcards.path,
#                         solversuite=wildcards.solversuite).output[0], 'r') as f:
#             preconditioners = f.read().split()
# 
#         df = pd.DataFrame(data=list(product(solvers, preconditioners, SIZES)),
#                           columns=["solver", "preconditioner", "size"])
# 
#         df.to_csv(output[0], index=False)

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

rule all_solvers:
    output:
        "results/{path}/{solversuite}/all_solvers.csv"
    input:
        "results/{path}/{solversuite}/solvers.txt",
        "results/{path}/{solversuite}/preconditioners.txt",
        get_solvers
    run:
        li = [pd.read_csv(fname, index=False) for fname in input]
        df = pd.concat(li, ignore_index=True)
        df.to_csv(output[0], index=False)

rule all_preconditioners:
    output:
        "results/{path}/{solversuite}/solver~{solver}/all_preconditioners.csv"
    input:
        "results/{path}/{solversuite}/preconditioners.txt",
        get_preconditioners
    run:
        li = [pd.read_csv(fname, index=False) for fname in input]
        df = pd.concat(li, ignore_index=True)
        df.to_csv(output[0], index=False)

rule all_sizes:
    output:
        "results/{path}/all_sizes.csv"
    input:
        get_sizes
    run:
        li = [pd.read_csv(fname, index=False) for fname in input]
        df = pd.concat(li, ignore_index=True)
        df.to_csv(output[0], index=False)

rule extract_times:
    output:
        "results/{path}/{solversuite}/solver~{solver}/preconditioner~{preconditioner}/size~{size}/solver.csv"
    input:
        "results/{path}/{solversuite}/solver~{solver}/preconditioner~{preconditioner}/size~{size}/solver.log"
    shell:
        "touch {output[0]}"
#     notebook:
#         "codes/notebooks/extract.py.ipynb"

rule solve:
    output:
        "results/{version}/{fipyversion}/{script}/{platform}/{solversuite}/"
        "solver~{solver}/preconditioner~{preconditioner}/size~{size}/solver.log"
    input:
        "results/{version}/{fipyversion}/{script}.py"
    conda:
        "benchmark_{solversuite}"
    shell:
        "FIPY_SOLVERS={wildcards.solversuite} python {input[0]}"
        " --solver={wildcards.solver}"
        " --preconditioner={wildcards.preconditioner}"
        " --numberOfElements={wildcards.size}"
        " --output=results/{wildcards.version}/{wildcards.fipyversion}/"
        "{wildcards.script}/{wildcards.platform}/{wildcards.solversuite}/"
        "solver~{wildcards.solver}/preconditioner~{wildcards.preconditioner}/"
        "size~{wildcards.size}"

rule plot:
    output:
        "results/{path}/{solversuite}/all.png"
    input:
        "results/{path}/{solversuite}/params.csv",
        get_params
