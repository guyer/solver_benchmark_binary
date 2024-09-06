from snakemake.utils import Paramspace
import pandas as pd

VERSIONS, FIPYVERSIONS = glob_wildcards("results/{version,[A-Za-z0-9]+}/{fipyversion,[^/]+}/")
SIZES = [9, 99, 999]

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

checkpoint make_params:
    output:
        "results/{path}/{solversuite}/params.csv"
    conda:
        "fipy_solver_benchmarking_311"
    shell:
        "FIPY_SOLVERS={wildcards.solversuite} python codes/scripts/params.py > {output}"

def get_params(wildcards):
    csv_file = checkpoints.make_params.get(path=wildcards.path,
                                           solversuite=wildcards.solversuite).output[0]
    paramspace = Paramspace(pd.read_csv(csv_file))
    return expand(f"results/{wildcards.path}/{wildcards.solversuite}/{{params}}/solver.log",
                  params=paramspace.instance_patterns)

# checkpoint preconditioners:
#     output:
#         "results/{version}/{fipyversion}/{platform}/{solversuite}/preconditioners.txt"
#     conda:
#         "fipy_solver_benchmarking_311"
# #         "benchmark_{solversuite}"
#     shell:
#         "FIPY_SOLVERS={wildcards.solversuite} python codes/scripts/preconditioners.py > {output}"
# 
# checkpoint solvers:
#     output:
#         "results/{version}/{fipyversion}/{platform}/{solversuite}/solvers.txt"
#     conda:
#         "fipy_solver_benchmarking_311"
#         # "benchmark_{solversuite}"
#     shell:
#         "FIPY_SOLVERS={wildcards.solversuite} python codes/scripts/solvers.py > {output}"

rule solve:
    output:
        "results/{version}/{fipyversion}/{script}/{platform}/{solversuite}/"
        "solver~{solver}/preconditioner~{preconditioner}/solver.log"
    input:
        "results/{version}/{fipyversion}/{script}.py"
    conda:
        "benchmark_{solversuite}"
    shell:
        "FIPY_SOLVERS={wildcards.solversuite} python {input[0]} "
        "--solver={wildcards.solver} "
        "--preconditioner={wildcards.preconditioner}"
#         "--numberOfElements={wildcards.size} "

# def get_logs(wildcards):
#     PATH = f"results/{wildcards.version}/{wildcards.fipyversion}/{wildcards.platform}/{wildcards.solversuite}"
#     with open(PATH + "/solvers.txt", 'r') as f:
#         SOLVERS = f.read().split()
#     with open(PATH + "/preconditioners.txt", 'r') as f:
#         PRECONDITIONERS = f.read().split()
#     return expand(PATH+"{solver}/{preconditioner}/{size}", zip, solver=SOLVERS, preconditioner=PRECONDITIONERS,size=SIZES)
    
rule plot:
    output:
        "results/{path}/{solversuite}/all.png"
    input:
        "results/{path}/{solversuite}/params.csv",
        get_params
