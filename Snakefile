from snakemake.utils import Paramspace
import pandas as pd

VERSIONS, FIPYVERSIONS = glob_wildcards("results/{version,[A-Za-z0-9]+}/{fipyversion,[^/]+}/")
SIZES = [9, 99, 999]

def versions():
    return shell("git rev-parse --short HEAD")

# VERSIONS += [versions()]
# FIPYVERSIONS += ["3.4.5+7.g3a4f82063"]

# rule versions:
#     input:
#     output: 
#     shell:
#         "git rev-parse --short HEAD > {output}"

rule all:
    input: expand("results/{version}/{fipyversion}/nucleation.py", zip, version=VERSIONS, fipyversion=FIPYVERSIONS)

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
        "results/{version}/{fipyversion}/{platform}/{script}/{solversuite}/params.csv"
    conda:
        "fipy_solver_benchmarking_311"
    shell:
        "FIPY_SOLVERS={wildcards.solversuite} python codes/scripts/params.py > {output}"

def get_params(wildcards):
    paramspace = Paramspace(pd.read_csv(checkpoints.make_params.get(version=wildcards.version,
                                                                    fipyversion=wildcards.fipyversion,
                                                                    script=wildcards.script,
                                                                    platform=wildcards.platform,
                                                                    solversuite=wildcards.solversuite).output[0]))
    return expand(f"results/{wildcards.version}/{wildcards.fipyversion}/{wildcards.script}/{wildcards.platform}/{wildcards.solversuite}/{{params}}/solver.log",
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
    shell:
        "FIPY_SOLVERS={wildcards.solversuite} python {input[0]} "
        "--solver={wildcards.solver} "
        "--preconditioner={wildcards.preconditioner}"
#         "--numberOfElements={wildcards.size} "

# rule plotty:
#     output:
#         "results/{version}/{fipyversion}/{platform}/{solversuite}/all.png"
#     input:
#         "results/{version}/{fipyversion}/{platform}/{solversuite}/"
#         "{solver}/{preconditioner}/{size}/solver.log",
#         solvers="results/{version}/{fipyversion}/{platform}/{solversuite}/solvers.txt"

# def get_logs(wildcards):
#     PATH = f"results/{wildcards.version}/{wildcards.fipyversion}/{wildcards.platform}/{wildcards.solversuite}"
#     with open(PATH + "/solvers.txt", 'r') as f:
#         SOLVERS = f.read().split()
#     with open(PATH + "/preconditioners.txt", 'r') as f:
#         PRECONDITIONERS = f.read().split()
#     return expand(PATH+"{solver}/{preconditioner}/{size}", zip, solver=SOLVERS, preconditioner=PRECONDITIONERS,size=SIZES)
    
rule plot:
    output:
        "results/{version}/{fipyversion}/{script}/{platform}/{solversuite}/all.png"
    input:
        "results/{version}/{fipyversion}/{script}/{platform}/{solversuite}/params.csv",
        get_params
