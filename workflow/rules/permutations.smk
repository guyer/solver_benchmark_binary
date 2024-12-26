rule all_permutations:
    output:
        "config/all_permutations.csv"
    input:
        benchmarks="config/benchmark_permutations.csv",
        fipy_revs="config/fipy_permutations.csv"
    run:
        benchmarks = pd.read_csv(input["benchmarks"], index_col=False)
        fipy_revs = pd.read_csv(input["fipy_revs"], index_col=False)
        df = benchmarks.join(fipy_revs, how="cross")
        df.to_csv(output[0], index=False)

rule benchmark_permutations:
    output:
        "config/benchmark_permutations.csv"
    conda:
        "../envs/snakemake.yml"
    log:
        "logs/benchmark_permutations.log"
    script:
        "../scripts/benchmark_permutations.py"

rule fipy_permutations:
    output:
        "config/fipy_permutations.csv"
    input:
        expand("config/fipy~{rev}/suite_permutations.csv",
               rev=config["fipy_revs"])
    log:
        "logs/all_permutations.log"
    run:
        concat_csv(input, output[0], log[0])

rule suite_permutations:
    output:
        "config/fipy~{rev}/suite_permutations.csv"
    input:
        expand("config/fipy~{{rev}}/suite~{suite}/permutations.csv",
               suite=config["suites"])
    log:
        "logs/fipy~{rev}/suite_permutations.csv"
    run:
        concat_csv(input, output[0], log[0])

rule solver_preconditioner_permutations:
    output:
        "config/fipy~{rev}/suite~{suite}/permutations.csv"
    input:
        solvers="config/fipy~{rev}/suite~{suite}/solvers.json",
        preconditioners="config/fipy~{rev}/suite~{suite}/preconditioners.json"
    log:
        "logs/fipy~{rev}/suite~{suite}/solver_preconditioner_permutations.log"
    run:
        solvers = pd.read_json(input["solvers"])
        preconditioners = pd.read_json(input["preconditioners"])

        permutations = solvers.join(preconditioners, how="cross")
        # it makes no sense to precondition LU
        permutations = permutations.query("(solver != 'LinearLUSolver')"
                                          "| (preconditioner == 'none')")

        permutations["suite"] = wildcards["suite"]
        permutations["fipy_rev"] = wildcards["rev"]
        permutations.to_csv(output[0], index=False)

rule get_solvers_preconditioners:
    output:
        solvers="config/fipy~{rev}/suite~{suite}/solvers.json",
        preconditioners="config/fipy~{rev}/suite~{suite}/preconditioners.json"
    input:
        "results/fipy~{rev}/suite~{suite}/environment.yml"
    conda:
        get_conda_environment
    log:
        "logs/fipy~{rev}/suite~{suite}/get_solvers_preconditioners.log"
    script:
        "../scripts/get_solvers_preconditioners.py"
