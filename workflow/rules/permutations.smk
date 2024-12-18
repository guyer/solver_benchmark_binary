rule aggregate_permutations:
    localrule: True
    output:
        "results/all_permutations.csv"
    input:
        expand("results/fipy~{rev}/suite~{suite}/permutations.csv",
               rev=config["fipy_revs"],
               suite=config["suites"])
    log:
        "logs/aggregate_permutations.log"
    run:
        concat_csv(input, output[0], log[0])

checkpoint rev_and_suite_permutations:
    output:
        "results/fipy~{rev}/suite~{suite}/permutations.csv"
    input:
        preconditioners="results/fipy~{rev}/suite~{suite}/preconditioners.txt",
        solvers="results/fipy~{rev}/suite~{suite}/solvers.txt"
    conda:
        "../envs/snakemake.yml"
    log:
        "logs/fipy~{rev}/suite~{suite}/rev_and_suite_permutations.log"
    script:
        "../scripts/rev_and_suite_permutations.py"
