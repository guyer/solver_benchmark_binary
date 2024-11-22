checkpoint aggregate_permutations:
    output:
        "config/all_permutations.csv"
    input:
        expand("config/fipy~{rev}/{suite}_permutations.csv",
               rev=config["fipy_revs"],
               suite=config["suites"])
    log:
        "logs/aggregate_permutations.log"
    run:
        concat_csv(input, output[0], log[0])

rule rev_and_suite_permutations:
    output:
        "config/fipy~{rev}/{suite}_permutations.csv"
    input:
        preconditioners=get_preconditioners,
        solvers=get_solvers,
        clone="resources/fipy~{rev}/repo/"
    conda:
        "../envs/snakemake.yml"
    log:
        "logs/fipy~{rev}/{suite}_permutations.log"
    script:
        "../scripts/rev_and_suite_permutations.py"
