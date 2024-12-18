checkpoint aggregate_permutations:
    output:
        "config/all_permutations.csv"
    input:
        "workflow/envs/fipy_benchmark_petsc.yml"
    conda:
        "../envs/snakemake.yml"
    log:
        "logs/fipy_petsc_permutations.log"
    script:
        "../scripts/rev_and_suite_permutations.py"
