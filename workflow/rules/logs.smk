rule extract_times:
    output:
        "results/{id}/solver.json"
    input:
        log="results/{id}/solver.log"
    conda:
        "../envs/snakemake.yml"
    log:
        "results/{id}/extract_times.log"
    script:
        "../scripts/extract_times.py"
