rule total_times:
    localrule: True
    output:
        "results/total_times.json"
    input:
        all="results/all.json",
    log:
        "logs/total_times.log"
    conda:
        "../envs/snakemake.yml"
    script:
        "../scripts/extract_total_times.py"
