rule extract_times:
    output:
        "results/fipy~{rev}/suite~{suite}/{id}/solver.json"
    input:
        log="results/fipy~{rev}/suite~{suite}/{id}/solver.log",
        config="results/fipy~{rev}/suite~{suite}/{id}/config.json"
    conda:
        "../envs/snakemake.yml"
    log:
        "logs/fipy~{rev}/suite~{suite}/{id}/extract_times.log"
    script:
        "../scripts/extract_times.py"
