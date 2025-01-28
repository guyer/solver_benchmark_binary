rule extract_times:
    output:
        "results/fipy~{rev}/suite~{suite}/{id}/solver.json"
    input:
        log="results/fipy~{rev}/suite~{suite}/{id}/solver.log",
    params:
        config=get_config_by_id,
    conda:
        "../envs/snakemake.yml"
    log:
        "logs/fipy~{rev}/suite~{suite}/{id}/extract_times.log"
    script:
        "../scripts/extract_times.py"
