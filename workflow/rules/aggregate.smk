import pandas as pd

rule plot_permutations:
    output:
        total="results/plots/fipy~{rev}/total.png",
        prepare="results/plots/fipy~{rev}/prepare.png",
        solve="results/plots/fipy~{rev}/solve.png",
        plots=report(
            directory("results/plots/fipy~{rev}"),
            patterns=["{name}.png"],
            category="{rev}",
            labels={"name": "{name}"}
        )
    input:
        "results/total_times.json"
    log:
        "logs/plot_permutations_{rev}.log"
    conda:
        "../envs/snakemake.yml"
    script:
        "../scripts/plot_permutations.py"

rule total_times:
    output:
        "results/total_times.json"
    input:
        all="results/all.json",
        script="workflow/scripts/extract_total_times.py"
    log:
        "logs/total_times.log"
    conda:
        "../envs/snakemake.yml"
    shell:
        r"""
        python {input.script} {input.all} \
            > {output:q} 2> {log:q}
        """

