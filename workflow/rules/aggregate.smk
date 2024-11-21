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
    run:
        df = pd.read_json(input[0])
        df = df.query(f"fipy_rev == '{wildcards.rev}'")
        plot_all(df, output.total, ymin=1e-2, ymax=1e2)
        plot_all(df, output.prepare, data_set="prepare_seconds", ylabel="prepare time", ymin=1e-2, ymax=1e2)
        plot_all(df, output.solve, data_set="solve_seconds", ylabel="solve time", ymin=1e-4, ymax=1e2)

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

