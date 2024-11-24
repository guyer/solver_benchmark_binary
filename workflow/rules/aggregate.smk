rule plot_permutations:
    output:
        total="results/plots/fipy~{rev}/{benchmark}/total.png",
        prepare="results/plots/fipy~{rev}/{benchmark}/prepare.png",
        solve="results/plots/fipy~{rev}/{benchmark}/solve.png",
        plots=report(
            directory("results/plots/fipy~{rev}/{benchmark}"),
            patterns=["{name}.png"],
            category="{rev}",
            subcategory="{benchmark}",
            labels={
                "benchmark": "{benchmark}",
                "name": "{name}"
            }
        )
    input:
        "results/total_times.json"
    log:
        "logs/{benchmark}/plot_permutations_{rev}.log"
    conda:
        "../envs/snakemake.yml"
    script:
        "../scripts/plot_permutations.py"

rule total_times:
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
