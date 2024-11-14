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
        df = df.query(f"fipy_version == '{wildcards.rev}'")
        plot_all(df, output.total, ymin=1e-2, ymax=1e2)
        plot_all(df, output.prepare, data_set="prepare_seconds", ylabel="prepare time", ymin=1e-2, ymax=1e2)
        plot_all(df, output.solve, data_set="solve_seconds", ylabel="solve time", ymin=1e-4, ymax=1e2)

def extract_total_times(input, output):
    df = pd.read_json(input, convert_dates=["time_stamp"])
    df["solve_time"] = pd.to_timedelta(df["solve_time"])
    df["prepare_time"] = pd.to_timedelta(df["prepare_time"])

    gb = df.groupby("simulation_id")
    prepare_time = gb["prepare_time"].sum()
    solve_time = gb["solve_time"].sum()

    df2 = df[df["state"].isin(["START", "END"])].copy()
    df2["time_delta"] = df2["time_stamp"].diff()
    df2 = df2[df2["state"] == "END"].set_index("simulation_id")
    df2["prepare_time"] = prepare_time
    df2["solve_time"] = solve_time

    second = pd.Timedelta("1 s")
    df2["elapsed_seconds"] = df2["time_delta"] / second
    df2["prepare_seconds"] = df2["prepare_time"] / second
    df2["solve_seconds"] = df2["solve_time"] / second

    df2.sort_values("numberOfElements").to_json(output, date_format="iso")

    df2.to_json(output, date_format="iso")

rule total_times:
    output:
        "results/total_times.json"
    input:
        "results/all.json"
    run:
        extract_total_times(input[0], output[0])

