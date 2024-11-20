import sys
import pandas as pd

assert len(sys.argv) == 1
input, = sys.argv

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

df2.to_json(sys.stdout, date_format="iso")
