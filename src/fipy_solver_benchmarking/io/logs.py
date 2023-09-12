import json
import numpy as np
import pandas as pd
import uuid

def read_events(fname):
    success_statuses = ['KSP_CONVERGED_RTOL', 'SCIPY_SUCCESS', 'AZ_normal', 
                        'KSP_CONVERGED_ITS', 'KSP_CONVERGED_RTOL', 
                        'Pysparse_CONVERGED_RTOL', 'AMGX_SOLVE_SUCCESS']
    events = []
    state = {}
    versions = {}

    with open(fname, 'r') as f:
        begin_time = np.nan
        solve_time = np.nan
        simulation_id = np.nan
        for line in f:
            entries = line.split("|")
            # Python version uses '|', too
            entries = entries[:4] + ["|".join(entries[4:])]

            (time_stamp,
             level,
             logger,
             function,
             msg) = [s.strip() for s in entries]

            if (logger, function) == ("fipy.solvers", "<module>"):
                pass
            elif (level, logger, function) == ("DEBUG", "fipy", "<module>"):
                versions = json.loads(msg)
                # # fix name collision
                # versions["package"]["solver_suite"] = versions["package"]["solver"]
                # del versions["package"]["solver"]
            elif (level, logger, function) == ("INFO", "fipy", "<module>"):
                # obsolete log format
                versions = {"package": json.loads(msg)}
            elif (level, function) == ("DEBUG", "<module>"):
                if msg.startswith("result stored in"):
                    continue
                state = json.loads(msg)
                state["logfile"] = fname
                state["time_stamp"] = time_stamp
                state["solver_class"] = logger.split('.')[-1]
                state["solve_time"] = str(solve_time)
                if state["state"] == "START":
                    # create unique id for each simulation run
                    simulation_id = uuid.uuid4()
                    converged = False
                state["simulation_id"] = simulation_id
                state["converged"] = converged
                state.update(versions)
                events.append(state.copy())
            elif (level, function) == ("DEBUG", "_setConvergence"):
                state["state"] = "SWEEP"
                event = json.loads(msg)
                event.update(state.copy())
                event["time_stamp"] = time_stamp
                event["solver_class"] = logger.split('.')[-1]
                event["solve_time"] = str(solve_time)
                converged = event["status_name"] in success_statuses
                events.append(event)
                solve_time = np.nan
            elif (level, function) == ("DEBUG", "_solve_"):
                if msg == "BEGIN solve":
                    begin_time = pd.to_datetime(time_stamp)
                    solve_time = np.nan
                elif msg == "END solve":
                    solve_time = pd.to_datetime(time_stamp) - begin_time
                    begin_time = np.nan

    return events

def events2df(events):
    df = pd.json_normalize(events)
    df["time_stamp"] = pd.to_datetime(df["time_stamp"])
    df["solve_time"] = pd.to_timedelta(df["solve_time"])    
    df.loc[df["preconditioner"].isna()
           | (df["preconditioner"] == "NoneType"), "preconditioner"] = "unpreconditioned"

    return df

def extract_total_times(df):
    solve_time = df.groupby("simulation_id")["solve_time"].sum()
    
    df2 = df[df["state"].isin(["START", "END"])].copy()
    df2["time_delta"] = df2["time_stamp"].diff()
    df2 = df2[df2["state"] == "END"].set_index("simulation_id")
    df2["solve_time"] = solve_time
    
    df2["elapsed_seconds"] = df2["time_delta"] / pd.Timedelta("00:00:01")
    df2["solve_seconds"] = df2["solve_time"] / pd.Timedelta("00:00:01")

    return df2.sort_values("numberOfElements")

def extract_sweep_times(df):
    df2 = df.copy()
    df2["time_delta"] = df["time_stamp"].diff()
    df2["elapsed_seconds"] = df2["time_delta"] / pd.Timedelta("00:00:01")

    return df2[~df2["state"].isin(["START", "END"])]