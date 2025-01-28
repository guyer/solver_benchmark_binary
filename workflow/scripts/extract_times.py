import json
import numpy as np
import pandas as pd
import uuid
import logging

def read_events(fname):
    success_statuses = ['KSP_CONVERGED_RTOL', 'SCIPY_SUCCESS', 'AZ_normal',
                        'KSP_CONVERGED_ITS', 'KSP_CONVERGED_RTOL',
                        'Pysparse_CONVERGED_RTOL', 'AMGX_SOLVE_SUCCESS']
    events = []
    state = {}
    versions = {}

    with open(fname, 'r') as f:
        begin_solve_time = np.nan
        begin_prepare_time = np.nan
        solve_time = np.nan
        prepare_time = np.nan
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

            # parse time format from logger %(asctime)
            # https://docs.python.org/3/library/logging.html#logging.Formatter.formatTime
            time_stamp = pd.to_datetime(time_stamp,
                                        format="%Y-%m-%d %H:%M:%S,%f")

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
                    # Create unique id for each simulation run.
                    # There could be more than one simulation
                    # in a given log file.
                    simulation_id = str(uuid.uuid4())
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
                event["prepare_time"] = str(prepare_time)
                converged = event["status_name"] in success_statuses
                events.append(event)
                solve_time = np.nan
            elif (level, function) == ("DEBUG", "log"):
                if logger.endswith("Convergence") or logger.endswith("Divergence"):
                    state["state"] = "SWEEP"
                    event = json.loads(msg)
                    event.update(state.copy())
                    event["time_stamp"] = time_stamp
                    event["solver_class"] = event["solver"].split("(")[0]
                    event["solve_time"] = str(solve_time)
                    event["prepare_time"] = str(prepare_time)
                    converged = event["status_name"] in success_statuses
                    events.append(event)
                    solve_time = np.nan
            elif (level, function) == ("DEBUG", "_solve_"):
                if msg == "BEGIN solve":
                    begin_solve_time = time_stamp
                    solve_time = np.nan
                elif msg.startswith("END solve - "):
                    solve_time = msg.split('-')[-1]
                    solve_time = pd.to_timedelta(solve_time)
                    begin_solve_time = np.nan
                elif msg == "END solve":
                    solve_time = time_stamp - begin_solve_time
                    begin_solve_time = np.nan
            elif (level, function) == ("DEBUG", "_prepareLinearSystem"):
                if msg == "BEGIN _prepareLinearSystem":
                    begin_prepare_time = time_stamp
                    prepare_time = np.nan
                elif msg.startswith("END _prepareLinearSystem - "):
                    prepare_time = msg.split('-')[-1]
                    prepare_time = pd.to_timedelta(prepare_time)
                    begin_prepare_time = np.nan
                elif msg == "END _prepareLinearSystem":
                    prepare_time = time_stamp - begin_prepare_time
                    begin_prepare_time = np.nan

    return events

if __name__ == "__main__":
    # https://stackoverflow.com/a/55849527/2019542
    logger = logging.getLogger('rev_and_suite_permutations')
    fh = logging.FileHandler(str(snakemake.log[0]))
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    try:
        events = read_events(snakemake.input.log)

        df = pd.json_normalize(events)

        if "preconditioner" in df.columns:
            df.loc[df["preconditioner"].isna()
                   | (df["preconditioner"] == "NoneType"),
                   "preconditioner"] = "unpreconditioned"

        # Add configuration to each record, tidy-style
        config = pd.DataFrame([snakemake.params.config]).set_index("id", drop=False)
        df = config.join(df, lsuffix="_L", how="cross")

        df.to_json(snakemake.output[0], date_format="iso")
    except Exception as e:
        logger.error(e, exc_info=True)
        raise e
