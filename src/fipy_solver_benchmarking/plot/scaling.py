import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from matplotlib.legend import Legend


def plot_all(df, color_by_suite=True,
             by=["package.solver", "solver_class", "preconditioner"],
             solve_time=False, title=None):
    color_map = {
        'no-pysparse': 'red',
        'trilinos': 'red',
        'petsc': 'blue',
        'scipy': 'green',
        'pysparse': 'orange',
        'pyamgx': 'cyan'
    }
    
    fig, ax = plt.subplots(figsize=(8,6))
    groups = df.groupby(by + ["numberOfElements"])
    groups = groups.agg(converged=("converged", "all"),
                        elapsed_count=("elapsed_seconds", "count"),
                        elapsed_mean=("elapsed_seconds", "mean"),
                        elapsed_std=("elapsed_seconds", "std"),
                        solve_count=("solve_seconds", "count"),
                        solve_mean=("solve_seconds", "mean"),
                        solve_std=("solve_seconds", "std")).reset_index()
    groups = groups.groupby(by)
    for key, group in groups:
        if color_by_suite:
            color = color_map[key[0]]
            group.mask(~group["converged"].astype(bool)).plot("numberOfElements", "elapsed_mean", loglog=True,
                       ax=ax, label=key, color=color, marker=".", markersize=1.5)
            group.mask(group["converged"].astype(bool)).plot("numberOfElements", "elapsed_mean", loglog=True,
                       ax=ax, label=None, color=color, marker="x", linestyle="")
            
            if solve_time:
                group.mask(~group["converged"].astype(bool)).plot("numberOfElements", "solve_mean", loglog=True,
                           ax=ax, label=key, color=color, marker=".", markersize=1.5, linestyle="--")
                group.mask(group["converged"].astype(bool)).plot("numberOfElements", "solve_mean", loglog=True,
                           ax=ax, label=None, color=color, marker="+", linestyle="")
        else:
            group.mask(~group["converged"].astype(bool)).plot("numberOfElements", "elapsed_mean", loglog=True,
                       ax=ax, label=key, marker=".", markersize=1.5)
            color = ax.lines[-1].get_color()
            group.mask(group["converged"].astype(bool)).plot("numberOfElements", "elapsed_mean", loglog=True,
                       ax=ax, label=None, color=color, marker="x", linestyle="")
            
            if solve_time:
                group.mask(~group["converged"].astype(bool)).plot("numberOfElements", "solve_mean", loglog=True,
                           ax=ax, label=key, color=color, marker=".", markersize=1.5, linestyle="--")
                group.mask(group["converged"].astype(bool)).plot("numberOfElements", "solve_mean", loglog=True,
                           ax=ax, label=None, color=color, marker="+", linestyle="")

        # plot uncertainty
        err = group["elapsed_std"] / np.sqrt(group["elapsed_count"])
        ax.fill_between(group["numberOfElements"],
                        group["elapsed_mean"] - err,
                        group["elapsed_mean"] + err,
                        color=color,
                        alpha=0.1)

    if color_by_suite:
        legend_elements = [Line2D([0], [0], color=c, label=s)
                           for s, c in color_map.items()
                           if s in df["package.solver"].unique()]
        ax.legend(handles=legend_elements, bbox_to_anchor=(1., 1.))
    else:
        # only label converged lines
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::2], labels[::2], bbox_to_anchor=(1., 1.))

    if ~df["converged"].astype(bool).all():
        legend_elements = [Line2D([0], [0], color="black", marker=".", markersize=1.5),
                           Line2D([0], [0], color="black", linewidth=0, marker="x")]
        labels = ["converged", "not converged"]

        if solve_time:
            legend_elements.append(Line2D([0], [0], color="black", marker=".", markersize=1.5, linestyle="--"))
            labels.append("solve time")

        leg = Legend(ax, handles=legend_elements, labels=labels,
                     loc='upper left', frameon=False)
        ax.add_artist(leg);

    ax.set_ylabel("elapsed time / s")
    ax.set_xlabel("number of elements")

    if title is not None:
        ax.set_title(title)

    plt.show()
    
    return ax

def plot_solve_fraction(df, color_by_suite=True,
             by=["package.solver", "solver_class", "preconditioner"]):
    color_map = {
        'no-pysparse': 'red',
        'trilinos': 'red',
        'petsc': 'blue',
        'scipy': 'green',
        'pysparse': 'orange',
        'pyamgx': 'cyan'
    }
    
    fig, ax = plt.subplots(figsize=(8,6))
    df = df.copy()
    df["solve_fraction"] = df["solve_seconds"] / df["elapsed_seconds"]
    groups = df.groupby(by + ["numberOfElements"])
    groups = groups.agg(converged=("converged", "all"),
                        elapsed_count=("elapsed_seconds", "count"),
                        elapsed_mean=("elapsed_seconds", "mean"),
                        elapsed_std=("elapsed_seconds", "std"),
                        solve_count=("solve_seconds", "count"),
                        solve_mean=("solve_seconds", "mean"),
                        solve_std=("solve_seconds", "std"),
                        solve_fraction_count=("solve_fraction", "count"),
                        solve_fraction_mean=("solve_fraction", "mean"),
                        solve_fraction_std=("solve_fraction", "std")).reset_index()
    
    groups = groups.groupby(by)
    for key, group in groups:
        if color_by_suite:
            color = color_map[key[0]]
            group.mask(~group["converged"].astype(bool)).plot("numberOfElements", "solve_fraction_mean", logx=True,
                       ax=ax, label=key, color=color, marker=".", markersize=1.5)
            group.mask(group["converged"].astype(bool)).plot("numberOfElements", "solve_fraction_mean", logx=True,
                       ax=ax, label=None, color=color, marker="x", linestyle="")
        else:
            group.mask(~group["converged"].astype(bool)).plot("numberOfElements", "solve_fraction_mean", logx=True,
                       ax=ax, label=key, marker=".", markersize=1.5)
            color = ax.lines[-1].get_color()
            group.mask(group["converged"].astype(bool)).plot("numberOfElements", "solve_fraction_mean", logx=True,
                       ax=ax, label=None, color=color, marker="x", linestyle="")
            
        # plot uncertainty
        err = group["solve_fraction_std"] / np.sqrt(group["solve_fraction_count"])
        ax.fill_between(group["numberOfElements"],
                        group["solve_fraction_mean"] - err,
                        group["solve_fraction_mean"] + err,
                        color=color,
                        alpha=0.1)

    if color_by_suite:
        legend_elements = [Line2D([0], [0], color=c, label=s)
                           for s, c in color_map.items()
                           if s in df["package.solver"].unique()]
        ax.legend(handles=legend_elements, bbox_to_anchor=(1., 1.))
    else:
        # only label converged lines
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::2], labels[::2], bbox_to_anchor=(1., 1.))

    if ~df["converged"].astype(bool).all():
        legend_elements = [Line2D([0], [0], color="black", marker=".", markersize=1.5),
                           Line2D([0], [0], color="black", linewidth=0, marker="x")]
        leg = Legend(ax, handles=legend_elements, labels=["converged", "not converged"],
                     loc='upper left', frameon=False)
        ax.add_artist(leg);

    ax.set_ylabel("solve fraction")
    ax.set_xlabel("number of elements")
    
    plt.show()

    return ax

def plot_by_solver(df):
    for solver_class, group1 in df.groupby(["solver_class"]):
        plot_all(group1, by=["package.solver", "preconditioner"], title=solver_class)
        
def plot_sweep_times(df):
    for numberOfElements, group1 in df.groupby("numberOfElements"):
        fig, ax = plt.subplots(figsize=(8,6))
        for label, group2 in group1.groupby(["package.solver", "solver_class", "preconditioner"]):
            group2 = group2.reset_index()
            group2.plot(y="elapsed_seconds", ax=ax, label=label, logy=True)
        _ = plt.xlim(xmax=9)
        plt.legend(loc="upper right", bbox_to_anchor=(1.7, 1.0))
        plt.title("size = {0}x{0}".format(int(np.sqrt(numberOfElements))))
        plt.ylabel("sweep time / s")
        plt.xlabel("sweep")  

        plt.show()
