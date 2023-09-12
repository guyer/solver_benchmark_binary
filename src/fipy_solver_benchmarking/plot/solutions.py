import matplotlib as mpl
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage

from ..io.data import load_field


def scale_and_diff(v1, v2):
    """Interpolate and calculate difference between `v1` and `v2`
    
    Parameters
    ----------
    v1, v2 : ndarray
        Square arrays of values.

    Returns
    -------
    diff : ndarray
        Difference between `v1` and `v2`
        (at resolution of whichever is lower).
    """
    if v1 is None or v2 is None:
        return None

    N1 = v1.shape[0]
    N2 = v2.shape[0]

    if N1 == N2:
        diff = v2 - v1
    elif N2 > N1:
        diff = ndimage.zoom(v2, N1 / N2) - v1
    else:
        diff = v2 - ndimage.zoom(v1, N2 / N1)
    
    return abs(diff)

def extract_and_diff(v1, v2):
    """Slice and calculate difference between `v1` and `v2`
    
    Parameters
    ----------
    v1, v2 : ndarray
        Square arrays of values.

    Returns
    -------
    diff : ndarray
        Difference between `v1` and `v2`
        (at size of whichever is smaller).
    """
    if v1 is None or v2 is None:
        return None

    N1 = v1.shape[0]
    N2 = v2.shape[0]

    if N1 == N2:
        diff = v2 - v1
    elif N2 > N1:
        lower, upper = int((N2 - N1) / 2), int((N2 + N1) / 2)
        diff = v2[lower:upper, lower:upper] - v1
    else:
        lower, upper = int((N1 - N2) / 2), int((N1 + N2) / 2)
        diff = v2 - v1[lower:upper, lower:upper]
    
    return abs(diff)

def value1(v1, v2):
    """Extract part of `v1` corresponding to `v2`.
    
    Parameters
    ----------
    v1, v2 : ndarray
        Square arrays of values.

    Returns
    -------
    value : ndarray
        Slice of `v1` that overlaps with `v2`.
    """
    if v1 is None or v2 is None:
        return None

    N1 = v1.shape[0]
    N2 = v2.shape[0]

    if N1 <= N2:
        value = v1
    else:
        lower, upper = (N1 - N2) // 2, (N1 + N2) // 2
        value = v1[lower:upper, lower:upper]
    
    return value

def value2(v1, v2):
    """Extract part of `v2` corresponding to `v1`.
    
    Parameters
    ----------
    v1, v2 : ndarray
        Square arrays of values.

    Returns
    -------
    value : ndarray
        Slice of `v1` that overlaps with `v2`.
    """
    if v1 is None or v2 is None:
        return None

    N1 = v1.shape[0]
    N2 = v2.shape[0]

    if N1 >= N2:
        value = v2
    else:
        lower, upper = int((N2 - N1) / 2), int((N2 + N1) / 2)
        value = v2[lower:upper, lower:upper]
    
    return value

def perpendicular_label(ax, baseline, color, label, rows):
    """Annotate text perpendicular to the outer edge of axes.
    
    Parameters
    ----------
    ax : ~matploblib.axes.Axes
        Where to apply label.
    baseline : float
        The distance from the axes to position the label.
    color : color
        The foreground color of the label.
    label : str
        The text of the label.
    rows : bool
        Whether to apply label to rows or columns.

    Returns
    -------
    ann : ~matplotlib.text.Annotation
        The displayed label.
    """
    kwargs = dict(text=label, size=8, color=color,
                  textcoords='offset points', xycoords='axes fraction')
    if rows:
        ann = ax.annotate(xy=(0.5, 1), xytext=(0, baseline), rotation=90,
                          ha='center', va='bottom', **kwargs)
    else:
        ann = ax.annotate(xy=(0, 0.5), xytext=(-baseline, 0), rotation=0,
                          ha='right', va='center', **kwargs)
    return ann

def parallel_label(ax, baseline, color, label, rows):
    """Annotate text parallel to the outer edge of axes.
    
    Parameters
    ----------
    ax : ~matploblib.axes.Axes
        Where to apply label.
    baseline : float
        The distance from the axes to position the label.
    color : color
        The foreground color of the label.
    label : str
        The text of the label.
    rows : bool
        Whether to apply label to rows or columns.

    Returns
    -------
    ann : ~matplotlib.text.Annotation
        The displayed label.
    """
    kwargs = dict(text=label, size=10, color=color,
                  textcoords='offset points', xycoords='axes fraction',
                  ha='center', va='center')
    if rows:
        ann = ax.annotate(xy=(0.5, 1), xytext=(0, baseline), rotation=0,
                          **kwargs)
    else:
        ann = ax.annotate(xy=(0, 0.5), xytext=(-baseline, 0), rotation=90,
                          **kwargs)
    return ann

def all_converged(df):
    """Determine whether all elements are converged.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Pandas object with a boolean "converged" field.

    Returns
    -------
    bool
        Whether all elements of `df` are converged.
    """
    try:
        # df is a Series
        converged = bool(df.all())
    except ValueError:
        # df is a DataFrame
        converged = bool(df.all().all())
        
    return converged

def label_width(ann, rows):
    """Calculate display width of label annotation.
    
    Parameters
    ----------
    ann : ~matplotlib.text.Annotation
        The displayed label.
    rows : bool
        Whether label was applied to rows or columns.
    
    Returns
    -------
    float
        Width of :class:`~matplotlib.transforms.Bbox` of label
        perpendicular to rows (columns).
    """
    renderer = ann.figure.canvas.get_renderer()
    bbox = ann.get_window_extent(renderer=renderer)
    if rows:
        width = bbox.y1 - bbox.y0
    else:
        width = bbox.x1 - bbox.x0

    return width

def annotate_label(df, ax, name, label, rows, baseline, perpendicular):
    """Apply label to subplot.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Pandas object with a boolean "converged" field.
    ax : ~matploblib.axes.Axes
        Where to display difference image.
    name : str
        Name for this level of index labels.
    label
        Value of index
    rows : bool
        Whether to apply labels to rows or columns.
    baseline : float
        The distance from the axes to position the label.
    perpendicular : bool
        Whether label should be perpendicular to rows (columns) or parallel.

    Returns
    -------
    ann : ~matplotlib.text.Annotation
        The displayed label.
    width : float
        Displayed width of :class:`~matplotlib.transforms.Bbox` of label
        perpendicular to rows (columns).
    """
    if name == "numberOfElements":
        label = '{0}x{0}'.format(int(np.sqrt(label)))

    if all_converged(df):
        color = "black"
    else:
        color = "red"

    if perpendicular:
        ann = perpendicular_label(ax, baseline, color, label, rows=rows)
    else:
        ann = parallel_label(ax, baseline, color, label, rows=rows)

    return ann, label_width(ann, rows)
    
def label_MultiIndex(df, axs, offset=0, rows=True, baseline=0):
    """Recursively label subplots according to a multilevel index.
    
    The labels closest to the :class:`~matploblib.axes.Axes` will be
    rotated perpendicular to the rows (columns). All remaining labels
    will be parallel.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Aggregated data with boolean "converged" field.
    axs : array-like (1-dimensional) of ~matploblib.axes.Axes
        Where to apply the labels.
    offset : int, default 0
        Initial row or column for labeling.
    rows : bool, default True
        Whether to apply labels to rows or columns.
    baseline : float, default 0
        The distance from the axes to position the labels.

    Returns
    -------
    perpendicular : bool
        Whether next level of labels should be perpendicular to rows (columns).
    count : int
        Number of rows or columns that were labeled.
    next_baseline : float
        The distance from the axes to position the next level of labels.
    """
    if isinstance(df.index, pd.MultiIndex):
        labels = df.index.levels[0]
        name = df.index.names[0]
    elif isinstance(df, pd.DataFrame):
        labels = df.index
        name = df.index.name
    else:
        # At lowest level of indices. Break recursion.
        return True, 1, baseline

    count = 0
    next_baseline = baseline
    for label in labels:
        try:
            labeldf = df.loc[label]
        except KeyError:
            continue

        (perpendicular,
         num,
         prev_baseline) = label_MultiIndex(df=labeldf,
                                           axs=axs, offset=offset + count,
                                           rows=rows, baseline=baseline)
        
        prev_baseline += 10
        ann, width = annotate_label(df=labeldf,
                                    ax=axs[offset + count + num // 2],
                                    name=name, label=label,
                                    rows=rows, baseline=prev_baseline + 10,
                                    perpendicular=perpendicular)
            
        next_baseline = max(next_baseline, prev_baseline + width)
            
        count += num

    return False, count, next_baseline

def plot_difference(ax, groupA, groupB, fname, key, diff_fn,
                    norm, converged_cmap, diverged_cmap):
    """Plot image showing error between two records.
    
    Parameters
    ----------
    ax : ~matploblib.axes.Axes
        Where to display difference image.
    groupA, groupB : pandas.DataFrame
        Records to compare. There may be multiple runs with the same parameters.
    fname : str
        Name of `.npz` file containing x-coordinates, y-coordinates, and values
        to load from directory hierarchy determined by "logfile",
        "package.solver", "solver_class", "preconditioner", "numberOfElements"
        fields in `dfA` and `dfB`.
    key : str
        Designation of the field of interest within `fname`.
    diff_fn : function
        Function to calculate difference between two ndarray objects.
    norm : 
        ???
    converged_cmap, diverged_cmap : 
    """
    fieldA = load_field(record=groupA.iloc[0], fname=fname, key=key)
    fieldB = load_field(record=groupB.iloc[0], fname=fname, key=key)
    diff = diff_fn(fieldA, fieldB)

    if diff is not None:
        if groupA["converged"].all() & groupB["converged"].all():
            cmap = converged_cmap
        else:
            cmap = diverged_cmap

        ax.imshow(diff, norm=norm, cmap=cmap)

def plot_error_matrix(dfA, dfB, by, fname, key="value",
                      diff_fn=scale_and_diff,
                      vmin=1e-16, vmax=1, logscale=True):
    """Plot matrix of images showing the error between two sets of records.
    
    Top and left axes will be labeled hierarchically by groups.
    
    Parameters
    ----------
    dfA, dfB : pandas.DataFrame
        Records to compare.
    by : list of str
        Keys to group `dfA` and `dfB` by.
    fname : str
        Name of `.npz` file containing x-coordinates, y-coordinates, and values
        to load from directory hierarchy determined by "logfile",
        "package.solver", "solver_class", "preconditioner", "numberOfElements"
        fields in `dfA` and `dfB`.
    key : str
        Designation of the field of interest within `fname`.
    diff_fn : function
        Function to calculate difference between two ndarray objects.
    vmin, vmax : float
        Display range of data.
    logscale : bool, default True
        Whether to display difference on logarithmic or linear scale.
    """
    groupsA = dfA.groupby(by)
    groupsB = dfB.groupby(by)
    
    fig, axs = plt.subplots(ncols=len(groupsA), nrows=len(groupsB),
                            figsize=(10, 10 * len(groupsB) / len(groupsA)))
    plt.setp(axs.flat, xticks=[], yticks=[], frame_on=False)
    
    if logscale:
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
    
    converged_cmap = mpl.colormaps['viridis'].copy()
    converged_cmap.set_bad(color='black')
    cax = plt.axes([0.95, 0.3, 0.025, 0.4])
    converged_cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=converged_cmap),
                                cax=cax)
    converged_cb.set_label("converged", labelpad=-70)

    diverged_cmap = mpl.colormaps['inferno'].copy()
    diverged_cmap.set_bad(color='black')
    cax = plt.axes([1.025, 0.3, 0.025, 0.4])
    diverged_cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=diverged_cmap),
                               cax=cax, ticklocation='left', format="")
    diverged_cb.set_label("diverged", labelpad=-45, color="red")
    
    for col, (keyA, groupA) in zip(axs.T, groupsA):
        for ax, (keyB, groupB) in zip(col, groupsB):
            plot_difference(ax=ax, groupA=groupA, groupB=groupB,
                            fname=fname, key=key,
                            diff_fn=diff_fn, norm=norm,
                            converged_cmap=converged_cmap,
                            diverged_cmap=diverged_cmap)

    label_MultiIndex(df=groupsA.agg({"converged": "all"}),
                     axs=axs[0], rows=True)
    
    label_MultiIndex(df=groupsB.agg({"converged": "all"}),
                     axs=axs[:, 0], rows=False)
    
    plt.show()
