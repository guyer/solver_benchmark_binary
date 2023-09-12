import os


def read_tsv(tsv):
    x, y, v = np.loadtxt(tsv, skiprows=1, unpack=True)
    N = int(np.sqrt(len(x)))
    return [field.reshape((N,N)) for field in (x, y, v)] + [N]

def read_npz(npz):
    data = np.load(npz)
    return data["x"], data["y"], data["value"], data["value"].shape[0]

def tsv2npz(tsv):
    x, y, v, N = read_tsv(tsv)
    np.savez(os.path.splitext(tsv)[0], x=x, y=y, value=v)

def load_field(record, fname, key):
    """Load and intorpolate `path2` to resolution of `path1`
    
    Parameters
    ----------
    record : pd.Series
        Pandas Series with information about particular results to load.
    fname : str
        Name of `.npz` file containing x-coordinates, y-coordinates, and values
        to load from directory hierarchy determined by "logfile",
        "package.solver", "solver_class", "preconditioner", "numberOfElements"
        fields in `record`.
    key : str
        Designation of the field of interest within `fname`.

    Returns
    -------
    field : ndarray
        2D array of values.
    """
    preconditioner = record["preconditioner"]
    if preconditioner == "unpreconditioned":
        preconditioner = "NoneType"

    path = os.path.join(os.path.dirname(record["logfile"]),
                        record["package.solver"], record["solver_class"],
                        preconditioner, str(record["numberOfElements"]),
                        fname)

    try:
        with np.load(path) as data:
            field = data[key]
    except:
        field = None

    return field
