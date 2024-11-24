"""Introspect FiPy for names of solvers in active suite
"""

import os

if "FIPY_SOLVERS" not in os.environ:
    os.environ["FIPY_SOLVERS"] = snakemake.wildcards.suite

import fipy as fp

if __name__ == "__main__":
    # https://stackoverflow.com/a/55849527/2019542
    logger = logging.getLogger('list_solvers')
    fh = logging.FileHandler(str(snakemake.log[0]))
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    try:
        solvers = []

        for k, v in fp.solvers.__dict__.items():
            if k not in ["Solver", "DefaultSolver", "DefaultAsymmetricSolver", "DummySolver", "GeneralSolver"]:
                if isinstance(v, type) and issubclass(v, fp.solver.Solver):
                    solvers.append(k)

        with open(snakemake.output[0], 'w') as f:
            f.write(" ".join(sorted(solvers)))
    except Exception as e:
        logger.error(e, exc_info=True)
        raise e
