"""Introspect FiPy for names of solvers and preconditioners in active suite
"""

import logging
import os
import json

if "FIPY_SOLVERS" not in os.environ:
    os.environ["FIPY_SOLVERS"] = snakemake.wildcards.suite

import fipy as fp

def get_classes(cls):
    classes = []

    for k, v in fp.solvers.__dict__.items():
        if (isinstance(v, type)
            and issubclass(v, cls)):
                classes.append(k)

    return classes

if __name__ == "__main__":
    # https://stackoverflow.com/a/55849527/2019542
    logger = logging.getLogger('get_solvers_preconditioners')
    fh = logging.FileHandler(str(snakemake.log[0]))
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    try:
        solvers = get_classes(cls=fp.solver.Solver)
        exclude = ["Solver", "DefaultSolver", "DefaultAsymmetricSolver",
                   "DummySolver", "GeneralSolver"]
        solvers = [solver for solver in solvers if solver not in exclude]
        
        solvers = ["LinearGMRESSolver"]

        with open(snakemake.output["solvers"], 'w') as f:
            json.dump({"solver": solvers}, f)

        preconditioners = get_classes(cls=fp.preconditioner.Preconditioner)
        preconditioners.append("none")

        preconditioners = ["JacobiPreconditioner", "LUPreconditioner"]

        with open(snakemake.output["preconditioners"], 'w') as f:
            json.dump({"preconditioner": preconditioners}, f)
    except Exception as e:
        logger.error(e, exc_info=True)
        raise e
