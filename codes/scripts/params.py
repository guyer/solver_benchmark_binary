"""Introspect FiPy for names of preconditioners in active suite
"""

import fipy as fp
import pandas as pd
from itertools import product

solvers = []

for k, v in fp.solvers.__dict__.items():
    if k not in ["Solver", "DefaultSolver", "DefaultAsymmetricSolver", "DummySolver", "GeneralSolver"]:
        if isinstance(v, type) and issubclass(v, fp.solver.Solver):
            solvers.append(k)

preconditioners = []

for k, v in fp.solvers.__dict__.items():
    if isinstance(v, type) and issubclass(v, fp.preconditioner.Preconditioner):
        preconditioners.append(k)

df = pd.DataFrame(data=list(product(solvers, preconditioners)),
                  columns=["solver", "preconditioner"])

print(df.to_csv(index=False))