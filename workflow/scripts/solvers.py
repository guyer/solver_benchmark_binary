"""Introspect FiPy for names of solvers in active suite
"""

import fipy as fp

solvers = []

for k, v in fp.solvers.__dict__.items():
    if k not in ["Solver", "DefaultSolver", "DefaultAsymmetricSolver", "DummySolver", "GeneralSolver"]:
        if isinstance(v, type) and issubclass(v, fp.solver.Solver):
            solvers.append(k)

print(" ".join(sorted(solvers)))
