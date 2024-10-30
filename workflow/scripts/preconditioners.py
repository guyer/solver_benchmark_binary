"""Introspect FiPy for names of preconditioners in active suite
"""

import fipy as fp

preconditioners = []

for k, v in fp.solvers.__dict__.items():
    if isinstance(v, type) and issubclass(v, fp.preconditioner.Preconditioner):
        preconditioners.append(k)

print(" ".join(sorted(preconditioners)))
