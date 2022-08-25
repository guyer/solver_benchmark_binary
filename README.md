Timing and solver convergence info developed while working on #701 (Add 
support for PETSc solvers)


---------------------
DATA & FILE OVERVIEW
---------------------

- NOTES.md
    - Information about first round results
    - Lost all context at this point?
- README.md
    - This file
- TODO.md
    - Tasks to work on
- arithmetic_diffusion/
    - `fast` queue on `ruth` 
    - `{fipy}/examples/benchmarking/solvers/diffusion.py`
    - diffusion coefficient equal to solution variable, evaluated as
    arithmetic mean.
- compare_1st/
    - `fast` queue on `ruth` (?)
    - `{fipy}/examples/benchmarking/solvers/diffusion.py`
    - diffusion coefficient equal to solution variable, evaluated as
    arithmetic mean.
    - single sweep, with output of solution.
- constant_diffusion
    - `fast` queue on `ruth`
    - `{fipy}/examples/benchmarking/solvers/diffusion.py`
    - constant diffusion coefficient.
- gpu/
    - `rgpu5` (?) via `mr-french`
    - `{fipy}/examples/benchmarking/solvers/diffusion.py`
    - diffusion coefficient equal to solution variable, evaluated as
    arithmetic mean.
- harmonic/
    - `fast` queue on `ruth`
    - `{fipy}/examples/benchmarking/solvers/diffusion.py`
    - diffusion coefficient equal to solution variable, evaluated as
    harmonic mean.
- harmonic_right/
    - `fast` queue on `ruth`
    - `{fipy}/examples/benchmarking/solvers/diffusion.py`
    - diffusion coefficient equal to solution variable, evaluated as
    harmonic mean.
    - right-hand boundary condition set to ???.
- linux/
    - `fast` queue on `ruth`
    - `{fipy}/examples/benchmarking/solvers/diffusion.py`
    - diffusion coefficient equal to solution variable, evaluated as
    arithmetic mean.
- linux2/
    - `fast` queue on `ruth`
    - `{fipy}/examples/benchmarking/solvers/diffusion.py`
    - diffusion coefficient equal to solution variable, evaluated as
    arithmetic mean.
    - Single sweep with no output.
    - @520c449b
- macOS/
    - macOS `PN129671`
    - `{fipy}/examples/benchmarking/solvers/diffusion.py`
    - diffusion coefficient equal to solution variable, evaluated as
    arithmetic mean.
- nucleation1/
    - `ruth`
    - Simulation results using Phase Field Benchmark 8(b)
    - ```
    OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc mpirun -np 32 python examples/benchmarking/solvers/nucleation.py --numberOfElements=10000 --solver=gmres --preconditioner=jacobi --output=examples/benchmarking/solvers/nucleation1
    ```
    - job 357096
- nucleation2/
    - `ruth`
    - Simulation results using Phase Field Benchmark 8(b)
    - ```
    OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc python examples/benchmarking/solvers/nucleation.py --numberOfElements=10000 --solver=gmres --preconditioner=jacobi --output=examples/benchmarking/solvers/nucleation2
    ```
    - job 357097
- nucleation3/
    - macOS `PN129671`
    - Simulation results using Phase Field Benchmark 8(b)
    - ```
    OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc python examples/benchmarking/solvers/nucleation.py --numberOfElements=10000 --solver=gmres --preconditioner=jacobi --output=examples/benchmarking/solvers/nucleation3
    ```
- nucleation4/
    - macOS `PN129671`
    - Simulation results using Phase Field Benchmark 8(b)
    - ```
    OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc python examples/benchmarking/solvers/nucleation.py --numberOfElements=10000 --solver=pcg --preconditioner=jacobi --output=examples/benchmarking/solvers/nucleation4
    ```
- nucleation5/
    - macOS `PN129671`
    - Simulation results using Phase Field Benchmark 8(b)
    - ```
    OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc python examples/benchmarking/solvers/nucleation.py --numberOfElements=10000 --output=examples/benchmarking/solvers/nucleation5
    ```
- nucleation6/
    - ???
- nucleation7/
    - ???
- nucleation8/
    - ???
- nucleation9/
    - ???
- nucleation10/
    - ???
- nucleation11/
    - ???
- nucleation12/
    - ???
- nucleation13/
    - `fast` queue on `ruth`
    - ```
    sbatch --time="12:00:00" --partition=fast --job-name="nucleation initialization" examples/benchmarking/solvers/setup.sh --env fipy310 -- OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc python examples/benchmarking/solvers/nucleation.py --numberOfElements=10000 --solver=gmres --preconditioner=jacobi --output=examples/benchmarking/solvers/nucleation13
    ```
    - job 357098
    - Never rsynced from 
    `ruth:/working/guyer/fipy/examples/benchmarking/solvers/`.
- nucleation14/
    - `fast` queue on `ruth`
    - ```
    sbatch --time="12:00:00" --partition=fast --job-name="nucleation initialization" --ntasks=32 --ntasks-per-core=2 examples/benchmarking/solvers/setup.sh --env fipy310 -- OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc mpirun -np 32 python examples/benchmarking/solvers/nucleation.py --numberOfElements=1000000 --solver=gmres --preconditioner=jacobi --output=examples/benchmarking/solvers/nucleation14
    ```
    - job 357099
- nucleation15/
    - `fast` queue on `ruth`
    - ```
    sbatch --time="12:00:00" --partition=fast --job-name="nucleation initialization" --ntasks=32 --ntasks-per-core=2 examples/benchmarking/solvers/setup.sh --env fipy310 -- OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc mpirun -np 32 python examples/benchmarking/solvers/nucleation.py --numberOfElements=1000000 --solver=lu --output=examples/benchmarking/solvers/nucleation15
    ```
    - Never rsynced from 
    `ruth:/working/guyer/fipy/examples/benchmarking/solvers/`.
- nucleation16/
    - ???
- nucleation17/
    - Precursor to `nucleation18/`.
    - Never rsynced from 
    `ruth:/working/guyer/fipy/examples/benchmarking/solvers/`.
- nucleation18/
    - `fast` queue on `ruth`
    - Phase Field Benchmark 8(b) for 600 time units.
    - ```
    sbatch --time="12:00:00" --partition=fast --job-name="nucleation initialization" --ntasks=32 --ntasks-per-core=2 examples/benchmarking/solvers/setup.sh --env fipy310 -- OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc mpirun -np 32 python examples/benchmarking/solvers/nucleation.py --numberOfElements=1000000 --solver=lu --output=examples/benchmarking/solvers/nucleation18
    ```
- nucleation19/
    - ???
- nucleation20/
    - ???
- nucleation21/
    - ???
- nucleation22/
    - `fast` queue on `ruth`
    - ```
    bash examples/benchmarking/solvers/dispatch.sh --sbatch fast --env fipy27 --solversuite pysparse --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation22/solver.log nucleation.py --preconditioner=^C-output /working/guyer/fipy/examples/benchmarking/solvers/nucleation22 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t\=300.0.npz --totaltime=301 --store_by_solver
    ```
- nucleation23/
    - `fast` queue on `ruth`
    - Phase Field Benchmark 8(b) with `nucleation18/t=300.0.npz` as initial
    condition, run for 1 time unit.
- nucleation23_24/
    - Temporary merger of two different simulation runs.
    - No longer needed as `nucleation24/` results have been merged into
    `nucleation23/`?
- nucleation25/
    - `gpu` queue on `mr-french`
    - Phase Field Benchmark 8(b) with nucleation18/t=300.0.npz as initial
    condition, run for 1 time unit.
- solver.1.log
- solver.log
- solver_diagnostics/
- test_solvers/
- test_solvers2/
- test_times.ipynb
- timings/


------------------------------------
DATA-SPECIFIC INFORMATION FOR: *.log
------------------------------------

[Python event log](https://docs.python.org/3/library/logging.html).

-----------------------------------
DATA-SPECIFIC INFORMATION FOR: *.md
-----------------------------------

[Markdown](https://daringfireball.net/projects/markdown/)-formatted text.
