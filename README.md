Timing and solver convergence info developed while working on #701 (Add 
support for PETSc solvers)


---------------------
DATA & FILE OVERVIEW
---------------------

```
./
├── NOTES.md
│    ⎧ Information about first round results
│    ⎩ Lost all context at this point?
├── README.md
│    ⎧ this file
│    ⎩ some text
├── TODO.md
│    { Tasks to work on
├── solver.1.log
├── solver.log
├── test_times.ipynb
├── results/
│   ├── diffusion/
│   │    { Simulation results from `{fipy}/examples/benchmarking/solvers/diffusion.py`.
│   │   ├── arithmetic_diffusion/
│   │   │    ⎧ `fast` queue on `ruth`
│   │   │    ⎨ diffusion coefficient equal to solution variable, evaluated as
│   │   │    ⎩ arithmetic mean.
│   │   ├── compare_1st/
│   │   │    ⎧ `fast` queue on `ruth` (?)
│   │   │    ⎨ diffusion coefficient equal to solution variable, evaluated as
│   │   │    ⎪ arithmetic mean.
│   │   │    ⎩ single sweep, with output of solution.
│   │   ├── constant_diffusion/
│   │   │    ⎧ `fast` queue on `ruth`
│   │   │    ⎩ constant diffusion coefficient.
│   │   ├── gpu/
│   │   │    ⎧ `rgpu5` (?) via `mr-french`
│   │   │    ⎨ diffusion coefficient equal to solution variable, evaluated as
│   │   │    ⎩ arithmetic mean.
│   │   ├── harmonic/
│   │   │    ⎧ `fast` queue on `ruth`
│   │   │    ⎨ diffusion coefficient equal to solution variable, evaluated as
│   │   │    ⎩ harmonic mean.
│   │   ├── harmonic_right/
│   │   │    ⎧ `fast` queue on `ruth`
│   │   │    ⎨ diffusion coefficient equal to solution variable, evaluated as
│   │   │    ⎪ harmonic mean.
│   │   │    ⎩ right-hand boundary condition set to ???.
│   │   ├── linux/
│   │   │    ⎧ - `fast` queue on `ruth`
│   │   │    ⎨ - diffusion coefficient equal to solution variable, evaluated as
│   │   │    ⎩ arithmetic mean.
│   │   ├── linux2/
│   │   │    ⎧ `fast` queue on `ruth`
│   │   │    ⎨ - diffusion coefficient equal to solution variable, evaluated as
│   │   │    ⎪ arithmetic mean.
│   │   │    ⎪ - Single sweep with no output.
│   │   │    ⎩ - @520c449b
│   │   └── macOS/
│   │        ⎧ - macOS `PN129671`
│   │        ⎨ - diffusion coefficient equal to solution variable, evaluated as
│   │        ⎩ arithmetic mean.
│   ├── nucleation/
│   │    { Simulation results using Phase Field Benchmark 8(b)
│   │   ├── nucleation1/
│   │   │    ⎧ - `ruth`
│   │   │    ⎨ - ```
│   │   │    ⎪ OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc mpirun -np 32 python examples/benchmarking/solvers/nucleation.py --numberOfElements=10000 --solver=gmres --preconditioner=jacobi --output=examples/benchmarking/solvers/nucleation1
│   │   │    ⎪ ```
│   │   │    ⎩ - job 357096
│   │   ├── nucleation10/
│   │   │    { - ???
│   │   ├── nucleation11/
│   │   │    { - ???
│   │   ├── nucleation12/
│   │   │    { - ???
│   │   ├── nucleation14/
│   │   │    ⎧ - `fast` queue on `ruth`
│   │   │    ⎨ - ```
│   │   │    ⎪ sbatch --time="12:00:00" --partition=fast --job-name="nucleation initialization" --ntasks=32 --ntasks-per-core=2 examples/benchmarking/solvers/setup.sh --env fipy310 -- OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc mpirun -np 32 python examples/benchmarking/solvers/nucleation.py --numberOfElements=1000000 --solver=gmres --preconditioner=jacobi --output=examples/benchmarking/solvers/nucleation14
│   │   │    ⎪ ```
│   │   │    ⎩ - job 357099
│   │   ├── nucleation16/
│   │   │    { - ???
│   │   ├── nucleation18/
│   │   │    ⎧ - `fast` queue on `ruth`
│   │   │    ⎨ - 600 time units.
│   │   │    ⎪ - ```
│   │   │    ⎪ sbatch --time="12:00:00" --partition=fast --job-name="nucleation initialization" --ntasks=32 --ntasks-per-core=2 examples/benchmarking/solvers/setup.sh --env fipy310 -- OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc mpirun -np 32 python examples/benchmarking/solvers/nucleation.py --numberOfElements=1000000 --solver=lu --output=examples/benchmarking/solvers/nucleation18
│   │   │    ⎩ ```
│   │   ├── nucleation19/
│   │   │    { - ???
│   │   ├── nucleation20/
│   │   │    { - ???
│   │   ├── nucleation21/
│   │   │    { - ???
│   │   ├── nucleation22/
│   │   │    ⎧ - `fast` queue on `ruth`
│   │   │    ⎨ - ```
│   │   │    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch fast --env fipy27 --solversuite pysparse --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation22/solver.log nucleation.py --preconditioner=^C-output /working/guyer/fipy/examples/benchmarking/solvers/nucleation22 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver
│   │   │    ⎩ ```
│   │   ├── nucleation23/
│   │   │    ⎧ - `fast` queue on `ruth`
│   │   │    ⎨ - Use `nucleation18/t=300.0.npz` as initial condition
│   │   │    ⎩ - run for 1 time unit.
│   │   ├── nucleation23_24/
│   │   │    ⎧ - Temporary merger of two different simulation runs.
│   │   │    ⎩ - No longer needed as `nucleation24/` results have been merged into `nucleation23/`?
│   │   ├── nucleation24/
│   │   ├── nucleation25/
│   │   │    ⎧ - `gpu` queue on `mr-french`
│   │   │    ⎨ - Use nucleation18/t=300.0.npz as initial condition
│   │   │    ⎩ - run for 1 time unit.
│   │   ├── nucleation3/
│   │   │    ⎧ - macOS `PN129671`
│   │   │    ⎨ - Simulation results using Phase Field Benchmark 8(b)
│   │   │    ⎪ - ```
│   │   │    ⎪ OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc python examples/benchmarking/solvers/nucleation.py --numberOfElements=10000 --solver=gmres --preconditioner=jacobi --output=examples/benchmarking/solvers/nucleation3
│   │   │    ⎩ ```
│   │   ├── nucleation4/
│   │   │    ⎧ - macOS `PN129671`
│   │   │    ⎨ - Simulation results using Phase Field Benchmark 8(b)
│   │   │    ⎪ - ```
│   │   │    ⎪ OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc python examples/benchmarking/solvers/nucleation.py --numberOfElements=10000 --solver=pcg --preconditioner=jacobi --output=examples/benchmarking/solvers/nucleation4
│   │   │    ⎩ ```
│   │   ├── nucleation5/
│   │   │    ⎧ - macOS `PN129671`
│   │   │    ⎨ - Simulation results using Phase Field Benchmark 8(b)
│   │   │    ⎪ - ```
│   │   │    ⎪ OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc python examples/benchmarking/solvers/nucleation.py --numberOfElements=10000 --output=examples/benchmarking/solvers/nucleation5
│   │   │    ⎩ ```
│   │   ├── nucleation6/
│   │   │    { - ???
│   │   ├── nucleation7/
│   │   │    { - ???
│   │   ├── nucleation8/
│   │   │    { - ???
│   │   └── nucleation9/
│   │        { - ???
│   └── old/
│        { old stuff
├── solver_diagnostics/
│   ├── CH2D_petsc.txt
│   ├── coupled.txt
│   ├── impingement_petsc.txt
│   ├── liquidVapor1D.txt
│   ├── liquidVapor1D_unscaled.txt
│   ├── petsc_Jacobi.txt
│   ├── petsc_Jacobi_LU_maxdiag.txt
│   ├── petsc_LU_maxdiag.txt
│   ├── petsc_default_precon.txt
│   ├── petsc_impingement.txt
│   ├── petsc_maxdiag.txt
│   ├── petsc_solve.txt
│   ├── pysparse_solve.txt
│   ├── tanh1D_petsc.txt
│   ├── tanh1D_trilinos.txt
│   ├── trilinos_impingement.txt
│   └── trilinos_solve.txt
└── timings/
    ├── petsc_1.tsv
    ├── petsc_2.tsv
    ├── petsc_3.tsv
    ├── petsc_3b.tsv
    ├── petsc_4.tsv
    ├── petsc_4b.tsv
    ├── pysparse_1.tsv
    ├── scipy_1.tsv
    ├── trilinos_1.tsv
    └── trilinos_2.tsv
```

------------------------------------
DATA-SPECIFIC INFORMATION FOR: *.log
------------------------------------

[Python event log](https://docs.python.org/3/library/logging.html).

-----------------------------------
DATA-SPECIFIC INFORMATION FOR: *.md
-----------------------------------

[Markdown](https://daringfireball.net/projects/markdown/)-formatted text.
