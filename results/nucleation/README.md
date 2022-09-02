Timing and solver convergence info developed while working on #701 (Add 
support for PETSc solvers)


---------------------
DATA & FILE OVERVIEW
---------------------

```
results/nucleation/
├── README.md
│    { Description of `results/nucleation/`.
├── nucleation18/
│    ⎧ - `fast` queue on `ruth`
│    ⎨ - 600 time units.
│    ⎪ - ```
│    ⎪ sbatch --time="12:00:00" --partition=fast --job-name="nucleation initialization" --ntasks=32 --ntasks-per-core=2 examples/benchmarking/solvers/setup.sh --env fipy310 -- OMP_NUM_THREADS=1 FIPY_SOLVERS=petsc mpirun -np 32 python examples/benchmarking/solvers/nucleation.py --numberOfElements=1000000 --solver=lu --output=examples/benchmarking/solvers/nucleation18
│    ⎪ ```
│    ⎩ - @a8881302 - Fix typo
├── nucleation23/
│    ⎧ - `fast` queue on `ruth`
│    ⎨ - Use `nucleation18/t=300.0.npz` as initial condition
│    ⎪ - run for 1 time unit.
│    ⎪ - ```
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch fast "12:00:00" --env fipy27 --solversuite pysparse --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation23/solver.log nucleation.py --preconditioner=none --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation23 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch fast "12:00:00" --env fipy27 --solversuite pysparse --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation23/solver.log nucleation.py --preconditioner=jacobi --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation23 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch fast "12:00:00" --env fipy37 --solversuite trilinos --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation23/solver.log nucleation.py --preconditioner=jacobi --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation23 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch fast "12:00:00" --env fipy37 --solversuite trilinos --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation23/solver.log nucleation.py --preconditioner=ilu --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation23 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch fast "12:00:00" --env fipy37 --solversuite trilinos --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation23/solver.log nucleation.py --preconditioner=none --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation23 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch fast "12:00:00" --env fipy310 --solversuite petsc --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation23/solver.log nucleation.py --preconditioner=none --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation23 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch fast "12:00:00" --env fipy310 --solversuite petsc --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation23/solver.log nucleation.py --preconditioner=jacobi --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation23 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch fast "12:00:00" --env fipy310 --solversuite petsc --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation23/solver.log nucleation.py --preconditioner=ilu --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation23 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch fast "12:00:00" --env fipy310 --solversuite scipy --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation23/solver.log nucleation.py --preconditioner=ilu --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation23 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch fast "12:00:00" --env fipy310 --solversuite scipy --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation23/solver.log nucleation.py --preconditioner=jacobi --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation23 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch fast "12:00:00" --env fipy310 --solversuite scipy --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation23/solver.log nucleation.py --preconditioner=none --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation23 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ ```
│    ⎩ - @d60f8560 - Add slurm timeout
├── nucleation25/
│    ⎧ - `gpu` queue on `mr-french`
│    ⎨ - Use nucleation18/t=300.0.npz as initial condition
│    ⎪ - run for 1 time unit.
│    ⎪ - ```
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch gpu "12:00:00" --env fipy310 --solversuite pyamgx --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation25/solver.log nucleation.py --preconditioner=none --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation25 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch gpu "12:00:00" --env fipy310 --solversuite pyamgx --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation25/solver.log nucleation.py --preconditioner=jacobi --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation25 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch gpu "12:00:00" --env fipy310 --solversuite pyamgx --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation25/solver.log nucleation.py --preconditioner=ilu --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation25 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
│    ⎪ ```
│    ⎩ - @70f844e0 - Make pyamgx preconditioners proper subclasses
└── nucleation26/
     ⎧ - `gpu` queue on `mr-french`
     ⎨ - Use nucleation18/t=300.0.npz as initial condition
     ⎪ - run for 1 time unit.
     ⎪ - ```
     ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch gpu "12:00:00" --env fipy310 --solversuite pyamgx --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation26/solver.log nucleation.py --preconditioner=ilu --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation26 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
     ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch gpu "12:00:00" --env fipy310 --solversuite pyamgx --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation26/solver.log nucleation.py --preconditioner=jacobi --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation26 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
     ⎪ bash examples/benchmarking/solvers/dispatch.sh --sbatch gpu "12:00:00" --env fipy310 --solversuite pyamgx --log examples/benchmarking/solvers/config_template.json /working/guyer/fipy/examples/benchmarking/solvers/nucleation26/solver.log nucleation.py --preconditioner=none --output /working/guyer/fipy/examples/benchmarking/solvers/nucleation26 --restart /working/guyer/fipy/examples/benchmarking/solvers/nucleation18/t=300.0.npz --totaltime=301 --store_by_solver --checkpoint_interval=1.
     ⎪ ```
     ⎩ - repeat of `nucleation25`
```

(generated with `tree --info --filesfirst -FL 1 results/nucleation`)
