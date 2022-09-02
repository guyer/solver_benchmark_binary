Timing and solver convergence info developed while working on #701 (Add 
support for PETSc solvers)


---------------------
DATA & FILE OVERVIEW
---------------------

```
results/diffusion/
├── README.md
│    { this file
├── arithmetic_diffusion
│    ⎧ `fast` queue on `ruth`
│    ⎨ diffusion coefficient equal to solution variable, evaluated as
│    ⎩ arithmetic mean.
│   ├── solver.NNNNNN.log
│   ├── solver.MMMMMM.log
│   ├── ...
│   ├── petsc
│   │   ├── LinearCGSSolver
│   │   │   ├── ILUPreconditioner
│   │   │   │   ├── 100
│   │   │   │   │   └── solution.npz
│   │   │   │   ├── 10000
│   │   │   │   │   └── ...
│   │   │   │   ├── 1000000
│   │   │   │   │   └── ...
│   │   │   │   ├── 9
│   │   │   │   │   └── ...
│   │   │   │   ├── 961
│   │   │   │   │   └── ...
│   │   │   │   └── 99856
│   │   │   │   │   └── ...
│   │   │   ├── JacobiPreconditioner
│   │   │   │   └── ...
│   │   │   └── NoneType
│   │   │       └── ...
│   │   ├── LinearGMRESSolver
│   │   │   └── ...
│   │   ├── LinearLUSolver
│   │   │   └── ...
│   │   └── LinearPCGSolver
│   │       └── ...
│   ├── pysparse
│   │   └── ...
│   ├── scipy
│   │   └── ...
│   └── trilinos
│       └── ...
├── compare_1st
│    ⎧ `fast` queue on `ruth` (?)
│    ⎨ diffusion coefficient equal to solution variable, evaluated as
│    ⎪ arithmetic mean.
│    ⎩ single sweep, with output of solution.
├── constant_diffusion
│    ⎧ `fast` queue on `ruth`
│    ⎩ constant diffusion coefficient.
├── gpu
│    ⎧ `rgpu5` (?) via `mr-french`
│    ⎨ diffusion coefficient equal to solution variable, evaluated as
│    ⎩ arithmetic mean.
├── harmonic
│    ⎧ `fast` queue on `ruth`
│    ⎨ diffusion coefficient equal to solution variable, evaluated as
│    ⎩ harmonic mean.
├── harmonic_right
│    ⎧ `fast` queue on `ruth`
│    ⎨ diffusion coefficient equal to solution variable, evaluated as
│    ⎪ harmonic mean.
│    ⎩ right-hand boundary condition set to 0.1.
├── linux
│    ⎧ - `fast` queue on `ruth`
│    ⎨ - diffusion coefficient equal to solution variable, evaluated as
│    ⎩ arithmetic mean.
├── linux2
│    ⎧ `fast` queue on `ruth`
│    ⎨ - diffusion coefficient equal to solution variable, evaluated as
│    ⎪ arithmetic mean.
│    ⎪ - Single sweep with no output.
│    ⎩ - @520c449b
└── macOS
     ⎧ - macOS `PN129671`
     ⎨ - diffusion coefficient equal to solution variable, evaluated as
     ⎩ arithmetic mean.
```