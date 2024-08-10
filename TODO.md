- [ ] Reimpliment using signac?
- [ ] Reimpliment using SnakeMake?
- [ ] Reimpliment using Sumatra?
- [ ] Reimpliment using CORR?

- [x] split out history from {fipy}/examples/benchmarking/solvers and merge with this
- [ ] Address suggestions from {fipy}/examples/benchmarking/solvers/2022-08-31.md

- [x] Make a common `Solver._solve()` with steps implemented by different 
      suites, instead of unique (but not really) `_solve()` for each.
- [?] Make all solvers handle units.
- [ ] Figure out matrix and RHS caching
- [x] What is `Solver._Lxb` about and why only used by `Solver._norms`?
      Appropriate for `Solver._solve()`?
- [x] Document `Solver._solve()`.
- [x] Document `Solver._solve_()`.
- [x] Switch to `criterion="default"`
    - [x] examples/cahnHilliard/mesh2D.py
    - [x] examples/cahnHilliard/tanh1D.py
    - [x] examples/convection/exponential1DSource/tri2D.py
    - [x] examples/diffusion/nthOrder/input4thOrder1D.py
    - [x] examples/elphf/diffusion/mesh1D.py
    - [x] examples/flow/stokesCavity.py
    - [x] examples/phase/binary.py
    - [x] examples/reactiveWetting/liquidVapor1D.py
        - `initialResidual` is 0.03338815978162811
        - sweep residual is 0.03338815978162811 in perpetuity
        - with `legacy`, second sweep residual jumps to ~351000000
          and adaptation occurs
        - if PETSc uses `KSP.NormType.PRECONDITIONED`, then it progresses 
          normally.
        - with `legacy`, norm is `KSP.NormType.DEFAULT`, which gets 
          translated to `KSP.NormType.PRECONDITIONED`, since the solve 
          gets an ILU preconditioner by default.
        - with `legacy`, SciPy normalizes by the initial residual,
          as opposed ot using the RHS norm without.
        - ditto for Trilinos
        - ditto for PySparse
        - Initial bnorm=3.4e10, rnorm = 3.3e-2
        - Setting tolerance = default / 1e12 = 1e-22 should work, and sort 
          of does, but PETSc throws lots of `KSP_DIVERGED_BREAKDOWN`
          and Trilinos throws lots of `AZ_maxits`. SciPy and PySparse both 
          use LU, so don't (can't) throw these errors.
        - Trilinos throws `AZ_maxits`, regardless, although not as many.
        - Ratio of time with tolerance=1e-22 vs legacy for 200 totalSteps:
            - PETSc GMRES +1.6%
            - PySparse LU +2.5%
            - SciPy LU +0.3%
            - Trilinos GMRES +7%
        - *** `criterion="initial"` resolves the issue. ***
    - [x] examples/reactiveWetting/liquidVapor2D.py
        - see examples/reactiveWetting/liquidVapor1D.py
    - These were raising `DivergenceWarning` before:
        - [x] fipy.meshes.cylindricalUniformGrid1D.CylindricalUniformGrid1D
        - [x] fipy.meshes.cylindricalNonUniformGrid1D.CylindricalNonUniformGrid1D
        - [x] examples.reactiveWetting.liquidVapor2D
- [ ] Other benchmarks
    - [ ] Coupled diffusion/phase transformation
    - [ ] Cahn-Hilliard
    - [ ] Coupled Cahn-Hilliard
    - [ ] Flow
    - [ ] Reactive wetting
    - [ ] Circular mesh
    - [ ] 3D problem
- [x] Solver tests of normalizations
- [ ] PyAMGX
- [x] Ensure Trilinos residual is a number and not an array
- [x] Fix off-by-one in PySparse iteration count
- [x] Check parallel
- [x] Log build and solve times at higher resolution
      (instead of relying on log times, which are limited to millisecond 
      resolution and which include system time)
- [x] Default tolerance -> 1e-5?
    - [x] fipy.terms.term.Term._test
    - [x] fipy.terms.transientTerm.TransientTerm
    - [x] examples.diffusion.steadyState.mesh1D.tri2Dinput
    - [x] examples.diffusion.steadyState.mesh20x20.tri2Dinput
    - [x] examples.diffusion.steadyState.mesh50x50.input
    - [x] examples.diffusion.steadyState.mesh50x50.tri2Dinput
    - [x] examples.diffusion.steadyState.otherMeshes.grid3Dinput
    - [x] examples.diffusion.mesh1D
    - [x] examples.diffusion.variable
    - [x] examples.phase.impingement.mesh40x1
    - [x] examples.convection.exponential1D.tri2D
    - [x] examples.convection.exponential2D.tri2D
- [x] Try turning divergence_tolerance back on
    - [x] examples.diffusion.steadyState.mesh1D.inputPeriodic
        - Lnorm: 4.00000001, bnorm: 4.0822787753900393e-08, rnorm: 1.414213562373095
        - Using `criterion="initial"` raises
          ```
          [0] Relative tolerance 346.427 must be non-negative and less than 1.0
          ```
          because (rnorm / bnorm) * tolerance = 346.427
        - This didn't happen before because 
          `ksp.setInitialGuessNonzero(True)` wasn't set, so `DTOL` was 
          never checked.
        - Use `.setConvergenceTest()`?
            - How to integrate with rest of solvers?
        - Raise `divergence_tolerance`?
            - How to integrate with rest of solvers?
    - [x] examples.diffusion.mesh1D
    - [x] examples.phase.binary
    - [x] examples.phase.quaternary
    - [x] examples.elphf.diffusion.mesh1D
    - [x] examples.elphf.phaseDiffusion
    - [x] examples.reactiveWetting.liquidVapor1D
- [ ] Convert uses of `.flat` and `.flatten()` to `.ravel()` or 
      `.reshape((-1,))`
- [ ] Consider "RHS" -> "rhs"
- [ ] Redo parallel scaling plots
