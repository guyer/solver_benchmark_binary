2019-12-19T11:58:22-05:00
=========================

If Krylov has default preconditioning (ilu)
-------------------------------------------

fipy.terms.binaryTerm._BinaryTerm._test
examples.levelSet.electroChem.howToWriteAScript
  ValueError: the array phi contains no zero contour (no zero level set)
examples.levelSet.electroChem.simpleTrenchSystem
  ValueError: the array phi contains no zero contour (no zero level set)
examples.levelSet.electroChem.gold
  ValueError: the array phi contains no zero contour (no zero level set)
examples.levelSet.electroChem.leveler
  ValueError: the array phi contains no zero contour (no zero level set)
  
If LU is normalized by maxdiag
------------------------------

fipy.terms.binaryTerm._BinaryTerm._test
examples.levelSet.electroChem.howToWriteAScript
  ValueError: the array phi contains no zero contour (no zero level set)
examples.levelSet.electroChem.simpleTrenchSystem
  ValueError: the array phi contains no zero contour (no zero level set)
examples.levelSet.electroChem.gold
  ValueError: the array phi contains no zero contour (no zero level set)
examples.levelSet.electroChem.leveler
  ValueError: the array phi contains no zero contour (no zero level set)
examples.cahnHilliard.mesh2D

If Krylov has jacobi preconditioning and LU is normalized by maxdiag
--------------------------------------------------------------------

fipy.terms.term.Term._test
fipy.terms.binaryTerm._BinaryTerm._test
examples.diffusion.steadyState.mesh50x50.tri2Dinput
examples.diffusion.nthOrder.input4thOrder1D
examples.diffusion.mesh1D
examples.diffusion.electrostatics
examples.phase.binary
examples.phase.quaternary
examples.convection.powerLaw1D.mesh1D
examples.convection.peclet
examples.convection.source
examples.elphf.diffusion.mesh2D
examples.elphf.phaseDiffusion
examples.elphf.poisson
examples.levelSet.electroChem.gold
examples.cahnHilliard.tanh1D
examples.cahnHilliard.mesh2D
examples.flow.stokesCavity

If Krylov has jacobi preconditioning
------------------------------------

fipy.terms.binaryTerm._BinaryTerm._test
examples.diffusion.steadyState.mesh50x50.tri2Dinput
examples.diffusion.nthOrder.input4thOrder1D
examples.diffusion.mesh1D
examples.diffusion.electrostatics
examples.phase.binary
examples.phase.quaternary
examples.convection.powerLaw1D.mesh1D
examples.convection.peclet
examples.convection.source
examples.elphf.diffusion.mesh2D
examples.elphf.phaseDiffusion
examples.elphf.poisson
examples.levelSet.electroChem.gold
examples.cahnHilliard.tanh1D
examples.flow.stokesCavity

If all are normalized by maxdiag
--------------------------------

examples.phase.impingement.mesh20x20
examples.levelSet.electroChem.gold
examples.levelSet.electroChem.adsorbingSurfactantEquation.AdsorbingSurfactantEquation
examples.reactiveWetting.liquidVapor1D

How could fipy.terms.binaryTerm._BinaryTerm._test be impacted by this 
change? Is some solver state changed?
