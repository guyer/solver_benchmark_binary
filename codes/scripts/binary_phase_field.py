#!/usr/bin/env python
# coding: utf-8

# # Binary Phase Field Benchmark
# 
# FiPy implementation of phase transormation in 2D with solutal driving force

# **Do not edit `binary_phase_field.py`**. Generate the batch-runnable file from the notebook with
# ```bash
# jupyter nbconvert binary_phase_field.ipynb --to python --output-dir=../scripts/
# ```

# ## Import Python modules

# In[1]:


import argparse
import json
import os
import re
import sys

import fipy as fp
from fipy.tools import numerix as nmx
from fipy.tools import parallelComm


# Jupyter notebook handles some things differently than from the commandline

# In[2]:


try:
    from IPython import get_ipython
    isnotebook = get_ipython() is not None
except:
    isnotebook = False


# ## Initialize
# ### Load parameters

# In[3]:


parser = argparse.ArgumentParser()
parser.add_argument("--output", help="directory to store results in",
                    default=None)
parser.add_argument("--store_by_solver",
                    help="store results in nested subdirectories based on solver,"
                    "preconditioner, and system size",
                    action='store_true')
parser.add_argument("--restart", help="solution to initialize from",
                    default=None)
parser.add_argument("--checkpoint_interval", help="frequency to save results",
                    type=float, default=6.)
parser.add_argument("--totaltime", help="duration of full simulation",
                    type=float, default=600.)
parser.add_argument("--numnuclei", help="number of nuclei",
                    type=int, default=100)
parser.add_argument("--factor", help="fraction of critical nucleus size for new nuclei",
                    type=float, default=1.1)
parser.add_argument("--nucleation_scale", help="size of domain for nuclei",
                    type=float, default=1000)
parser.add_argument("--numberOfElements", help="number of total cells in a Grid2D",
                    type=int, default=1000000)
parser.add_argument("--solver", help="solver class to use",
                    choices=("pcg", "cgs", "gmres", "lu"), default="pcg")
parser.add_argument("--preconditioner", help="preconditioner class to use",
                    choices=("jacobi", "ilu", "ssor", "icc", "none"), default="none")
parser.add_argument("--sweeps", help="number of nonlinear sweeps to take",
                    type=int, default=5)
parser.add_argument("--iterations", help="maximum number of linear iterations to take for each sweep",
                    type=int, default=1000)
parser.add_argument("--tolerance", help="linear solver tolerance",
                    type=float, default=1e-10)
parser.add_argument("--store_matrix",
                    help="store the matrix and RHS vector along with other output",
                    action='store_true')


# ### Set any parameters for interactive notebook

# In[81]:


if isnotebook:
    # argv = ["--numberOfElements=10000", "--totaltime=1.2", "--checkpoint_interval=0.12",
    #         "--nucleation_scale=100", "--output=nucleation6"]
    argv = ["--numberOfElements=10000", "--nucleation_scale=1000", "--output=solidification1",
           "--restart=../../results/nucleation/nucleation18/t=300.0.npz",
           "--store_matrix"]
else:
    argv = None


# In[82]:


args, unknowns = parser.parse_known_args(args=argv)


# ### Initialize mesh and solution variables
# 
# Either restart from some `path/to/restart/t={time}.npz`, where the time is assigned to `elapsed`
# 
# or
# 
# Create a mesh based on parameters. Set
# >  the computational domain is ... 1000×1000 

# In[94]:


nx = ny = int(nmx.sqrt(args.numberOfElements))
mesh = fp.Grid2D(nx=nx, ny=ny)
phi = fp.CellVariable(mesh=mesh, name="$\phi$", value=0., hasOld=True)
elapsed = 0.


# In[95]:


if args.restart is not None:
    data = nmx.load(args.restart)
    lower, upper = int((1000 - nx) / 2), int((1000 + nx) / 2)
    phi.setValue(data["phi"][lower:upper, lower:upper].flat)

    # scanf("%g") simulator
    # https://docs.python.org/3/library/re.html#simulating-scanf
    scanf_g = "[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
    pattern = ".*t=({g})\.npz".format(g=scanf_g)
    elapsed = float(re.match(pattern, args.restart).group(1))


# In[96]:


x, y = mesh.cellCenters[0], mesh.cellCenters[1]
X, Y = mesh.faceCenters[0], mesh.faceCenters[1]


# In[97]:


C = fp.CellVariable(mesh=mesh, name="$C$", value=0.3 + 0.4 * phi, hasOld=True)


# In[98]:


if isnotebook:
    viewer = fp.Viewer(vars=C, datamin=0., datamax=1.)
    viewer.plot()


# ## Create solver

# In[99]:


precon = None

if args.preconditioner == "jacobi":
    precon = fp.JacobiPreconditioner()
elif args.preconditioner == "ilu":
    precon = fp.ILUPreconditioner()
elif args.preconditioner == "ssor":
    precon = fp.SSORPreconditioner()
elif args.preconditioner == "icc":
    precon = fp.ICPreconditioner()
elif args.preconditioner == "none":
    precon = None

if args.solver == "cgs":
    solver_class = fp.LinearCGSSolver
elif args.solver == "gmres":
    solver_class = fp.LinearGMRESSolver
elif args.solver == "lu":
    if args.preconditioner != "none":
        # preconditioned lu doesn't make any sense
        exit()

    solver_class = fp.LinearLUSolver
elif args.solver == "pcg":
    solver_class = fp.LinearPCGSolver

solver = solver_class(tolerance=args.tolerance, criterion="initial",
                      iterations=args.iterations, precon=precon)


# ## Define governing equation

# ## Binary Alloy Phase Field
# 
# binary diffusion with frozen phase; diffusion with convection due to phase boundary

# \begin{align}
# \frac{1}{M_\phi}\frac{\partial \phi}{\partial t} &= \nabla\cdot\left(\kappa_\phi\nabla \phi\right)
# - \left[
#     \left(1-C\right)\frac{L_A\left(T - T_m^A\right)}{T_m^A}
#     + C\frac{L_B\left(T - T_m^B\right)}{T_m^B}
# \right] p'(\phi)
# - \left[
#     \left(1-C\right)\frac{W_A}{2}
#     + C\frac{W_B}{2}
# \right] g'(\phi)
# \\
# \frac{\partial C}{\partial t} &= \nabla\cdot\left(D_C\nabla C\right)
# \\ &\quad {} + \nabla\cdot\left(
#     D_C\frac{C\left(1-C\right)V_m}{RT}
#     \left\{
#         \left[
#             \frac{L_B\left(T - T_m^B\right)}{T_m^B}
#             - \frac{L_A\left(T - T_m^A\right)}{T_m^A}
#         \right] p'(\phi)
#         - \frac{W_B - W_A}{2} g'(\phi)
#     \right\}
#     \nabla \phi\right)
# \end{align}

# non-dimensionalize

# \begin{align}
# \frac{1}{M_\phi}\frac{\partial \phi}{\partial t \tau} &= \frac{1}{\lambda}\nabla\cdot\left(\kappa_\phi\frac{1}{\lambda}\nabla \phi\right)
# - \left[
#     \left(1-C\right)L_A\left(\frac{T}{T_m^A} - 1\right)
#     + C L_B\left(\frac{T}{T_m^B} - 1\right)
# \right] p'(\phi)
# - \left[
#     \left(1-C\right)\frac{W_A}{2}
#     + C\frac{W_B}{2}
# \right] g'(\phi)
# \\
# \frac{\partial C}{\partial t \tau} &= \frac{1}{\lambda}\nabla\cdot\left(D_C\frac{1}{\lambda}\nabla C\right)
# \\ &\quad {} + \frac{1}{\lambda} \nabla\cdot\left(
#     D_C\frac{C\left(1-C\right)V_m}{RT}
#     \left\{
#         \left[
#             L_B\left(\frac{T}{T_m^B} - 1\right)
#             - L_A\left(\frac{T}{T_m^A} - 1\right)
#         \right] p'(\phi)
#         - \frac{W_B - W_A}{2} g'(\phi)
#     \right\}
#     \frac{1}{\lambda} \nabla \phi\right)
# \end{align}

# \begin{align}
# \frac{\partial \phi}{\partial t} &= \nabla\cdot\left(\frac{M_\phi\tau\kappa_\phi}{\lambda^2}\nabla \phi\right)
# - \left[
#     \left(1-C\right)M_\phi\tau L_A\left(\frac{T}{T_m^A} - 1\right)
#     + C M_\phi\tau L_B\left(\frac{T}{T_m^B} - 1\right)
# \right] p'(\phi)
# - \left[
#     \left(1-C\right)\frac{M_\phi\tau W_A}{2}
#     + C\frac{M_\phi\tau W_B}{2}
# \right] g'(\phi)
# \\
# \frac{\partial C}{\partial t} &= \nabla\cdot\left(\frac{D_C\tau}{\lambda^2}\nabla C\right)
# \\ &\quad {} + \nabla\cdot\left(
#     \frac{D_C\tau}{\lambda^2}C\left(1-C\right)
#     \left\{
#         \left[
#             \frac{L_B V_m}{RT}\left(\frac{T}{T_m^B} - 1\right)
#             - \frac{L_A V_m}{RT}\left(\frac{T}{T_m^A} - 1\right)
#         \right] p'(\phi)
#         - \frac{\frac{W_B V_m}{RT} - \frac{W_A V_m}{RT}}{2} g'(\phi)
#     \right\}
#     \nabla \phi\right)
# \end{align}

# \begin{align}
# M_\phi &= \frac{T_m^A \beta_A}{6 L_A \delta_A}
# \\
# &\sim \mathrm{\frac{K \frac{cm}{K s}}{\frac{J}{cm^3} cm}}
# \\
# &\sim \mathrm{\frac{cm^3}{J\,s}}
# \end{align}

# \begin{align}
# \tau M_\phi &\sim \mathrm{s \frac{cm^3}{J\,s}}
# \\
# &\sim \mathrm{\frac{cm^3}{J}}
# \end{align}

# \begin{align}
# \tau M_\phi \frac{\kappa_\phi}{\lambda^2} &\sim \mathrm{\frac{cm^3}{J}\frac{J}{cm}\frac{1}{cm^2}}
# \\
# &\sim 1
# \end{align}

# \begin{align}
# \tau M_\phi L_? \sim \tau M_\phi W_? &\sim \mathrm{\frac{cm^3}{J}\frac{J}{cm^3}}
# \\
# &\sim 1
# \end{align}

# \begin{align}
# \frac{D_C \tau}{\lambda^2} &\sim \mathrm{\frac{cm^2}{s}\frac{s}{cm^2}}
# \\
# &\sim 1
# \end{align}

# \begin{align}
# \frac{L_? V_m}{R T} \sim \frac{W_? V_m}{R T} &\sim \mathrm{\frac{J}{cm^3}\frac{cm^3}{mol}\frac{mol\,K}{J\,K}}
# \\
# &\sim 1
# \end{align}

# \begin{align}
# \frac{\partial \phi}{\partial t} &= \nabla\cdot\left(\tilde{\kappa}_\phi\nabla \phi\right)
# - \left[
#     \left(1-C\right)\tilde{L}_A\left(\frac{T}{T_m^A} - 1\right)
#     + C \tilde{L}_B\left(\frac{T}{T_m^B} - 1\right)
# \right] p'(\phi)
# - \left[
#     \left(1-C\right)\frac{\tilde{W}_A}{2}
#     + C\frac{\tilde{W}_B}{2}
# \right] g'(\phi)
# \\
# \frac{\partial C}{\partial t} &= \nabla\cdot\left(\tilde{D}_C\nabla C\right)
# \\ &\quad {} + \nabla\cdot\left(
#     \tilde{D}_C C\left(1-C\right)
#     \left\{
#         \left[
#             \hat{L}_B\left(\frac{T}{T_m^B} - 1\right)
#             - \hat{L}_A\left(\frac{T}{T_m^A} - 1\right)
#         \right] p'(\phi)
#         - \frac{\hat{W}_B - \hat{W}_A}{2} g'(\phi)
#     \right\}
#     \nabla \phi\right)
# \end{align}

# Let $V_m / (R T) = 100 \tau M_\phi$, such that $\hat{L}_? = 100 \tilde{L}_?$ and $\hat{W}_? = 100 \tilde{W}_?$:
# \begin{align}
# \frac{\partial \phi}{\partial t} &= \nabla\cdot\left(\tilde{\kappa}_\phi\nabla \phi\right)
# - \left[
#     \left(1-C\right)\tilde{L}_A\left(\frac{T}{T_m^A} - 1\right)
#     + C \tilde{L}_B\left(\frac{T}{T_m^B} - 1\right)
# \right] p'(\phi)
# - \left[
#     \left(1-C\right)\frac{\tilde{W}_A}{2}
#     + C\frac{\tilde{W}_B}{2}
# \right] g'(\phi)
# \\
# \frac{\partial C}{\partial t} &= \nabla\cdot\left(\tilde{D}_C\nabla C\right)
# \\ &\quad {} + \nabla\cdot\left(
#     100
#     \tilde{D}_C C\left(1-C\right)
#     \left\{
#         \left[
#             \tilde{L}_B\left(\frac{T}{T_m^B} - 1\right)
#             - \tilde{L}_A\left(\frac{T}{T_m^A} - 1\right)
#         \right] p'(\phi)
#         - \frac{\tilde{W}_B - \tilde{W}_A}{2} g'(\phi)
#     \right\}
#     \nabla \phi\right)
# \end{align}

# Comparing the phase equation with the [nucleation benchmark](nucleation.ipynb) equation
# $$
# \begin{aligned}
# \frac{\partial\phi}{\partial t} &= \nabla^2 \phi - g'(\phi) + \Delta f p'(\phi)
# \\
# \Delta f &= \frac{1}{6\sqrt{2}}
# \end{aligned}
# $$

# \begin{align}
# \tilde{\kappa}_\phi &=1
# \\
# \tilde{W}_A = \tilde{W}_B &= 1
# \\
# \tilde{L}_A \left(\frac{T}{T_m^A} - 1\right) = \tilde{L}_B \left(\frac{T}{T_m^B} - 1\right)
# &= -\frac{1}{6\sqrt{2}}
# \end{align}
# which would result in no convection in the diffusion equation. So, let
# \begin{align}
# \tilde{L}_A \left(\frac{T}{T_m^A} - 1\right) &= -\frac{1.1}{6\sqrt{2}}
# \\
# \tilde{L}_B \left(\frac{T}{T_m^B} - 1\right) &= -\frac{0.9}{6\sqrt{2}}
# \\
# \tilde{D}_C &= 1
# \end{align}

# \begin{align}
# \frac{\partial \phi}{\partial t} &= \nabla^2 \phi
# + \left[
#     1.1\left(1-C\right)
#     + 0.9 C
# \right] \frac{1}{6\sqrt{2}} p'(\phi)
# - g'(\phi)
# \\
# \frac{\partial C}{\partial t} &= \nabla^2 C
# + \nabla\cdot\left(
#     C\left(1-C\right) \frac{20}{6\sqrt{2}} p'(\phi)
#     \nabla \phi\right)
# \\
# &= \nabla^2 C
# + \nabla\cdot\left(
#     C\left(1-C\right) \frac{20}{6\sqrt{2}} \nabla p(\phi)\right)
# \end{align}

# Freeze the phase field and solve the diffusion equation

# In[100]:


def p(phi):
    return phi**3 * (6 * phi**2 - 15 * phi + 10)

# DeltaDelta_f = 20. / (6 * nmx.sqrt(2.))
DeltaDelta_f = 20. / (6 * nmx.sqrt(2.))
phaseTransformationVelocity = (1 - C).harmonicFaceValue * DeltaDelta_f * p(phi).faceGrad

eq = (fp.TransientTerm(var=C)
      == fp.DiffusionTerm(var=C)
      + fp.ConvectionTerm(coeff=phaseTransformationVelocity, var=C))


# ## Setup output

# ### Setup ouput storage

# In[105]:


if (args.output is not None) and (parallelComm.procID == 0):
    if args.store_by_solver:
        suite = solver.__module__.split('.')[2]
        if args.preconditioner is None:
            preconditioner_name = "none"
        else:
            preconditioner_name = precon.__class__.__name__
        path = os.path.join(args.output, suite,
                            solver.__class__.__name__,
                            preconditioner_name,
                            str(nx * ny))
    else:
        path = args.output

    os.makedirs(path)

if parallelComm.procID == 0:
    print("storing results in {0}".format(path))


# ### Define output routines

# In[104]:


def saveC(elapsed):
    C_value = C.globalValue.reshape((nx, ny))
    if parallelComm.procID == 0:
        fname = os.path.join(path, "t={}.npz".format(elapsed))
        nmx.savez(fname, C=C_value)

def saveMatrix(elapsed):
    mtxname = os.path.join(path, "t={}.mtx".format(elapsed))
    eq.matrix.exportMmf(mtxname)
    
    rhs_value = eq.RHSvector
    if parallelComm.procID == 0:
        rhsname = os.path.join(path, "t={}.rhs.npz".format(elapsed))
        nmx.savez(rhsname, rhs=rhs_value)

def checkpoint_data(elapsed, store_matrix=False):
    saveC(elapsed)
    if store_matrix:
        saveMatrix(elapsed)


# ### Figure out when to save

# In[101]:


checkpoints = (fp.numerix.arange(int(elapsed / args.checkpoint_interval),
                                 int(args.totaltime / args.checkpoint_interval))
               + 1) * args.checkpoint_interval

checkpoints.sort()


# ## Solve and output

# In[102]:


times = checkpoints
times = times[(times > elapsed) & (times <= args.totaltime)]


# In[103]:


from steppyngstounes import CheckpointStepper, FixedStepper

eq.cacheMatrix()
eq.cacheRHSvector()

C.updateOld()
for checkpoint in CheckpointStepper(start=elapsed,
                                    stops=times,
                                    stop=args.totaltime):

    for step in FixedStepper(start=checkpoint.begin,
                             stop=checkpoint.end,
                             size=100.):

        state = dict(state="START", numberOfElements=mesh.numberOfCells, sweeps=args.sweeps)
        if precon is None:
            state["preconditioner"] = None
        else:
            state["preconditioner"] = precon.__class__.__name__

        solver._log.debug(json.dumps(state))

        for sweep in range(args.sweeps):
            res = eq.sweep(var=C, dt=step.size, solver=solver)

        state["state"] = "END"
        solver._log.debug(json.dumps(state))

        C.updateOld()
        # stats.append(current_stats(step.end))

        _ = step.succeeded(error=res / 1e-3)

    if checkpoint.end in checkpoints:
        # don't save nucleation events?
        checkpoint_data(checkpoint.end, store_matrix=args.store_matrix)

    if isnotebook:
        viewer.plot()
        # labelViewer.plot()

    _ = checkpoint.succeeded()

