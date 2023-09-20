#!/bin/bash

# E.g.,
# bash codes/scripts/dispatch.sh --env fipy27 --solversuite pysparse --log codes/loggers/config_template.json solver.log --output=results/problem/platform codes/notebooks/diffusion.ipynb --preconditioners=ilu --store_by_solver

USAGE="usage: $0 [-h] [OPTIONS] [--] NOTEBOOK [ARGS]

Iterates over solvers and mesh sizes by calling setup.sh, which activates
the appropriate conda environment and calls python on SCRIPT

positional arguments:
  SCRIPT      Python script to launch (expected to be in same
              directory as $0)

optional arguments:
  -h, --help  show this help message and exit
  --qsub      Invoke SCRIPT using 'qsub -cwd' for Sun grid engine
              (default: invoke using bash)
  --sbatch PARTITION SLURMTIME
              Invoke SCRIPT using 'sbatch --partition=${PARTITION} --time=${SLURMTIME}'
              for slurm (default: invoke using bash)
  --env ENV   Conda environment to activate before invoking SCRIPT
              (default: fipy)
  --np NP     Number of processes to invoke SCRIPT with (default: 1)
  --mprof     Whether to run mprof profiler (default: False)
  --output OUTPUT   Directory to store results in
  --log CONFIG LOGNAME  Path to log configuration file template and
                        name for log file.
  --solversuite SUITE   Solver package to use (default: petsc)
  --powermin POWERMIN   Power of ten for minimum size, minsize = 10**POWERMIN (default: 1)
  --powermin POWERMAX   Power of ten for maximum size, maxsize = 10**POWERMAX (default: 6)
  --powerstep POWERSTEP Increment in power of ten for size (default: 1)
  --preconditioners PRECONDITIONERS  Names of preconditioners (separated by spaces) (default: none)"

QSUB=0
SBATCH=0
PARTITION=""
SLURMTIME="1:00:00"
ENV=fipy
NP=1
LOGCONFIG=""
LOGNAME=""
OUTPUT="."
PYTHON=python
SOLVERSUITE=petsc
PRECONDITIONERS=none
POWERMIN=1
POWERMAX=6
POWERSTEP=1

while [[ $# > 0 ]] && [[ $1 == -* ]]
do
    case "$1" in
        --qsub)
            QSUB=1
            ;;
        --sbatch)
            SBATCH=1
            PARTITION="$2"
            SLURMTIME="$3"
            shift # option has two parameters
            shift
            ;;
        --env)
            ENV="$2"
            shift # option has parameter
            ;;
        --np)
            NP="$2"
            shift # option has parameter
            ;;
        --log)
            LOGCONFIG="$2"
            LOGNAME="$3"
            shift # option has two parameters
            shift
            ;;
        --mprof)
            PYTHON="mprof run"
            ;;
        --output)
            OUTPUT="$2"
            shift # option has parameter
            ;;
        --solversuite)
            SOLVERSUITE="$2"
            shift # option has parameter
            ;;
        --powermin)
            POWERMIN="$2"
            shift # option has parameter
            ;;
        --powermax)
            POWERMAX="$2"
            shift # option has parameter
            ;;
        --powerstep)
            POWERMAX="$2"
            shift # option has parameter
            ;;
        --preconditioners)
            PRECONDITIONERS="$2"
            shift # option has parameter
            ;;
        -h|--help)
            echo "$USAGE"
            exit 0
            ;;
        --)
            # end of options
            shift
            break
            ;;
        -*)
            # unknown option
            echo Unknown option: $1>&2
            exit 10
            ;;
    esac
    shift # option(s) fully processed, proceed to next input argument
done

if [[ "$#" < 1 ]]; then
    echo "$USAGE"
    exit 1
fi

SCRIPT=$1
shift

if [[ $NP > 1 ]]; then
    MPI="mpirun -np ${NP}"
else
    MPI=""
fi

for (( POWER=${POWERMIN}; POWER<=${POWERMAX}; POWER+=${POWERSTEP} ))
do
    size=$((10**${POWER}))
    for solver in pcg cgs gmres lu
    do
        for preconditioner in $PRECONDITIONERS
        do
            INVOCATION="OMP_NUM_THREADS=1 FIPY_SOLVERS=${SOLVERSUITE} ${LOG_CONFIG} \
              ${MPI} ${PYTHON} ${BASH_SOURCE%/*}/${SCRIPT} \
              --numberOfElements=${size} --solver=${solver} --preconditioner=${preconditioner} $@"

            JOBNAME="${SCRIPT}-${SOLVERSUITE}-${size}-${preconditioner}-${solver}"

            if [[ $QSUB == 1 ]]; then
                qsub -cwd -pe nodal ${NP} -q "wide64" -o "${dir}" -e "${dir}" \
                  "${BASH_SOURCE%/*}/setup.sh" \
                  --env "${ENV}" --output "${OUTPUT}" \
                  -- ${INVOCATION}
            elif [[ $SBATCH == 1 ]]; then
                sbatch --partition=${PARTITION} --job-name=${JOBNAME} --ntasks=${NP} \
                  --ntasks-per-core=2 --time=${SLURMTIME} \
                  "${BASH_SOURCE%/*}/setup.sh" \
                  --log ${LOGCONFIG} ${LOGNAME} --env "${ENV}" --output "${OUTPUT}" \
                  -- ${INVOCATION}
            else
                echo bash "${BASH_SOURCE%/*}/setup.sh" \
                  --log ${LOGCONFIG} ${LOGNAME} --env "${ENV}" --output "${OUTPUT}" \
                  -- ${INVOCATION}
                bash "${BASH_SOURCE%/*}/setup.sh" \
                  --log ${LOGCONFIG} ${LOGNAME} --env "${ENV}" --output "${OUTPUT}" \
                  -- ${INVOCATION}
            fi
        done
    done
done
