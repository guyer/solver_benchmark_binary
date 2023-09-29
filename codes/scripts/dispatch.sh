#!/bin/bash

# E.g.,
# bash codes/scripts/dispatch.sh --scheduler bash --env fipy27 --solversuite pysparse --log codes/loggers/config_template.json --output results/problem/platform --preconditioners "jacobi none" --codes/notebooks/diffusion.ipynb --store_by_solver
# or
# bash codes/scripts/dispatch.sh --scheduler sbatch --queue fast --skedargs "--time=12:00:00"  --env fipy_solver_benchmarking_11 --solversuite trilinos --log codes/loggers/config_template.json --output results/binary/linux -- codes/notebooks/binary_phase_field.ipynb --store_by_solver

USAGE="usage: $0 [-h] [OPTIONS] [--] NOTEBOOK [ARGS]

Converts NOTEBOOK to SCRIPT,
iterates over solvers, preconditionrs, and mesh sizes by calling setup.sh,
which activates the appropriate conda environment and calls python on SCRIPT.

positional arguments:
  NOTEBOOK    Jupyter notebook to convert to python and launch.
              Can be absolute path or relative to directory where
              dispatch.sh is executed.

optional arguments:
  -h, --help  show this help message and exit
  --scheduler SCHEDULER  Tool used to launch job.
              Can be 'sbatch' for Slurm, 'qsub' for Sun Grid Engine, or 'bash'
              (default: invoke using bash)
  --skedargs SKEDARGS  Space-separated string of argments to pass to the SCHEDULER.
  --queue QUEUE  Name of queue or partition on the SCHEDULER.
  --env ENV   Conda environment to activate before invoking SCRIPT
              (default: fipy)
  --np NP     Number of processes to invoke SCRIPT with (default: 1)
  --mprof     Whether to run mprof profiler (default: False)
  --output OUTPUT   Directory to store results in
  --log CONFIG  Path to log configuration file template.
  --solversuite SUITE   Solver package to use (default: petsc)
  --solvers SOLVERS  Names of solvers (separated by spaces)
              (default: 'pcg cgs gmres lu')
  --preconditioners PRECONDITIONERS  Names of preconditioners (separated by spaces)
              (default: 'jacobi ilu ssor icc none')
  --powermin POWERMIN   Power of ten for minimum size, minsize = 10**POWERMIN (default: 1)
  --powermax POWERMAX   Power of ten for maximum size, maxsize = 10**POWERMAX (default: 6)
  --powerstep POWERSTEP Increment in power of ten for size (default: 1)"

SCHEDULER="bash"
SKEDARGS=""
QUEUE=""
SLURMTIME="1:00:00"
ENV=fipy
MPI=""
NP=1
LOGCONFIG=""
LOGNAME=""
OUTPUT="."
PYTHON=python
SOLVERSUITE=petsc
SOLVERS="pcg cgs gmres lu"
PRECONDITIONERS="jacobi ilu ssor icc none"
POWERMIN=1
POWERMAX=6
POWERSTEP=1

function make_absolute () {
    case $1 in
        /*)
            # absolute path
            echo $1
            ;;
        *)
            # relative path, make absolute
            echo "$(pwd)/$1"
    esac
}

while [[ $# > 0 ]] && [[ $1 == -* ]]
do
    case "$1" in
        --scheduler)
            SCHEDULER="$2"
            shift # option has parameter
            ;;
        --skedargs)
            SKEDARGS="$2"
            shift # option has parameter
            ;;
        --queue)
            QUEUE="$2"
            shift # option has parameter
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
            shift # option has parameter
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
        --solvers)
            SOLVERS="$2"
            shift # option has parameter
            ;;
        --preconditioners)
            PRECONDITIONERS="$2"
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

if [[ $NP > 1 ]]; then
    MPI="mpirun -np ${NP}"
fi

OUTPUT="$(make_absolute ${OUTPUT})"
readonly OUTPUT

mkdir -p "${OUTPUT}"

NOTEBOOK=$1
shift

nbpath=${NOTEBOOK%/*}
nbbase=${NOTEBOOK##*/}
nbpref=${nbbase%.*}
nbfext=${nbbase##*.}

# https://stackoverflow.com/a/56155771/2019542
eval "$(conda shell.bash hook)"
conda activate $ENV

set -x
jupyter nbconvert ${NOTEBOOK} --to python --output-dir=${OUTPUT}
set +x

conda deactivate

script="${nbpref}.py"

declare -a setup

setup=($(make_absolute "${BASH_SOURCE%/*}/setup.sh"))
if [[ -n ${LOGCONFIG} ]]; then
    setup+=(--log "$(make_absolute ${LOGCONFIG})")
fi
setup+=(--env "${ENV}")

for (( power=${POWERMIN}; power<=${POWERMAX}; power+=${POWERSTEP} ))
do
    size=$((10**power))
    for solver in ${SOLVERS}
    do
        for preconditioner in ${PRECONDITIONERS}
        do
            job=(OMP_NUM_THREADS=1 FIPY_SOLVERS=${SOLVERSUITE})
            job+=(${MPI} ${PYTHON} ${script})
            job+=(--numberOfElements=${size})
            job+=(--solver="${solver}")
            job+=(--preconditioner="${preconditioner}")
            job+=($@)

            case "${SCHEDULER}" in
                qsub)
                    options=(-cwd -pe nodal "${NP}" -q "${QUEUE}")
                    options+=(-o "${OUTPUT}" -e "${OUTPUT}")
                    options+=(change working directory to ${OUTPUT})
                    options+=(${SKEDARGS})

                    set -x
                    qsub ${options[@]} -- ${setup[@]} -- ${job[@]}
                    set +x
                    ;;
                sbatch)
                    options=(--partition="${QUEUE}")
                    options+=(--job-name="${script##*/}-${SOLVERSUITE}-${size}-${preconditioner}-${solver}")
                    options+=(--ntasks="${NP}" --ntasks-per-core=2)
                    options+=(--chdir="${OUTPUT}")
                    options+=(${SKEDARGS})

                    set -x
                    sbatch ${options[@]} -- ${setup[@]} -- ${job[@]}
                    set +x
                    ;;
                bash)
                    # default
                    set -x
                    (pushd "${OUTPUT}" && ${setup[@]} -- ${job[@]}; popd)
                    set +x
                    ;;
                *)
                    echo Unknown scheduler: "${SCHEDULER}">&2
                    exit 20
                    ;;
            esac
        done
    done
done
