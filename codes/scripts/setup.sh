#!/bin/bash

USAGE="usage: $0 [-h] [OPTIONS] [--] SCRIPT [ARGS]

activates the appropriate conda environment and calls python on SCRIPT
with ARGS

optional arguments:
  -h, --help  show this help message and exit
  --env ENV   Conda environment to activate before invoking SCRIPT
              (default: fipy)
  --log CONFIG  Path to log configuration file template.
"

ENV=fipy
LOGFILE="solver.log"

while [[ $# > 0 ]] && [[ $1 == -* ]]
do
    case "$1" in
        --env)
            ENV="$2"
            shift # option has parameter
            ;;
        --log)
            LOGCONFIG="$2"
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

if [[ -n $LOGCONFIG ]]; then
    configbase=${LOGCONFIG##*/}
    configpref=${configbase%.*}
    configfext=${configbase##*.}
    
    JOBLOGCONFIG="${configpref}.${configfext}"

    if [[ -n $SLURM_JOB_ID ]]; then
        logbase=${LOGFILE##*/}
        logpref=${logbase%.*}
        logfext=${logbase##*.}

        LOGFILE="${logpref}.${SLURM_JOB_ID}.${logfext}"
        JOBLOGCONFIG="${configpref}.${SLURM_JOB_ID}.${configfext}"
    fi

    sed -e "s:%LOGFILE%:${LOGFILE}:g" "${LOGCONFIG}" > "${JOBLOGCONFIG}"

    LOGCONFIGENV="FIPY_LOG_CONFIG=${JOBLOGCONFIG}"
fi

# https://stackoverflow.com/a/56155771/2019542
eval "$(conda shell.bash hook)"
conda activate $ENV

set -x
env ${LOGCONFIGENV} "$@" --output=.
set +x
