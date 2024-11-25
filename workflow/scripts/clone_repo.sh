#!/usr/bin/env bash

exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

git clone --filter=blob:none "${snakemake_params[fipy_repo]}" \
  "${snakemake_output[repo]}"

pushd "${snakemake_output[repo]}"
git checkout "${snakemake_wildcards[rev]}"
popd
