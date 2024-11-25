#!/usr/bin/env bash

exec 2> "${snakemake_log[0]}"  # send all stderr from this script to the log file

cp "${snakemake_input[env]}" "${snakemake_output[env]}"

cat <<EOF >> "${snakemake_output[post]}"
#!/usr/bin/env bash

pip install --editable "resources/fipy~${snakemake_wildcards[rev]}/repo"
EOF
