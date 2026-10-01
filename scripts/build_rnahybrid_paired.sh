#!/bin/bash
set -euo pipefail

# Builds the batched RNAhybrid replacement described in
# ../RNAHYBRID_BATCH_SPEEDUP.md. Requires the RNAhybrid-2.1.2 source to be
# unpacked alongside this script first (see Install.txt).

script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
source_dir=${1:-${script_dir}/RNAhybrid-2.1.2}
output=${2:-${script_dir}/rnahybrid_paired}
paired_source=${3:-${script_dir}/rnahybrid_paired.c}
source_root=$(cd -- "$source_dir" && pwd)
source_files=${source_root}/src

if [[ ! -f "${source_root}/config.h" || ! -f "${source_files}/hybrid_core.c" ]]; then
    echo "invalid RNAhybrid source directory: $source_dir" >&2
    exit 2
fi
if [[ ! -f "$paired_source" ]]; then
    echo "missing paired RNAhybrid source: $paired_source" >&2
    exit 2
fi

output_dir=$(dirname -- "$output")
mkdir -p "$output_dir"
output_dir=$(cd -- "$output_dir" && pwd)
output=${output_dir}/$(basename -- "$output")
temporary=${output}.tmp.$$
trap 'rm -f -- "$temporary"' EXIT

${CC:-gcc} \
    -DHAVE_CONFIG_H \
    -I"${source_root}" \
    -I"${source_files}" \
    -O2 \
    -fcommon \
    -o "$temporary" \
    "$paired_source" \
    "${source_files}/hybrid_core.c" \
    "${source_files}/energy.c" \
    "${source_files}/input.c" \
    "${source_files}/plot.c" \
    -lm

chmod 755 "$temporary"
mv -f -- "$temporary" "$output"
trap - EXIT
printf 'Built %s\n' "$output"
