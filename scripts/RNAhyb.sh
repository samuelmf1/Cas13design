#!/bin/bash
set -euo pipefail

# Batched replacement for the original per-pair RNAhybrid loop.
script_dir=$(cd -- "$(dirname -- "$0")" && pwd)
if [[ $# -ne 1 ]]; then
    echo "usage: $0 PAIRS.txt" >&2
    exit 2
fi
if [[ ! -x "${script_dir}/rnahybrid_paired" ]]; then
    echo "missing paired RNAhybrid executable: ${script_dir}/rnahybrid_paired" >&2
    exit 2
fi
if [[ ! -r "$1" ]]; then
    echo "pair input is not readable: $1" >&2
    exit 2
fi
exec "${script_dir}/rnahybrid_paired" < "$1"
