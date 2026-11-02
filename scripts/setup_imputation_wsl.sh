#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)
export CONDA_PKGS_DIRS="$HOME/.cache/tsa2024-imputation-conda"
mkdir -p "$ROOT/data/02-tools" "$ROOT/data/02-reference/dog10k"
TASK_CONDA="${CONDA_EXE:-$HOME/anaconda3/bin/conda}"
TASK_PREFIX="$HOME/anaconda3/envs/tsa2024-imputation"
"$TASK_CONDA" create -y --override-channels -c conda-forge -c bioconda \
  -p "$TASK_PREFIX" python=3.12 pip bcftools samtools openjdk=17 ucsc-liftover
"$TASK_PREFIX/bin/python" -m pip install -r "$ROOT/scripts/requirements-imputation.txt"
"$TASK_CONDA" list -p "$TASK_PREFIX" --explicit > "$ROOT/data/02-tools/imputation-wsl-explicit.txt"
