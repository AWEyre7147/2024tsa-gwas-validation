#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
task_root=$(pwd)
task_conda="${CONDA_EXE:-$HOME/anaconda3/bin/conda}"
task_prefix="$HOME/anaconda3/envs/tsa2024-association"
if [ ! -d "$task_prefix/conda-meta" ]; then
  "$task_conda" create -y -p "$task_prefix" --file scripts/association-linux-64.lock.txt
fi
"$task_conda" run --no-capture-output -p "$task_prefix" Rscript scripts/install_gmmat.R "$task_root"
"$task_conda" list -p "$task_prefix" --explicit > data/tools/association/conda-explicit.private.txt
