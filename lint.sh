#!/usr/bin/env bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
FILES=(
    "automol/inchi_key.py"
    "automol/error.py"
    "automol/util/_util.py"
    "automol/util/__init__.py"
    "automol/util/heuristic.py"
    "automol/util/matrix.py"
    "automol/util/ring.py"
    "automol/util/tensor.py"
    "automol/util/vector.py"
)

(
    cd $SCRIPT_DIR
    pre-commit run black --files ${FILES[@]}
    pre-commit run ruff --files ${FILES[@]}
    pre-commit run mypy --files ${FILES[@]}
)
