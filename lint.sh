#!/usr/bin/env bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
FILES=(
    # "automol/inchi_key.py"
    # "automol/error.py"
    # "automol/util/_util.py"
    # "automol/util/__init__.py"
    # "automol/util/heuristic.py"
    # "automol/util/matrix.py"
    # "automol/util/ring.py"
    # "automol/util/tensor.py"
    # "automol/util/vector.py"
    # "automol/util/zmat.py"
    # "automol/util/dict_/_dict_.py"
    # "automol/util/dict_/multi.py"
    # "automol/form/_form.py"
    # "automol/form/reac.py"
    # "automol/mult/_mult.py"
    # "automol/mult/ts.py"
    # "automol/const.py"
    # "automol/vmat.py"
    # "automol/prop/_wfn.py"
    # "automol/prop/freq.py"
    # "automol/embed/_cleanup.py"
    # "automol/embed/_dgeom.py"
    "autoreact/*"
    "autochem/*"
)

(
    cd $SCRIPT_DIR
    pre-commit run black --files ${FILES[@]}
    pre-commit run ruff --files ${FILES[@]}
    pre-commit run mypy --files ${FILES[@]}
)
