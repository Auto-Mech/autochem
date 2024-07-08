#!/usr/bin/env bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
FILES=(
    "automol/inchi_key.py"
    "automol/error.py"
    "automol/vmat.py"
)

(
    cd $SCRIPT_DIR
    echo "New linting:"
    pre-commit run black --files ${FILES[@]}
    pre-commit run ruff --files ${FILES[@]}
    pre-commit run mypy --files ${FILES[@]}

    echo "Old linting:"
    pylint --rcfile=.pylintrc automol
    pylint --rcfile=.pylintrc phydat
    pylint --rcfile=.pylintrc autoreact
)
