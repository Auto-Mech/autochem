#!/usr/bin/env bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

(
    cd $SCRIPT_DIR
    echo "New linting:"
    bash ./lint.sh

    echo "Old linting:"
    pylint --rcfile=.pylintrc automol
    pylint --rcfile=.pylintrc phydat
    pylint --rcfile=.pylintrc autoreact
)
