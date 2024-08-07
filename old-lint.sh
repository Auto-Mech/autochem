#!/usr/bin/env bash
set -e  # if any command fails, quit
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

(
    cd $SCRIPT_DIR
    echo "New linting:"
    source ./lint.sh

    echo "Old linting:"
    pylint --rcfile=.pylintrc automol
    pylint --rcfile=.pylintrc phydat
    pylint --rcfile=.pylintrc autoreact
)
