#!/usr/bin/env bash
# ruff check automol
# ruff check phydat
# ruff check autoreact
pylint --rcfile=.pylintrc automol
pylint --rcfile=.pylintrc phydat
pylint --rcfile=.pylintrc autoreact

