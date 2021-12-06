#!/usr/bin/env bash
flake8 --exit-zero automol
pylint --rcfile=.pylintrc automol
flake8 --exit-zero phydat
pylint --rcfile=.pylintrc phydat
flake8 --exit-zero autoreact
pylint --rcfile=.pylintrc autoreact
