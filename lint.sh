#!/usr/bin/env bash
flake8 --max-line-length=88 --exit-zero automol
pylint --rcfile=.pylintrc automol
flake8 --max-line-length=88 --exit-zero phydat
pylint --rcfile=.pylintrc phydat
flake8 --max-line-length=88 --exit-zero autoreact
pylint --rcfile=.pylintrc autoreact
