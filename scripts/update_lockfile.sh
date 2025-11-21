#!/usr/bin/env bash

set -e

pixi lock 2>&1
git add pixi.lock

exit 0