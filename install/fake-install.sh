# run with: . /path/to/fake-install.sh
export THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export PROJECT_ROOT=$THIS_DIR/..
export PYTHONPATH=$PROJECT_ROOT:$PYTHONPATH
