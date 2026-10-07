#!/usr/bin/env bash

function help() {
    cat << EOF
Usage: $0 [-clean] [-serve] [-report]

Without any flag run the watcher.ipynb notebook.

Options:
- clean    clear the output of the watcher.ipynb notebook.
- serve    start a local server for interactive use of notebooks.
- report   run the report.ipynb notebook and generate a report.
EOF
    exit 0
}

case "$1" in
    "")
        # Default action when no options are provided
        uv run jupyter nbconvert --execute --inplace watcher.ipynb
        ;;
    -clean)
        # Action for the clean flag
        uv run jupyter nbconvert --clear-output --inplace watcher.ipynb
        ;;
    -serve)
        # Just start a local server for interactive use
        uv run jupyter notebook --no-browser --ServerApp.token='' --ServerApp.password=''
        ;;
    -report)
        uv run jupyter nbconvert --execute --inplace report.ipynb
        uv run majordome-build-qmd --file report.ipynb
        ;;
    -h|--help)
        help
        ;;
    *)
        echo "Error: unknown option '$1'" >&2
        help
        exit 1
        ;;
esac
