#!/usr/bin/env bash

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
    *)
        echo "Error: unknown option '$1'" >&2
        exit 1
        ;;
esac
