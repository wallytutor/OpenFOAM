#!/usr/bin/env bash

function help() {
    cat << EOF
Usage: $0 [-clean|-serve|-report|<dirname>]

Without any flag run the watcher.ipynb notebook for the latest case.
If a case directory is provided as an argument, use that case.
Only a single option can be provided at any time.

Options:
- clean    clear the output of the watcher.ipynb notebook.
- serve    start a local server for interactive use of notebooks.
- report   run the report.ipynb notebook and generate a report.
EOF
    exit 0
}

function clear_output() {
    uv run jupyter nbconvert \
        --ClearMetadataPreprocessor.enabled=True \
        --ClearMetadataPreprocessor.clear_cell_metadata=True \
        --clear-output --inplace watcher.ipynb
}

function start_server() {
    uv run jupyter notebook --no-browser \
        --ServerApp.token='' --ServerApp.password='' &

    echo "Started jupyter notebook server in background."
    echo "Press Ctrl+C to stop."
    trap 'kill $!' SIGINT SIGTERM
    wait $!
}

function generate_report() {
    if [ ! -d "solution.1" ] || [ ! -d "solution.2" ]; then
        echo "Error: solution.1 or solution.2 not found."
        exit 1
    fi

    uv run jupyter nbconvert --execute --inplace report.ipynb
    uv run majordome-build-qmd --file report.ipynb
}

function run_watcher() {
    if [ $# -eq 1 ] && [ -d "$1" ]; then
        case_name="$1"
    else
        case_name=$(command ls -1d solution.* | tail -n -1)
        echo "Using latest available case: $case_name"
    fi

    CASE_NAME=$case_name uv run jupyter nbconvert --execute --inplace watcher.ipynb
}

function main() {
    if [ $# -ne 1 ] && [ $# -ne 0 ]; then
        help
    fi

    case "$1" in
        -clean)
            clear_output
            ;;
        -serve)
            start_server
            ;;
        -report)
            generate_report
            ;;
        -h|--help)
            help
            ;;
        *)
            run_watcher "$1"
            ;;
    esac
}

main "$@"
