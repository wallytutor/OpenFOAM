#!/usr/bin/env bash
# ==============================================================================
# etc/test.sh - Test suite for OpenFOAM shell scripts under etc/rc/
#
# Comprehensive, edge-case test suite covering functions across:
#   - etc/rc/utilsrc   (formatting, environmental helpers, bootstrap)
#   - etc/rc/coresrc   (core detection, safe parallel core allocation)
#   - etc/rc/meshingrc (STL file inspection and renaming)
#
# Line length is strictly maintained at <= 80 characters.
# Tests execute in the parent shell with temporary redirection to ensure
# accurate tracking of test counters and side effects.
#
# Usage:
#   ./etc/test.sh
#   bash etc/test.sh
# ==============================================================================

# Script directory and repository root resolution
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"
RC_DIR="${SCRIPT_DIR}/rc"

# Create a temporary workspace isolated from system and repo files
TEST_TMP=$(mktemp -d "/tmp/openfoam_test_XXXXXX")
trap 'rm -rf "${TEST_TMP}"' EXIT

# Source the target rc files
source "${RC_DIR}/utilsrc"
source "${RC_DIR}/coresrc"
source "${RC_DIR}/meshingrc"

# ------------------------------------------------------------------------------
# Test Framework and Assertions
# ------------------------------------------------------------------------------

TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Print section banner
test_section() {
    printf "\n\033[1;34m=== %s ===\033[0m\n" "$1"
}

# Record a passing assertion
pass_test() {
    local desc="$1"
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    PASSED_TESTS=$((PASSED_TESTS + 1))
    printf "  \033[32m[PASS]\033[0m %s\n" "$desc"
}

# Record a failing assertion
fail_test() {
    local desc="$1"
    local detail="$2"
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    FAILED_TESTS=$((FAILED_TESTS + 1))
    printf "  \033[31m[FAIL]\033[0m %s\n" "$desc" >&2
    if [ -n "$detail" ]; then
        printf "         \033[33m%s\033[0m\n" "$detail" >&2
    fi
}

# Assert string equality
assert_eq() {
    local expected="$1"
    local actual="$2"
    local desc="$3"
    if [ "$expected" = "$actual" ]; then
        pass_test "$desc"
    else
        fail_test "$desc" "Expected: '${expected}', Got: '${actual}'"
    fi
}

# Assert exit status code
assert_status() {
    local expected="$1"
    local actual="$2"
    local desc="$3"
    if [ "$expected" -eq "$actual" ]; then
        pass_test "$desc"
    else
        fail_test "$desc" "Expected status: ${expected}, Got: ${actual}"
    fi
}

# Assert haystack string contains needle substring
assert_contains() {
    local needle="$1"
    local haystack="$2"
    local desc="$3"
    if [[ "$haystack" == *"$needle"* ]]; then
        pass_test "$desc"
    else
        fail_test "$desc" "Expected '${needle}' in output"
    fi
}

# Assert haystack string does not contain needle substring
assert_not_contains() {
    local needle="$1"
    local haystack="$2"
    local desc="$3"
    if [[ "$haystack" != *"$needle"* ]]; then
        pass_test "$desc"
    else
        fail_test "$desc" "Did not expect '${needle}' in output"
    fi
}

# Assert file exists on disk
assert_file_exists() {
    local path="$1"
    local desc="$2"
    if [ -f "$path" ]; then
        pass_test "$desc"
    else
        fail_test "$desc" "File '${path}' does not exist"
    fi
}

# Assert file does not exist on disk
assert_file_not_exists() {
    local path="$1"
    local desc="$2"
    if [ ! -f "$path" ]; then
        pass_test "$desc"
    else
        fail_test "$desc" "File '${path}' should not exist"
    fi
}

# Execute in the current shell with file redirection to avoid subshell
# isolation while capturing stdout, stderr, and exit status.
run_in_shell() {
    local out_file="${TEST_TMP}/run_out.tmp"
    local err_file="${TEST_TMP}/run_err.tmp"
    "$@" > "${out_file}" 2> "${err_file}"
    CMD_STATUS=$?
    CMD_OUT=$(<"${out_file}")
    CMD_ERR=$(<"${err_file}")
    /bin/rm -f "${out_file}" "${err_file}"
}

# ==============================================================================
# Test Suite: etc/rc/utilsrc - Terminal Message Functions
# ==============================================================================

test_section "Testing utilsrc: Terminal Message Functions"

# Test message_red: text content, color code, and reset code
run_in_shell message_red "Test red text"
assert_contains "Test red text" "$CMD_OUT" \
    "message_red outputs message text"
assert_contains $'\033[31m' "$CMD_OUT" \
    "message_red applies red color escape sequence"
assert_contains $'\033[0m' "$CMD_OUT" \
    "message_red resets terminal color sequence"

# Test message_yellow: text content, color code, and reset code
run_in_shell message_yellow "Test yellow text"
assert_contains "Test yellow text" "$CMD_OUT" \
    "message_yellow outputs message text"
assert_contains $'\033[33m' "$CMD_OUT" \
    "message_yellow applies yellow color escape sequence"
assert_contains $'\033[0m' "$CMD_OUT" \
    "message_yellow resets terminal color sequence"

# Test message_bold_green: text content, color code, and reset code
run_in_shell message_bold_green "Test bold green text"
assert_contains "Test bold green text" "$CMD_OUT" \
    "message_bold_green outputs message text"
assert_contains $'\033[1;32m' "$CMD_OUT" \
    "message_bold_green applies bold green escape sequence"
assert_contains $'\033[0m' "$CMD_OUT" \
    "message_bold_green resets terminal color sequence"

# Test message_bold_yellow: text content, color code, and reset code
run_in_shell message_bold_yellow "Test bold yellow text"
assert_contains "Test bold yellow text" "$CMD_OUT" \
    "message_bold_yellow outputs message text"
assert_contains $'\033[1;33m' "$CMD_OUT" \
    "message_bold_yellow applies bold yellow escape sequence"
assert_contains $'\033[0m' "$CMD_OUT" \
    "message_bold_yellow resets terminal color sequence"

# Test error_message: text content with prefix and red color
run_in_shell error_message "Something failed"
assert_contains "Error: Something failed" "$CMD_OUT" \
    "error_message prefixes message with 'Error:'"
assert_contains $'\033[31m' "$CMD_OUT" \
    "error_message applies red color escape sequence"

# Test warning_message: text content with prefix and yellow color
run_in_shell warning_message "Attention required"
assert_contains "Warning: Attention required" "$CMD_OUT" \
    "warning_message prefixes message with 'Warning:'"
assert_contains $'\033[33m' "$CMD_OUT" \
    "warning_message applies yellow color escape sequence"

# Test header_bold_green: multi-line banner formatting
run_in_shell header_bold_green "Header Title"
assert_contains "========================================" "$CMD_OUT" \
    "header_bold_green outputs top/bottom separator bars"
assert_contains "Header Title" "$CMD_OUT" \
    "header_bold_green outputs header title"
assert_contains $'\033[1;32m' "$CMD_OUT" \
    "header_bold_green applies bold green formatting"

# ==============================================================================
# Test Suite: etc/rc/utilsrc - Environment & Setup Helpers
# ==============================================================================

test_section "Testing utilsrc: Environment & Setup Helpers"

# Test handleOpenFOAM: detects existing directory and sources bashrc
SAVED_FOAM_TARGET="${FOAM_TARGET:-}"
SAVED_WM_PROJECT="${WM_PROJECT:-}"

mock_foam="${TEST_TMP}/mock_openfoam"
mkdir -p "${mock_foam}/etc"
echo 'export WM_PROJECT="OpenFOAM-Mock"' > "${mock_foam}/etc/bashrc"

FOAM_TARGET="${mock_foam}"
unset WM_PROJECT

run_in_shell handleOpenFOAM
assert_contains "OpenFOAM found at ${mock_foam}" "$CMD_OUT" \
    "handleOpenFOAM reports existing OpenFOAM installation"
assert_eq "OpenFOAM-Mock" "${WM_PROJECT:-}" \
    "handleOpenFOAM sources bashrc when WM_PROJECT unset"

# Test handleOpenFOAM: skips re-sourcing when WM_PROJECT is already set
mock_foam2="${TEST_TMP}/mock_openfoam2"
mkdir -p "${mock_foam2}/etc"
echo 'export WM_PROJECT="Overwritten"' > "${mock_foam2}/etc/bashrc"

FOAM_TARGET="${mock_foam2}"
export WM_PROJECT="AlreadySet"

run_in_shell handleOpenFOAM
assert_contains "OpenFOAM found at ${mock_foam2}" "$CMD_OUT" \
    "handleOpenFOAM reports installation when WM_PROJECT set"
assert_eq "AlreadySet" "${WM_PROJECT:-}" \
    "handleOpenFOAM preserves already active WM_PROJECT"

# Test handleOpenFOAM: installs OpenFOAM when FOAM_TARGET does not exist
mock_foam_install="${TEST_TMP}/mock_openfoam_installed"
FOAM_TARGET="${mock_foam_install}"
unset WM_PROJECT
APT_INSTALL_CALLED=false
APT_REPO_CALLED=false
sudo() {
    if [ "$1" = "add-apt-repository" ]; then
        APT_REPO_CALLED=true
    elif [ "$1" = "apt-get" ] && [ "$2" = "install" ]; then
        APT_INSTALL_CALLED=true
        mkdir -p "${mock_foam_install}/etc"
        echo 'export WM_PROJECT="OpenFOAM-AutoInstalled"' > \
            "${mock_foam_install}/etc/bashrc"
    fi
    return 0
}

run_in_shell handleOpenFOAM
assert_eq "true" "$APT_INSTALL_CALLED" \
    "handleOpenFOAM triggers installation when FOAM_TARGET missing"
assert_eq "true" "$APT_REPO_CALLED" \
    "handleOpenFOAM registers package repository when missing"
assert_eq "OpenFOAM-AutoInstalled" "${WM_PROJECT:-}" \
    "handleOpenFOAM sources installed OpenFOAM bashrc"

unset -f sudo
FOAM_TARGET="${SAVED_FOAM_TARGET}"
if [ -n "$SAVED_WM_PROJECT" ]; then
    export WM_PROJECT="${SAVED_WM_PROJECT}"
else
    unset WM_PROJECT
fi

# Test handleQuarto: recognizes installed quarto in environment
SAVED_PATH="$PATH"
mock_bin_quarto="${TEST_TMP}/mock_bin_quarto"
mkdir -p "${mock_bin_quarto}"
cat << 'EOF' > "${mock_bin_quarto}/quarto"
#!/bin/sh
if [ "$1" = "--version" ]; then
    echo "1.9.37"
    exit 0
fi
exit 0
EOF
chmod +x "${mock_bin_quarto}/quarto"

PATH="${mock_bin_quarto}:${SAVED_PATH}"
run_in_shell handleQuarto
assert_contains "Already installed: quarto 1.9.37" "$CMD_OUT" \
    "handleQuarto detects installed quarto binary version"
PATH="${SAVED_PATH}"

# Test handleQuarto: downloads and installs when quarto is missing
SAVED_DIR="$PWD"
mock_quarto_dir="${TEST_TMP}/mock_quarto_dir"
mkdir -p "${mock_quarto_dir}"
cd "${mock_quarto_dir}" || exit 1

QUARTO_WGET_CALLED=false
QUARTO_APT_CALLED=false
QUARTO_TINYTEX_CALLED=false

wget() {
    QUARTO_WGET_CALLED=true
    touch quarto-1.9.37-linux-amd64.deb
    return 0
}
sudo() {
    if [ "$1" = "apt-get" ] && [ "$2" = "install" ]; then
        QUARTO_APT_CALLED=true
    fi
    return 0
}
quarto() {
    if [ "$1" = "install" ] && [ "$2" = "tinytex" ]; then
        QUARTO_TINYTEX_CALLED=true
    fi
    return 0
}
command() {
    if [ "$1" = "-v" ] && [ "$2" = "quarto" ]; then
        return 1
    fi
    builtin command "$@"
}

run_in_shell handleQuarto

cd "${SAVED_DIR}" || exit 1
unset -f wget sudo quarto command

assert_eq "true" "$QUARTO_WGET_CALLED" \
    "handleQuarto downloads package deb when quarto is missing"
assert_eq "true" "$QUARTO_APT_CALLED" \
    "handleQuarto installs package deb when quarto is missing"
assert_eq "true" "$QUARTO_TINYTEX_CALLED" \
    "handleQuarto installs tinytex after package installation"

# Test handleUvInstall: recognizes installed uv in environment
mock_bin_uv="${TEST_TMP}/mock_bin_uv"
mkdir -p "${mock_bin_uv}"
cat << 'EOF' > "${mock_bin_uv}/uv"
#!/bin/sh
echo "uv 0.5.1"
EOF
chmod +x "${mock_bin_uv}/uv"

PATH="${mock_bin_uv}:${SAVED_PATH}"
run_in_shell handleUvInstall
assert_contains "Already installed: uv 0.5.1" "$CMD_OUT" \
    "handleUvInstall detects installed uv package manager"
PATH="${SAVED_PATH}"

# Test handleUvInstall: downloads and installs uv when missing
SAVED_HOME="$HOME"
mock_home="${TEST_TMP}/mock_home"
mkdir -p "${mock_home}/.local/bin"
echo 'export UV_INSTALLED="uv_active"' > "${mock_home}/.local/bin/env"

curl() {
    touch "${mock_home}/curl_called.tmp"
    return 0
}
sh() { return 0; }
command() {
    if [ "$1" = "-v" ] && [ "$2" = "uv" ]; then
        return 1
    fi
    builtin command "$@"
}

HOME="${mock_home}"
run_in_shell handleUvInstall
HOME="${SAVED_HOME}"

assert_file_exists "${mock_home}/curl_called.tmp" \
    "handleUvInstall invokes curl installer when uv is missing"
assert_eq "uv_active" "${UV_INSTALLED:-}" \
    "handleUvInstall sources environment file after installation"

unset -f command curl sh
unset UV_INSTALLED

# Test handleUvMode: non-WSL environment does not set UV_LINK_MODE
unset UV_LINK_MODE
grep() {
    if [ "$1" = "-qi" ] && [ "$2" = "microsoft" ]; then
        return 1
    fi
    command grep "$@"
}

run_in_shell handleUvMode
assert_eq "" "${UV_LINK_MODE:-}" \
    "handleUvMode does not set UV_LINK_MODE outside WSL"
unset -f grep

# Test handleUvMode: simulated WSL outside HOME sets UV_LINK_MODE=copy
grep() {
    if [ "$1" = "-qi" ] && [ "$2" = "microsoft" ]; then
        return 0
    fi
    command grep "$@"
}
realpath() {
    if [ "$1" = "$PWD" ]; then
        echo "/mnt/c/external/workspace"
    elif [ "$1" = "$HOME" ]; then
        echo "/home/user"
    else
        command realpath "$@"
    fi
}

unset UV_LINK_MODE
run_in_shell handleUvMode
assert_eq "copy" "${UV_LINK_MODE:-}" \
    "handleUvMode sets UV_LINK_MODE=copy in WSL outside HOME"

# Test handleUvMode: simulated WSL inside HOME preserves default
realpath() {
    if [ "$1" = "$PWD" ]; then
        echo "/home/user/workspace"
    elif [ "$1" = "$HOME" ]; then
        echo "/home/user"
    else
        command realpath "$@"
    fi
}

unset UV_LINK_MODE
run_in_shell handleUvMode
assert_eq "" "${UV_LINK_MODE:-}" \
    "handleUvMode preserves default UV_LINK_MODE in WSL inside HOME"

unset -f grep realpath
unset UV_LINK_MODE

# Test handleEnvironment: activates existing virtual environment
SAVED_DIR="$PWD"
SAVED_VENV="${VIRTUAL_ENV:-}"
mock_env_dir="${TEST_TMP}/test_venv_dir"
mkdir -p "${mock_env_dir}/.venv/bin"
echo 'export VIRTUAL_ENV_ACTIVATED="true"' > \
    "${mock_env_dir}/.venv/bin/activate"

cd "${mock_env_dir}" || exit 1
unset VIRTUAL_ENV
unset VIRTUAL_ENV_ACTIVATED

run_in_shell handleEnvironment
assert_contains "Activating existing virtual environment" "$CMD_OUT" \
    "handleEnvironment activates existing .venv"
assert_eq "true" "${VIRTUAL_ENV_ACTIVATED:-}" \
    "handleEnvironment sources the .venv/bin/activate script"

# Test handleEnvironment: deactivates conflicting virtual environment
cd "${mock_env_dir}" || exit 1
export VIRTUAL_ENV="/different/path/venv"
DEACTIVATE_CALLED=false
deactivate() {
    DEACTIVATE_CALLED=true
    unset VIRTUAL_ENV
}

run_in_shell handleEnvironment
assert_contains "Deactivating current virtual environment" "$CMD_OUT" \
    "handleEnvironment detects and deactivates conflicting venv"
assert_eq "true" "${DEACTIVATE_CALLED}" \
    "handleEnvironment calls deactivate function on conflict"

unset -f deactivate

# Test handleEnvironment: creates and initializes new virtual environment
mock_empty_env="${TEST_TMP}/test_new_venv_dir"
mkdir -p "${mock_empty_env}"
cd "${mock_empty_env}" || exit 1
unset VIRTUAL_ENV
unset NEW_VENV_CREATED
unset NEW_VENV_PACKAGES_INSTALLED

uv() {
    if [ "$1" = "venv" ]; then
        local target="${@: -1}"
        mkdir -p "${target}/bin"
        echo 'export NEW_VENV_CREATED="true"' > "${target}/bin/activate"
        return 0
    elif [ "$1" = "pip" ] && [ "$2" = "install" ]; then
        NEW_VENV_PACKAGES_INSTALLED="true"
        return 0
    fi
    return 1
}

run_in_shell handleEnvironment
assert_contains "Creating new virtual environment with uv" "$CMD_OUT" \
    "handleEnvironment creates new venv when .venv does not exist"
assert_eq "true" "${NEW_VENV_CREATED:-}" \
    "handleEnvironment invokes uv venv and activates it"
assert_eq "true" "${NEW_VENV_PACKAGES_INSTALLED:-}" \
    "handleEnvironment runs uv pip install for requirements.txt"

unset -f uv
cd "${SAVED_DIR}" || exit 1
if [ -n "$SAVED_VENV" ]; then
    export VIRTUAL_ENV="${SAVED_VENV}"
else
    unset VIRTUAL_ENV
fi

# Test bootstrapToolbox: calls top-level banner and component handlers
BOOTSTRAP_OPENFOAM_CALLED=false
BOOTSTRAP_QUARTO_CALLED=false
BOOTSTRAP_UV_INSTALL_CALLED=false
BOOTSTRAP_UV_MODE_CALLED=false

handleOpenFOAM()  { BOOTSTRAP_OPENFOAM_CALLED=true; }
handleQuarto()    { BOOTSTRAP_QUARTO_CALLED=true; }
handleUvInstall() { BOOTSTRAP_UV_INSTALL_CALLED=true; }
handleUvMode()    { BOOTSTRAP_UV_MODE_CALLED=true; }

run_in_shell bootstrapToolbox
assert_contains "Bootstraping tools..." "$CMD_OUT" \
    "bootstrapToolbox prints startup header banner"
assert_eq "true" "$BOOTSTRAP_OPENFOAM_CALLED" \
    "bootstrapToolbox invokes handleOpenFOAM"
assert_eq "true" "$BOOTSTRAP_QUARTO_CALLED" \
    "bootstrapToolbox invokes handleQuarto"
assert_eq "true" "$BOOTSTRAP_UV_INSTALL_CALLED" \
    "bootstrapToolbox invokes handleUvInstall"
assert_eq "true" "$BOOTSTRAP_UV_MODE_CALLED" \
    "bootstrapToolbox invokes handleUvMode"

# Restore original utilsrc functions
source "${RC_DIR}/utilsrc"

# ==============================================================================
# Test Suite: etc/rc/coresrc - Physical Core Count Detection
# ==============================================================================

test_section "Testing coresrc: Physical Core Count Detection"

# Real host hardware core detection test
run_in_shell get_physical_cores
assert_status 0 "$CMD_STATUS" \
    "get_physical_cores exits with 0 on host system"
if [[ "$CMD_OUT" =~ ^[0-9]+$ ]] && [ "$CMD_OUT" -gt 0 ]; then
    pass_test "get_physical_cores returns a positive integer (${CMD_OUT})"
else
    fail_test "get_physical_cores returns a positive integer" \
        "Got: '${CMD_OUT}'"
fi

# Edge case: lscpu command is missing from PATH
command() {
    if [ "$1" = "-v" ] && [ "$2" = "lscpu" ]; then
        return 1
    fi
    builtin command "$@"
}

run_in_shell get_physical_cores
unset -f command

assert_status 1 "$CMD_STATUS" \
    "get_physical_cores fails when lscpu is not available"
assert_contains "lscpu command is not found" "$CMD_ERR" \
    "get_physical_cores reports missing lscpu to stderr"

# Edge case: lscpu produces unparseable or comments-only output
lscpu() { echo "# comments only"; }
run_in_shell get_physical_cores
assert_status 1 "$CMD_STATUS" \
    "get_physical_cores fails when output parsing fails"
assert_contains "Failed to parse physical core count" "$CMD_ERR" \
    "get_physical_cores reports parse error to stderr"
unset -f lscpu

# Multi-socket hardware parsing verification
lscpu() {
    cat << 'EOF'
# Core,Socket
0,0
1,0
0,1
1,1
0,0
1,0
EOF
}
run_in_shell get_physical_cores
assert_status 0 "$CMD_STATUS" \
    "get_physical_cores parses multi-socket configuration"
assert_eq "4" "$CMD_OUT" \
    "get_physical_cores accurately counts unique socket/core pairs (4)"
unset -f lscpu

# ==============================================================================
# Test Suite: etc/rc/coresrc - safe_run_ncores Argument Validation
# ==============================================================================

test_section "Testing coresrc: safe_run_ncores Argument Validation"

# Missing arguments
run_in_shell safe_run_ncores
assert_status 1 "$CMD_STATUS" \
    "safe_run_ncores returns 1 when no arguments provided"
assert_contains "Missing REQUESTED_CORES argument" "$CMD_ERR" \
    "safe_run_ncores reports missing argument error"

# Non-integer primary argument
run_in_shell safe_run_ncores "abc"
assert_status 1 "$CMD_STATUS" \
    "safe_run_ncores rejects non-numeric core request"
assert_contains "must be positive integers" "$CMD_ERR" \
    "safe_run_ncores reports positive integer requirement"

# Zero primary argument
run_in_shell safe_run_ncores 0
assert_status 1 "$CMD_STATUS" \
    "safe_run_ncores rejects 0 cores requested"
assert_contains "must be positive integers" "$CMD_ERR" \
    "safe_run_ncores reports error for 0 requested cores"

# Negative primary argument
run_in_shell safe_run_ncores -4
assert_status 1 "$CMD_STATUS" \
    "safe_run_ncores rejects negative cores requested"
assert_contains "must be positive integers" "$CMD_ERR" \
    "safe_run_ncores reports error for negative requested cores"

# Non-integer fallback argument
run_in_shell safe_run_ncores 128 "invalid_fallback"
assert_status 1 "$CMD_STATUS" \
    "safe_run_ncores rejects non-numeric fallback argument"
assert_contains "must be positive integers" "$CMD_ERR" \
    "safe_run_ncores validates fallback arguments are integers"

# Zero fallback argument
run_in_shell safe_run_ncores 128 0
assert_status 1 "$CMD_STATUS" \
    "safe_run_ncores rejects 0 in fallback arguments"
assert_contains "must be positive integers" "$CMD_ERR" \
    "safe_run_ncores validates fallbacks are strictly positive"

# Negative fallback argument
run_in_shell safe_run_ncores 128 -8
assert_status 1 "$CMD_STATUS" \
    "safe_run_ncores rejects negative values in fallback arguments"
assert_contains "must be positive integers" "$CMD_ERR" \
    "safe_run_ncores validates fallbacks cannot be negative"

# ==============================================================================
# Test Suite: etc/rc/coresrc - safe_run_ncores Fallback Allocation
# ==============================================================================

test_section "Testing coresrc: safe_run_ncores Fallback Allocation"

# Mock get_physical_cores = 20 to test deterministic allocation scenarios
get_physical_cores() { echo "20"; }

# Case 1: Primary requested cores is feasible (8 <= 20)
run_in_shell safe_run_ncores 8
assert_status 0 "$CMD_STATUS" \
    "safe_run_ncores returns 0 when requested is feasible"
assert_eq "8" "$CMD_OUT" \
    "safe_run_ncores returns requested cores (8) directly"
assert_eq "" "$CMD_ERR" \
    "safe_run_ncores produces no stderr warning when feasible"

# Case 2: Primary requested cores exactly equals physical cores (20 == 20)
run_in_shell safe_run_ncores 20
assert_status 0 "$CMD_STATUS" \
    "safe_run_ncores returns 0 when requested equals physical cores"
assert_eq "20" "$CMD_OUT" \
    "safe_run_ncores returns physical cores (20) directly"
assert_eq "" "$CMD_ERR" \
    "safe_run_ncores produces no stderr warning when exact match"

# Case 3: Primary requested is feasible even with fallbacks supplied
run_in_shell safe_run_ncores 16 8 4
assert_status 0 "$CMD_STATUS" \
    "safe_run_ncores succeeds with primary cores when fallbacks given"
assert_eq "16" "$CMD_OUT" \
    "safe_run_ncores selects primary 16 without trying fallbacks"
assert_eq "" "$CMD_ERR" \
    "safe_run_ncores does not print fallback warnings when primary fits"

# Case 4: Exceeds physical cores, no fallbacks provided -> half cores (10)
run_in_shell safe_run_ncores 128
assert_status 0 "$CMD_STATUS" \
    "safe_run_ncores succeeds when falling back to half cores"
assert_eq "10" "$CMD_OUT" \
    "safe_run_ncores 128 on 20-core host returns half cores (10)"
assert_contains "Requested cores (128) exceed physical cores (20)" \
    "$CMD_ERR" \
    "safe_run_ncores warns that requested cores exceed physical count"
assert_contains "Using half of physical cores (10)" "$CMD_ERR" \
    "safe_run_ncores notifies fallback to half cores"

# Case 5: Exceeds physical cores, single fallback also exceeds -> half (10)
run_in_shell safe_run_ncores 128 48
assert_status 0 "$CMD_STATUS" \
    "safe_run_ncores succeeds when all options exceed physical cores"
assert_eq "10" "$CMD_OUT" \
    "safe_run_ncores 128 48 on 20-core host returns half cores (10)"
assert_contains "Requested cores (128) exceed physical cores (20)" \
    "$CMD_ERR" \
    "safe_run_ncores reports primary request exceeded"
assert_contains "Fallback cores (48) exceed physical cores (20)" \
    "$CMD_ERR" \
    "safe_run_ncores reports fallback cores exceeded"
assert_contains "Using half of physical cores (10)" "$CMD_ERR" \
    "safe_run_ncores reports fallback to half cores after failed fallback"

# Case 6: Fallback chain where middle fallback fits (128 64 32 16 on 20 -> 16)
run_in_shell safe_run_ncores 128 64 32 16
assert_status 0 "$CMD_STATUS" \
    "safe_run_ncores succeeds with valid fallback in chain"
assert_eq "16" "$CMD_OUT" \
    "safe_run_ncores 128 64 32 16 on 20-core host returns 16"
assert_contains "Fallback cores (64) exceed physical cores (20)" \
    "$CMD_ERR" \
    "safe_run_ncores logs exceeded fallback 64"
assert_contains "Fallback cores (32) exceed physical cores (20)" \
    "$CMD_ERR" \
    "safe_run_ncores logs exceeded fallback 32"
assert_contains "Using fallback cores (16)" "$CMD_ERR" \
    "safe_run_ncores logs selection of feasible fallback 16"
assert_not_contains "half of physical cores" "$CMD_ERR" \
    "safe_run_ncores does not use half cores when a fallback succeeds"

# Case 7: Fallback chain where first fallback fits (128 16 8 -> 16)
run_in_shell safe_run_ncores 128 16 8
assert_status 0 "$CMD_STATUS" \
    "safe_run_ncores succeeds on first fallback in chain"
assert_eq "16" "$CMD_OUT" \
    "safe_run_ncores selects first feasible fallback (16)"
assert_not_contains "Fallback cores (8)" "$CMD_ERR" \
    "safe_run_ncores does not evaluate fallbacks after a match"

# ==============================================================================
# Test Suite: etc/rc/coresrc - safe_run_ncores Edge Cases
# ==============================================================================

test_section "Testing coresrc: safe_run_ncores Edge Cases (Odd & 1 Core)"

# Odd physical cores count: 7 cores -> integer half is 3
get_physical_cores() { echo "7"; }
run_in_shell safe_run_ncores 64
assert_status 0 "$CMD_STATUS" \
    "safe_run_ncores succeeds on odd physical cores"
assert_eq "3" "$CMD_OUT" \
    "safe_run_ncores calculates integer half of 7 as 3"

# Single physical core: 1 core -> half is 0 clamped to 1 minimum
get_physical_cores() { echo "1"; }
run_in_shell safe_run_ncores 16
assert_status 0 "$CMD_STATUS" \
    "safe_run_ncores succeeds on single-core host"
assert_eq "1" "$CMD_OUT" \
    "safe_run_ncores ensures minimum core fallback is at least 1"

# Upstream failure: get_physical_cores returns an error
get_physical_cores() { return 1; }
run_in_shell safe_run_ncores 8
assert_status 1 "$CMD_STATUS" \
    "safe_run_ncores returns error when get_physical_cores fails"

# Restore genuine get_physical_cores
source "${RC_DIR}/coresrc"

# Integration check on the actual hardware host
run_in_shell safe_run_ncores 1
assert_status 0 "$CMD_STATUS" \
    "safe_run_ncores 1 succeeds on host hardware"
assert_eq "1" "$CMD_OUT" \
    "safe_run_ncores 1 returns 1 on host hardware"

# ==============================================================================
# Test Suite: etc/rc/meshingrc - Help and Argument Handling
# ==============================================================================

test_section "Testing meshingrc: Help and Argument Handling"

# Help option -h
run_in_shell rename_stl_files_from_solid_names -h
assert_status 0 "$CMD_STATUS" \
    "rename_stl_files_from_solid_names -h exits with status 0"
assert_contains "Usage: rename_stl_files_from_solid_names" "$CMD_OUT" \
    "rename_stl_files_from_solid_names -h prints usage guide"

# Help option --help
run_in_shell rename_stl_files_from_solid_names --help
assert_status 0 "$CMD_STATUS" \
    "rename_stl_files_from_solid_names --help exits with status 0"
assert_contains "Usage: rename_stl_files_from_solid_names" "$CMD_OUT" \
    "rename_stl_files_from_solid_names --help prints usage guide"

# Non-existent target directory
run_in_shell rename_stl_files_from_solid_names "/nonexistent/path/for/stl"
assert_status 1 "$CMD_STATUS" \
    "rename_stl_files_from_solid_names fails on invalid directory"
assert_contains "Cannot access directory" "$CMD_OUT" \
    "rename_stl_files_from_solid_names reports directory error"

# ==============================================================================
# Test Suite: etc/rc/meshingrc - Empty Directory & Dry-run Mode
# ==============================================================================

test_section "Testing meshingrc: Empty Directory & Dry-run Mode"

# Empty directory handling
empty_stl_dir="${TEST_TMP}/empty_stl_dir"
mkdir -p "${empty_stl_dir}"

run_in_shell rename_stl_files_from_solid_names "${empty_stl_dir}"
assert_status 0 "$CMD_STATUS" \
    "rename_stl_files_from_solid_names succeeds on empty directory"
assert_contains "Renamed: 0" "$CMD_OUT" \
    "rename_stl_files_from_solid_names reports 0 renamed files"
assert_contains "Skipped: 0" "$CMD_OUT" \
    "rename_stl_files_from_solid_names reports 0 skipped files"
assert_contains "Warnings: 0" "$CMD_OUT" \
    "rename_stl_files_from_solid_names reports 0 warnings"

# Non-STL files in directory ignored
non_stl_dir="${TEST_TMP}/non_stl_dir"
mkdir -p "${non_stl_dir}"
touch "${non_stl_dir}/mesh.obj" "${non_stl_dir}/notes.txt"

run_in_shell rename_stl_files_from_solid_names "${non_stl_dir}"
assert_status 0 "$CMD_STATUS" \
    "rename_stl_files_from_solid_names ignores non-stl files"
assert_contains "Renamed: 0" "$CMD_OUT" \
    "Non-stl directory reports 0 renames"

# Dry run mode (-n flag)
dry_dir="${TEST_TMP}/dry_run_stl_dir"
mkdir -p "${dry_dir}"

cat << 'EOF' > "${dry_dir}/tunnel(1).stl"
solid "environment-bottom"
  facet normal 0 0 1
  endfacet
endsolid "environment-bottom"
EOF

run_in_shell rename_stl_files_from_solid_names -n "${dry_dir}"
assert_status 0 "$CMD_STATUS" \
    "rename_stl_files_from_solid_names -n exits with 0"
assert_contains "DRY RUN MODE" "$CMD_OUT" \
    "rename_stl_files_from_solid_names reports DRY RUN MODE"
assert_contains \
    "[DRY RUN] Would rename: tunnel(1).stl -> environment-bottom.stl" \
    "$CMD_OUT" \
    "rename_stl_files_from_solid_names logs planned rename in dry run"
assert_contains "Would rename: 1" "$CMD_OUT" \
    "rename_stl_files_from_solid_names reports correct dry-run count"

# Verify original file was NOT renamed during dry run
assert_file_exists "${dry_dir}/tunnel(1).stl" \
    "Original file is preserved during dry run"
assert_file_not_exists "${dry_dir}/environment-bottom.stl" \
    "Target file is not created during dry run"

# Dry run mode (--dry-run long flag)
run_in_shell rename_stl_files_from_solid_names --dry-run "${dry_dir}"
assert_status 0 "$CMD_STATUS" \
    "rename_stl_files_from_solid_names --dry-run exits with 0"
assert_contains "Would rename: 1" "$CMD_OUT" \
    "rename_stl_files_from_solid_names --dry-run reports rename count"

# ==============================================================================
# Test Suite: etc/rc/meshingrc - File Renaming, Skipping & Warnings
# ==============================================================================

test_section "Testing meshingrc: File Renaming, Skipping & Warnings"

# Full file renaming workflow with mixed cases
work_dir="${TEST_TMP}/active_stl_dir"
mkdir -p "${work_dir}"

# File 1: Needs renaming (tunnel(1).stl -> wall.stl)
cat << 'EOF' > "${work_dir}/tunnel(1).stl"
solid "wall"
  facet normal 0 0 1
  endfacet
endsolid "wall"
EOF

# File 2: Already has matching solid name (inlet.stl -> inlet)
cat << 'EOF' > "${work_dir}/inlet.stl"
solid "inlet"
  facet normal 0 1 0
  endfacet
endsolid "inlet"
EOF

# File 3: Missing quotes around solid name (unparseable solid name)
cat << 'EOF' > "${work_dir}/broken.stl"
solid unquoted_solid_name
  facet normal 1 0 0
  endfacet
endsolid unquoted_solid_name
EOF

# File 4: Solid name containing spaces ("base plate")
cat << 'EOF' > "${work_dir}/plate.stl"
solid "base plate"
  facet normal 0 0 -1
  endfacet
endsolid "base plate"
EOF

run_in_shell rename_stl_files_from_solid_names "${work_dir}"
assert_status 0 "$CMD_STATUS" \
    "rename_stl_files_from_solid_names executes rename batch with status 0"

# Verify rename operations
assert_contains "Renaming: tunnel(1).stl -> wall.stl" "$CMD_OUT" \
    "rename_stl_files_from_solid_names logs active rename"
assert_file_exists "${work_dir}/wall.stl" \
    "Target file 'wall.stl' was created"
assert_file_not_exists "${work_dir}/tunnel(1).stl" \
    "Source file 'tunnel(1).stl' was moved"

# Verify solid name with space
assert_contains "Renaming: plate.stl -> base plate.stl" "$CMD_OUT" \
    "rename_stl_files_from_solid_names handles solid name with spaces"
assert_file_exists "${work_dir}/base plate.stl" \
    "Target file 'base plate.stl' was created"
assert_file_not_exists "${work_dir}/plate.stl" \
    "Source file 'plate.stl' was moved"

# Verify skipped file
assert_contains "Skipping: inlet.stl" "$CMD_OUT" \
    "rename_stl_files_from_solid_names skips already matching name"
assert_file_exists "${work_dir}/inlet.stl" \
    "Skipped file 'inlet.stl' remains untouched"

# Verify warning on unparseable solid name
assert_contains "Warning: No solid name found in broken.stl" \
    "$CMD_OUT" \
    "rename_stl_files_from_solid_names warns about missing solid name"
assert_file_exists "${work_dir}/broken.stl" \
    "Unparseable file 'broken.stl' remains untouched"

# Verify summary counts
assert_contains "Renamed: 2" "$CMD_OUT" \
    "Summary reports 2 renamed files"
assert_contains "Skipped: 1" "$CMD_OUT" \
    "Summary reports 1 skipped file"
assert_contains "Warnings: 1" "$CMD_OUT" \
    "Summary reports 1 warning"

# Test execution in current directory (directory argument omitted)
curr_dir="${TEST_TMP}/current_dir_stl"
mkdir -p "${curr_dir}"

cat << 'EOF' > "${curr_dir}/temp.stl"
solid "outlet"
endfacet
endsolid "outlet"
EOF

cd "${curr_dir}" || exit 1
run_in_shell rename_stl_files_from_solid_names
cd "${REPO_ROOT}" || exit 1

assert_status 0 "$CMD_STATUS" \
    "rename_stl_files_from_solid_names succeeds when directory omitted"
assert_file_exists "${curr_dir}/outlet.stl" \
    "rename_stl_files_from_solid_names renamed temp.stl in current dir"
assert_file_not_exists "${curr_dir}/temp.stl" \
    "rename_stl_files_from_solid_names moved temp.stl in current dir"

# ==============================================================================
# Final Test Summary and Exit
# ==============================================================================

DIV="================================================================"
printf "\n\033[1;36m%s\033[0m\n" "$DIV"
printf "\033[1;36mTest Summary: Total: %d | Passed: %d | Failed: %d\033[0m\n" \
    "$TOTAL_TESTS" "$PASSED_TESTS" "$FAILED_TESTS"
printf "\033[1;36m%s\033[0m\n" "$DIV"

if [ "$FAILED_TESTS" -eq 0 ]; then
    printf "\033[1;32mALL TESTS PASSED SUCCESSFULLY!\033[0m\n\n"
    exit 0
else
    printf "\033[1;31mSOME TESTS FAILED! (%d failures)\033[0m\n\n" \
        "$FAILED_TESTS" >&2
    exit 1
fi
