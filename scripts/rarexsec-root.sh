#!/usr/bin/env bash

set -euo pipefail

SCRIPT_NAME="$(basename "$0")"

usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} <macro-file> [function-name]

Provide the macro path as the first argument and (optionally) the function name
without parentheses. For example:
  ${SCRIPT_NAME} analysis/event_display.C event_display_detector
EOF
}

if [[ ${#} -ge 1 && ( $1 == *"("* || $1 == *")"* ) ]]; then
    echo "Macro argument should not include a function call; pass the macro and function separately." >&2
    usage
    exit 1
fi

if [[ ${#} -ge 2 && ( $2 == *"("* || $2 == *")"* ) ]]; then
    echo "Function argument should be provided without parentheses." >&2
    usage
    exit 1
fi

TOPDIR="${RAREXSEC:-$(cd "$(dirname "$0")"/.. && pwd)}"
export RAREXSEC="$TOPDIR"

LIBDIR="$TOPDIR/build/lib"; [[ -d "$LIBDIR" ]] || LIBDIR="$TOPDIR/lib"
INCDIR="$TOPDIR/include"; [[ -d "$INCDIR/rarexsec" ]] || INCDIR="$TOPDIR/src"
SETUP="$TOPDIR/scripts/setup_rarexsec.C"; [[ -f "$SETUP" ]] || SETUP="$TOPDIR/setup_rarexsec.C"
LIB="$LIBDIR/librarexsec.so"
ROOT_OPTS=(-l -b -q)
while [[ $# -gt 0 && $1 == -* ]]; do
  ROOT_OPTS+=("$1")
  shift
done

MACRO="${1:-$TOPDIR/analysis/main.C}"
# If no function name is given, use the macro stem: foo.C -> foo()
FUNC="${2:-$(basename "${MACRO%.*}")}" 

: "${RAREXSEC_CFG:=$TOPDIR/data/samples.test-new-analysis.json}"
: "${RAREXSEC_CONFIG:=$RAREXSEC_CFG}"
: "${RAREXSEC_TREE:=events}"
: "${RAREXSEC_BEAMLINE:=numi-fhc}"
: "${RAREXSEC_PERIODS:=run1}"

export RAREXSEC_CFG RAREXSEC_CONFIG RAREXSEC_TREE RAREXSEC_BEAMLINE RAREXSEC_PERIODS

export LD_LIBRARY_PATH="$LIBDIR:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="$INCDIR:${ROOT_INCLUDE_PATH:-}"

root "${ROOT_OPTS[@]}" -e "gROOT->LoadMacro(\"$SETUP\"); setup_rarexsec(\"$LIB\",\"$INCDIR\"); gROOT->LoadMacro(\"$MACRO\"); ${FUNC}();"
