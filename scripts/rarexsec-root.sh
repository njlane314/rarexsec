#!/usr/bin/env bash

set -euo pipefail

SCRIPT_NAME="$(basename "$0")"

usage() {
    cat <<EOF_USAGE
Usage: ${SCRIPT_NAME} <macro-file> [function-name]

Provide the macro path as the first argument and (optionally) the function name
without parentheses. For example:
  ${SCRIPT_NAME} analysis/event_display.C event_display_detector

You may also pass a combined ROOT-style argument, such as:
  ${SCRIPT_NAME} 'analysis/event_display.C("event_display_detector()")'
EOF_USAGE
}

if [[ ${#} -lt 1 ]]; then
    usage
    exit 1
fi

if [[ ${#} -ge 2 && ( $2 == *"("* || $2 == *")"* ) ]]; then
    echo "Function argument should be provided without parentheses." >&2
    usage
    exit 1
fi

macro_arg="${1}"; shift || true
func_arg="${1-}" || true

if [[ -z "$func_arg" && "$macro_arg" == *"("* && "$macro_arg" == *")"* ]]; then
    if [[ "$macro_arg" =~ ^(.+)\((.*)\)$ ]]; then
        macro_arg="${BASH_REMATCH[1]}"
        func_arg="${BASH_REMATCH[2]}"
        func_arg="${func_arg%\)}"
        func_arg="${func_arg#\(}"
        func_arg="${func_arg%\"}"
        func_arg="${func_arg#\"}"
        func_arg="${func_arg%\'}"; func_arg="${func_arg#\'}"
    else
        echo "Unable to parse combined macro/function argument." >&2
        usage
        exit 1
    fi
fi

TOPDIR="${RAREXSEC:-$(cd "$(dirname "$0")"/.. && pwd)}"
export RAREXSEC="$TOPDIR"

LIBDIR="$TOPDIR/build/lib"; [[ -d "$LIBDIR" ]] || LIBDIR="$TOPDIR/lib"
INCDIR="$TOPDIR/include"; [[ -d "$INCDIR/rarexsec" ]] || INCDIR="$TOPDIR/src"
SETUP="$TOPDIR/scripts/setup_rarexsec.C"; [[ -f "$SETUP" ]] || SETUP="$TOPDIR/setup_rarexsec.C"
LIB="$LIBDIR/librarexsec.so"
MACRO="${macro_arg:-$TOPDIR/analysis/main.C}"
# If no function name is given, use the macro stem: foo.C -> foo()
FUNC="${func_arg:-$(basename "${MACRO%.*}")}" 

: "${RAREXSEC_CFG:=$TOPDIR/data/samples.test-new-analysis.json}"
: "${RAREXSEC_CONFIG:=$RAREXSEC_CFG}"
: "${RAREXSEC_TREE:=events}"
: "${RAREXSEC_BEAMLINE:=numi-fhc}"
: "${RAREXSEC_PERIODS:=run1}"

export RAREXSEC_CFG RAREXSEC_CONFIG RAREXSEC_TREE RAREXSEC_BEAMLINE RAREXSEC_PERIODS

export LD_LIBRARY_PATH="$LIBDIR:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="$INCDIR:${ROOT_INCLUDE_PATH:-}"

root -l -b -q -e "gROOT->LoadMacro(\"$SETUP\"); setup_rarexsec(\"$LIB\",\"$INCDIR\"); gROOT->LoadMacro(\"$MACRO\"); ${FUNC}();"
