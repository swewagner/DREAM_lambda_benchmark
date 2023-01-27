#!/usr/bin/env bash
set -e

IN_FILE=$1
OUT_FILE=$2
BINARY_DIR="../lib/lambda_data_simulation/build/bin"

$BINARY_DIR/translate_ref --in $IN_FILE \
                          --out $OUT_FILE \
                          &> /dev/null