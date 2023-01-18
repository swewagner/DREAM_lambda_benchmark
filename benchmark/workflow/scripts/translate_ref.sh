#!/usr/bin/env bash
set -e

IN_FILE=$1
OUT_FILE=$2
BINARY_DIR="/home/darklyght/Documents/Studium/Bioinformatik/Master/Praddidum/Benchmarking/DREAM_lambda_benchmark/lib/lambda_data_simulation/build"

$BINARY_DIR/translate_ref --in $IN_FILE \
                          --out $OUT_FILE \
                          &> /dev/null