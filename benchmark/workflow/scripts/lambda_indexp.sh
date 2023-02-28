#!/usr/bin/env bash
set -e

IN_FILE=$1
OUT_FILE=$2
LOG_FILE=$3
BINARY_DIR="../build/lambda3/src/lambda3-build/bin"

$BINARY_DIR/lambda3 mkindexp -d $IN_FILE -i $OUT_FILE &> $LOG_FILE