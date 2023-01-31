#!/usr/bin/env bash
set -e

IN_FILE=$1
OUT_FILE=$2
BINARY_DIR="../build/lambda/src/lambda3-build/bin"

$BINARY_DIR/lambda3 mkindexp -d $IN_FILE -i $OUT_FILE &> /dev/null