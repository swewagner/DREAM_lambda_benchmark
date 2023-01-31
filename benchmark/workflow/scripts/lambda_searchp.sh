#!/usr/bin/env bash
set -e

QUERY_FILE=$1
INDEX_FILE=$2
OUT_FILE=$3
BINARY_DIR="../build/lambda/src/lambda3-build/bin"

$BINARY_DIR/lambda3 searchp -q $QUERY_FILE -i $INDEX_FILE -o $OUT_FILE &> /dev/null