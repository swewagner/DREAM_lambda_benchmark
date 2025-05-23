#!/usr/bin/env bash
set -e

BIN_FILE=$1
QUERY_FILE=$2
OUT_DIR=$3
LOG_FILE=$4
DUMMY_FILE=$5
BINARY_DIR="../build/iota"

$BINARY_DIR/ibf_magic -r $BIN_FILE \
                        -q $QUERY_FILE \
                        -o $OUT_DIR \
                        nucleotide \
                        -M 2 \
                        &> $LOG_FILE

touch $DUMMY_FILE