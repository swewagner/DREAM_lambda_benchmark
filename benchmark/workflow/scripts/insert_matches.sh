#!/usr/bin/env bash
set -e

IN_FILE_Q=$1
IN_FILE_M=$2
OUT_FILE=$3
BINARY_DIR="../build/bin"

$BINARY_DIR/insert_query_matches --in_queries $IN_FILE_Q \
                                 --in_matches $IN_FILE_M \
                                 --out $OUT_FILE \
                                   &> /dev/null