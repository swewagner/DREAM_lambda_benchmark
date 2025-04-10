#!/usr/bin/env bash
set -e

IN_FILE=$1
OUT_FILE=$2
TRUTH_FILE=$3
NUM=$4
LEN=$5
R_NUM=$6
ERR_RATE=$7
BINARY_DIR="../build/bin"

$BINARY_DIR/generate_query_matches --in $IN_FILE \
                                   --out $OUT_FILE \
                                   --ground_truth $TRUTH_FILE \
                                   --num_of_queries $NUM \
                                   --len_of_match $LEN \
                                   --num_of_references $R_NUM \
                                   --max_error_rate $ERR_RATE\
                                   #&> /dev/null
