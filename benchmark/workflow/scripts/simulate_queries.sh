#!/usr/bin/env bash
set -e

IN_FILE=$1
OUT_FILE=$2
TRUTH_FILE=$3
NUM=$4
LEN=$5
R_NUM=$6
BINARY_DIR="../lib/lambda_data_simulation/build/bin"

$BINARY_DIR/generate_query_matches --in $IN_FILE \
                                   --out $OUT_FILE \
                                   --ground_truth $TRUTH_FILE \
                                   --num_of_queries $NUM \
                                   --len_of_queries $LEN \
                                   --num_of_references $R_NUM \
                                   &> /dev/null
