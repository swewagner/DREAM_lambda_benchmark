#!/usr/bin/env bash
#set -ex

IN_FILE_Q=$1
IN_FILE_R=$2
OUT_FILE_G=$3
OUT_DIR=$4
MATCH_LEN=$5
R_NUM=$6
B_NUM=$7
ERR_RATE=$8
LOG_FILE=$9
BINARY_DIR="../build/bin"
mkdir -p $OUT_DIR

$BINARY_DIR/distribute_matches_prot --in_queries $IN_FILE_Q \
                               --in_refs $IN_FILE_R \
                               --out $OUT_DIR \
                               --ground_truth $OUT_FILE_G \
                               --bin_num $B_NUM \
                               --len_of_match $MATCH_LEN \
                               --num_of_references $R_NUM \
                               --max_error_rate $ERR_RATE \
                               --verbose-ids
                               &> $LOG_FILE