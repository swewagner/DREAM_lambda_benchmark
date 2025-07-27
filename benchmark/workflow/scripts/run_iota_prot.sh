#!/usr/bin/env bash
set -e

BIN_FILE=$1
QUERY_FILE=$2
OUT_DIR=$3
K_MER=$4
MAX_ER=$5
LOG_FILE=$6
DUMMY_FILE=$7
BINARY_DIR="../build/iota"

$BINARY_DIR/ibf_magic --reference $BIN_FILE \
                        --query $QUERY_FILE \
                        --output_dir $OUT_DIR \
                        --kmer_size $K_MER \
                        --max_error $MAX_ER \
                        -t 0.2 \
                        protein \
                        --subject_domain auto \
                        --query_domain auto \
                        &> $LOG_FILE

touch $DUMMY_FILE
for file in $(ls $OUT_DIR/*.fasta);
do
  echo $file >> $DUMMY_FILE
done