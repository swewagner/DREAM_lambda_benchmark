#!/usr/bin/env bash
set -e

BIN_FILE=$1
QUERY_FILE=$2
OUT_DIR=$3
LOG_FILE=$4
DUMMY_FILE=$5
BINARY_DIR="../build/iota"

$BINARY_DIR/ibf_magic --reference $BIN_FILE \
                        --query $QUERY_FILE \
                        --output_dir $OUT_DIR \
                        --kmer_size 12 \
                        --bin_threshold 0.7 \
                        protein \
                        --subject_domain auto \
                        --query_domain auto \
                        &> $LOG_FILE

touch $DUMMY_FILE
for file in $(ls $OUT_DIR/*.fasta);
do
  echo $file >> $DUMMY_FILE
done