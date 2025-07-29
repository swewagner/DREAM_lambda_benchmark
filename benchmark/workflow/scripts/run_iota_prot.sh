#!/usr/bin/env bash
set -e

BIN_FILE=$1
QUERY_FILE=$2
OUT_DIR=$3
K_MER=$4
MAX_ER=$5
RED_ALPH=$6
LOG_FILE=$7
DUMMY_FILE=$8
BINARY_DIR="../build/iota"

$BINARY_DIR/ibf_magic --reference $BIN_FILE \
                        --query $QUERY_FILE \
                        --output_dir $OUT_DIR \
                        --kmer_size $K_MER \
                        --max_error $MAX_ER \
                        protein \
                        --subject_domain auto \
                        --query_domain auto \
                        --reduced_alphabet $RED_ALPH \
                        &> $LOG_FILE

touch $DUMMY_FILE
for file in $(ls $OUT_DIR/*.fasta);
do
  echo $file >> $DUMMY_FILE
done