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
                        protein \
                        -S auto \
                        -Q auto \
                        &> $LOG_FILE

touch $DUMMY_FILE
for file in $(ls $OUT_DIR/*.fasta);
do
  echo $file >> $DUMMY_FILE
done