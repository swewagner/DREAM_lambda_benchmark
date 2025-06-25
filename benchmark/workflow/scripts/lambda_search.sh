#!/usr/bin/env bash
set -e

QUERY_FILE=$1
INDEX_FILE=$2
OUT_FILE=$3
BLAST_M=$4
LOG_FILE=$5
BINARY_DIR="../build/lambda3/src/lambda3-build/bin"

if [ $BLAST_M = "blastN" ]; then
$BINARY_DIR/lambda3 searchn -q $QUERY_FILE -i $INDEX_FILE -o $OUT_FILE &> $LOG_FILE
else
$BINARY_DIR/lambda3 searchp -q $QUERY_FILE -i $INDEX_FILE -o $OUT_FILE &> $LOG_FILE
fi