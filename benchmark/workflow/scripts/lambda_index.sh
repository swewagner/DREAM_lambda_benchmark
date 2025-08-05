#!/usr/bin/env bash
set -e

IN_FILE=$1
OUT_FILE=$2
BLAST_M=$3
LOG_FILE=$4
BINARY_DIR="../build/lambda3/src/lambda3-build/bin"

if [ $BLAST_M = "blastN" ]; then
$BINARY_DIR/lambda3 mkindexn -d $IN_FILE -i $OUT_FILE -v 2 &> $LOG_FILE
else
$BINARY_DIR/lambda3 mkindexp -d $IN_FILE -i $OUT_FILE -v 2 &> $LOG_FILE
fi