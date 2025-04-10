#!/usr/bin/env bash
set -e

IN_QUERY_FILE=$1
IN_INDEX_FILE=$2
OUT_FILE=$3
LOG_FILE=$4
PATTERN=$5
ERROR=$6
THREADS=$7

valik search --index $IN_INDEX_FILE \
             --query $IN_QUERY_FILE \
             --threads $THREADS \
             --pattern $PATTERN \
             --error $ERROR \
             --output $OUT_FILE \
             &> $LOG_FILE