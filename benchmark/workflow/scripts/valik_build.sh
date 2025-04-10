#!/usr/bin/env bash
set -e

IN_FILE=$1
OUT_FILE=$2
LOG_FILE=$3
WINDOW=$4
KMER=$5
SIZE=$6

valik build $IN_FILE \
            --window $WINDOW \
            --kmer $KMER \
            --size $SIZE \
            --output $OUT_FILE \
             &> $LOG_FILE