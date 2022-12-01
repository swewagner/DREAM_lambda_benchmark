#!/usr/bin/env bash
set -e

OUT_FILE=$1
NUM=$2
LEN=$3
BINARY_DIR="/home/darklyght/Documents/Studium/Bioinformatik/Master/Praddidum/Benchmarking/DREAM-stellar-benchmark/lib/raptor_data_simulation/build/src/mason2/src/mason2-build/bin"

truncate -s 0 $OUT_FILE

$BINARY_DIR/mason_genome $(eval echo | awk "{ for (counter=1; counter <= $NUM; counter++) printf \"-l $LEN \"; }") -o $OUT_FILE &> /dev/null
