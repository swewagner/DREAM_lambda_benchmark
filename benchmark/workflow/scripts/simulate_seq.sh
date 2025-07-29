#!/usr/bin/env bash
set -e

OUT_FILE=$1
NUM=$2
LEN=$3
LOG_FILE=$4
BINARY_DIR="../build/mason2/src/mason2-build/bin"

truncate -s 0 $OUT_FILE

for i in $(seq 1 $NUM )
do
    # random seed generation
    seed=$(python3 <<EOF
import random
random.seed(None)
print(random.randint(0, int(1e6)))
EOF
    )
    # create a sequence of length ref_len
    $BINARY_DIR/mason_genome -s $seed -l $LEN -o $OUT_FILE.temp.fasta &> $LOG_FILE
    # convert multi line fasta to one line fasta
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $OUT_FILE.temp.fasta > $OUT_FILE.temp2.fasta
    sed -i '1,2d' $OUT_FILE.temp2.fasta
    # write sequence to output file
    echo ">$i" >> $OUT_FILE
    cat $OUT_FILE.temp2.fasta >> $OUT_FILE
    rm $OUT_FILE.temp.fasta $OUT_FILE.temp2.fasta
done