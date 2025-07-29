#!/usr/bin/env bash
set -e

rm -f results/pre_ref_seqs_nuc.fasta
rm -f results/pre_ref_seqs_prot.fasta
rm -f results/queries_nuc.fasta
rm -f results/queries_prot.fasta

rm -rf results/er_0/nuc
rm -rf results/er_0/prot
rm -rf results/er_0/blastN/iota
rm -rf results/er_0/blastN/lambda
rm -rf results/er_0/blastP/iota
rm -rf results/er_0/blastP/lambda
rm -rf results/er_0/blastX/iota
rm -rf results/er_0/blastX/lambda
rm -rf results/er_0/tBlastN/iota
rm -rf results/er_0/tBlastN/lambda

rm -rf results/er_0.05/nuc
rm -rf results/er_0.05/prot
rm -rf results/er_0.05/blastN/iota
rm -rf results/er_0.05/blastN/lambda
rm -rf results/er_0.05/blastP/iota
rm -rf results/er_0.05/blastP/lambda
rm -rf results/er_0.05/blastX/iota
rm -rf results/er_0.05/blastX/lambda
rm -rf results/er_0.05/tBlastN/iota
rm -rf results/er_0.05/tBlastN/lambda

rm -rf results/er_0.025/nuc
rm -rf results/er_0.025/prot
rm -rf results/er_0.025/blastN/iota
rm -rf results/er_0.025/blastN/lambda
rm -rf results/er_0.025/blastP/iota
rm -rf results/er_0.025/blastP/lambda
rm -rf results/er_0.025/blastX/iota
rm -rf results/er_0.025/blastX/lambda
rm -rf results/er_0.025/tBlastN/iota
rm -rf results/er_0.025/tBlastN/lambda
