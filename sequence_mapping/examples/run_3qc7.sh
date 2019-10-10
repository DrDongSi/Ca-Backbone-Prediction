#!/bin/bash
#--------------------------------------------------------------------------------
#  SBATCH CONFIG
#--------------------------------------------------------------------------------
#SBATCH -J  3qc7
#SBATCH -o 3qc7-%j.out
#SBATCH --partition Lewis,hpc5,hpc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --time 2-00:00
#--------------------------------------------------------------------------------

outputdir=/home/jh7x3/CaTrace2Seq/test/3qc7_out

mkdir -p /home/jh7x3/CaTrace2Seq/test/3qc7_out

cd /home/jh7x3/CaTrace2Seq/test/3qc7_out

printf "perl /home/jh7x3/CaTrace2Seq/scripts/CaTrace2Seq.pl /home/jh7x3/CaTrace2Seq/examples/3qc7/3qc7_fragment.pdb /home/jh7x3/CaTrace2Seq/examples/3qc7/3qc7.fasta /home/jh7x3/CaTrace2Seq/test/3qc7_out 50 10\n\n"

perl /home/jh7x3/CaTrace2Seq/scripts/CaTrace2Seq.pl /home/jh7x3/CaTrace2Seq/examples/3qc7/3qc7_fragment.pdb /home/jh7x3/CaTrace2Seq/examples/3qc7/3qc7.fasta /home/jh7x3/CaTrace2Seq/test/3qc7_out 50 10

