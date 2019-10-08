#!/bin/bash
#--------------------------------------------------------------------------------
#  SBATCH CONFIG
#--------------------------------------------------------------------------------
#SBATCH -J  8410
#SBATCH -o 8410-%j.out
#SBATCH --partition Lewis,hpc5,hpc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --time 2-00:00
#--------------------------------------------------------------------------------

outputdir=/home/jh7x3/CaTrace2Seq/test/8410_out

mkdir -p /home/jh7x3/CaTrace2Seq/test/8410_out

cd /home/jh7x3/CaTrace2Seq/test/8410_out

printf "perl /home/jh7x3/CaTrace2Seq/scripts/CaTrace2Seq.pl /home/jh7x3/CaTrace2Seq/examples/8410/8410_fragment.pdb /home/jh7x3/CaTrace2Seq/examples/8410/8410.fasta /home/jh7x3/CaTrace2Seq/test/8410_out 50 10\n\n"

perl /home/jh7x3/CaTrace2Seq/scripts/CaTrace2Seq.pl /home/jh7x3/CaTrace2Seq/examples/8410/8410_fragment.pdb /home/jh7x3/CaTrace2Seq/examples/8410/8410.fasta /home/jh7x3/CaTrace2Seq/test/8410_out 50 10

