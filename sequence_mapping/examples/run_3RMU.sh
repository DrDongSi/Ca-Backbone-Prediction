#!/bin/bash
#--------------------------------------------------------------------------------
#  SBATCH CONFIG
#--------------------------------------------------------------------------------
#SBATCH -J  3RMU
#SBATCH -o 3RMU-%j.out
#SBATCH --partition Lewis,hpc5,hpc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --time 2-00:00
#--------------------------------------------------------------------------------

outputdir=/home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/test/3RMU_out

mkdir -p /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/test/3RMU_out

cd /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/test/3RMU_out

printf "perl /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/scripts/CaTrace2Seq.pl /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/examples/3RMU/3RMU_fragment.pdb /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/examples/3RMU/3RMU.fasta /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/test/3RMU_out 50 10\n\n"

perl /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/scripts/CaTrace2Seq.pl /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/examples/3RMU/3RMU_fragment.pdb /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/examples/3RMU/3RMU.fasta /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/test/3RMU_out 50 10

