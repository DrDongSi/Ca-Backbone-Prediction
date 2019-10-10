#!/bin/bash
#--------------------------------------------------------------------------------
#  SBATCH CONFIG
#--------------------------------------------------------------------------------
#SBATCH -J  3n2t
#SBATCH -o 3n2t-%j.out
#SBATCH --partition Lewis,hpc5,hpc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --time 2-00:00
#--------------------------------------------------------------------------------

outputdir=/home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/test/3n2t_out

mkdir -p /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/test/3n2t_out

cd /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/test/3n2t_out

printf "perl /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/scripts/CaTrace2Seq.pl /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/examples/3n2t/3n2t_fragment.pdb /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/examples/3n2t/3n2t.fasta /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/test/3n2t_out 50 10\n\n"

perl /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/scripts/CaTrace2Seq.pl /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/examples/3n2t/3n2t_fragment.pdb /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/examples/3n2t/3n2t.fasta /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/test/3n2t_out 50 10

