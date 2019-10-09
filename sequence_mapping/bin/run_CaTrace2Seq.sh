#!/bin/bash
#--------------------------------------------------------------------------------
#  SBATCH CONFIG
#--------------------------------------------------------------------------------
#SBATCH -J  CaTrace2Seq
#SBATCH -o CaTrace2Seq-%j.out
#SBATCH --partition Lewis,hpc5,hpc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --time 2-00:00
#--------------------------------------------------------------------------------



if [ $# != 5 ]; then
	echo "$0 <path of fasta sequence> <output-directory> <length threshold for fragment> <number of cpus>"
	exit
fi

fasta_file=$1
Ca_trace_file=$2
outputdir=$3
threshold=$4
cpu_num=$5


printf "perl /scripts/CaTrace2Seq.pl $Ca_trace_file $fasta_file $outputdir $threshold $cpu_num\n\n"

perl /home/jh7x3/Ca-Backbone-Prediction/sequence_mapping/scripts/CaTrace2Seq.pl $Ca_trace_file $fasta_file $outputdir $threshold $cpu_num





  

