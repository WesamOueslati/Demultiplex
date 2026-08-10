#!/bin/bash
#SBATCH --account=bgmp   ### change this to your actual account for charging
#SBATCH --partition=bgmp     ### queue to submit to
#SBATCH --job-name=demult ### job name
#SBATCH --output=demult.out   ### file in which to store job stdout
#SBATCH --error=demult.err    ### file in which to store job stderr


DATA=/projects/bgmp/shared/2017_sequencing/

R1=${DATA}1294_S1_L008_R1_001.fastq.gz
R2=${DATA}1294_S1_L008_R2_001.fastq.gz
R3=${DATA}1294_S1_L008_R3_001.fastq.gz
R4=${DATA}1294_S1_L008_R4_001.fastq.gz


/usr/bin/time -v pixi run ./demultiplex.py -f $R1 -1 $R2 -2 $R3 -r $R4