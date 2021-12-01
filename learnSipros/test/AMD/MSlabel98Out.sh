#!/bin/bash
#SBATCH --partition=omicsbio
#SBATCH --ntasks=40
#SBATCH --time=60:00:00
#SBATCH --job-name=plink
#SBATCH --mem=150G
#SBATCH --chdir=/scratch/yixiong/learnSipros/test/AMD
#SBATCH --output=MSlabel98Out_stdout.txt
#SBATCH --error=MSlabel98Out_stderr.txt
#SBATCH --mail-user=yixiong@ou.edu
#SBATCH --mail-type=ALL

cd /scratch/yixiong/learnSipros/test/AMD
../Sipros \
  -g configFiles \
  -w MSlabel98 \
  -o MSlabel98Out
