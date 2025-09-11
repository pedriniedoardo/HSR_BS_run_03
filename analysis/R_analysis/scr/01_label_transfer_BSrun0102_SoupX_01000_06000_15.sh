#!/bin/bash
#SBATCH --job-name=scProcess
#SBATCH --account pedrini.edoardo
#SBATCH --mem=64GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=6  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="01_label_transfer_BSrun0102_SoupX_01000_06000_15.err"
#SBATCH --output="01_label_transfer_BSrun0102_SoupX_01000_06000_15.out"

echo "my job strart now" > 01_label_transfer_BSrun0102_SoupX_01000_06000_15.log;

date >> 01_label_transfer_BSrun0102_SoupX_01000_06000_15.log;

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate env_R4.2_renv;

Rscript /beegfs/scratch/ric.cosr/pedrini.edoardo/scRNAseq_Brainsphere_Absinta/BS_batch_2023/analysis/R_analysis/scr/01_label_transfer_BSrun0102_SoupX_01000_06000_15.R

date >> 01_label_transfer_BSrun0102_SoupX_01000_06000_15.log;
echo "all done!!" >> 01_label_transfer_BSrun0102_SoupX_01000_06000_15.log