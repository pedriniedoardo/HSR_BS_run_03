#!/bin/bash
#SBATCH --job-name=DE_pbulk_BS
#SBATCH --account pedrini.edoardo
#SBATCH --mem=128GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=6  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="06_DE_pbulk_fullObject.err"
#SBATCH --output="06_DE_pbulk_fullObject.out"

PROJECT_ROOT="/beegfs/scratch/ric.cosr/pedrini.edoardo/scRNAseq_Brainsphere_Absinta/BS_batch_2023/analysis/R_analysis_R45_pixi"

echo "my job strart now" > 06_DE_pbulk_fullObject.log;

date >> 06_DE_pbulk_fullObject.log;

# pixi selects R first, then renv (activated from .Rprofile) manages the package library:
# cd into PROJECT_ROOT (where .Rprofile/renv.lock/pixi.toml live) before `pixi run Rscript`.
# the relative paths in the R script (../../out/...) are resolved from PROJECT_ROOT.
cd "$PROJECT_ROOT" || exit 1

pixi run Rscript "$PROJECT_ROOT/scr/06_DE_pbulk_fullObject.R"

date >> 06_DE_pbulk_fullObject.log;
echo "all done!!" >> 06_DE_pbulk_fullObject.log
