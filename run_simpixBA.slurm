#!/bin/bash
#SBATCH --job-name=simpixBA
#SBATCH --time=4:00:00
#SBATCH -A phys56xx
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=simpixBA_%j.log

# Setup environment (same as interactive)
source ~/.bashrc
conda activate phys56xx

# Go to submit directory
cd $SLURM_SUBMIT_DIR

# ---- Inputs ----
SRC="imageB.png"
TGT="imageA.png"
TRIALS=20000000

# ---- Run Simpix B -> A ----
./simpix "$SRC" "$TGT" out_BA.png collageBA.png "$TRIALS"

echo "Done: imageB.png -> imageA.png"
