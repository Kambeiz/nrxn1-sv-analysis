#!/bin/bash
#SBATCH --job-name=nrxn1_analysis
#SBATCH --account=def-yourlab
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/nrxn1_%j.out
#SBATCH --error=logs/nrxn1_%j.err

# NRXN1 Structural Variant Analysis Pipeline
# Compute Canada SLURM submission script

module load python/3.10
module load samtools/1.15
module load bcftools/1.15

# Activate virtual environment
source ~/envs/nrxn1/bin/activate

# Set variables
INPUT_VCF=${1:-"data/input.vcf.gz"}
OUTPUT_DIR=${2:-"results"}
CONFIG=${3:-"config/default.yaml"}

echo "Starting NRXN1 analysis pipeline"
echo "Input: ${INPUT_VCF}"
echo "Output: ${OUTPUT_DIR}"
echo "Job ID: ${SLURM_JOB_ID}"

# Create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Run the pipeline
python main.py \
    --input ${INPUT_VCF} \
    --output ${OUTPUT_DIR} \
    --config ${CONFIG} \
    --log-level INFO \
    --log-file logs/nrxn1_${SLURM_JOB_ID}.log

echo "Pipeline completed"
