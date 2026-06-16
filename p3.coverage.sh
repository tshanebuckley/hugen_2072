#!/bin/bash
#SBATCH -A hugen2072-2026s
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --output=Log_%x.out
#SBATCH --time=8:00:00
#SBATCH --mail-user=tsb31@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL,TIMEOUT

# decided to calculate read depth by getting the 'meandepth' field from following these docs:
# https://www.htslib.org/doc/samtools-coverage.html

set -euo pipefail
set -v

# load all modules
module load gcc/8.2.0
module load bwa/0.7.17
module load samtools/1.22.1

# initial setup of variables
scratch=$SLURM_SCRATCH
outdir=$SLURM_SUBMIT_DIR/p3/prod
mkdir -p $outdir
p3=/ix1/hugen2072-2026s/p3
cram1h=$outdir/cram1.cram
cram2h=$outdir/cram2.cram
cram1=$scratch/cram1.cram
cram2=$scratch/cram2.cram

# move the cram files to scratch and index them
cp $cram1h $cram1
cp $cram2h $cram2

# index the files now in scratch
samtools index $cram1
samtools index $cram2

# calculate the qc-passed read depth
samtools coverage $cram1 -o $outdir/qc.coverage

# calculate the mapped read depth
samtools coverage $cram2 -o $outdir/mapped.coverage
