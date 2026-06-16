#!/bin/bash
#SBATCH -A hugen2072-2026s
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J align_rna_seq
#SBATCH --output=align_rna_seq.out
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process.
#SBATCH --mem=64g # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --mail-user=tsb31@pitt.edu
#SBATCH --mail-type=ALL

##array should start from zero

# Load modules
module load cutadapt/5.1
module load fastqc/0.12.1
module load star/2.7.11b
module load gcc/15.1.0

DATA="/ix1/hugen2072-2026s/p5"
OUT="$HOME/p5"
F1="HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read"
F2="HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read"
F3="UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read"
F4="UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read"
inputs=($F1 $F2 $F3 $F4)

set -euo pipefail
set -v

# Aligns RNA-seq data given a file path and outputs a file name with a given prefix.
# The first argument is the input file name prefix across both reads.
# The second argument is the directory where the files live.
# The third argument is the output directory.
# The prefix is pulled from the first 8 characters of the input file prefix (argument 1).
function align_rna_seq() {
  # Use cutadapt to trim adapter sequence
  # Right now the script is set up to take in HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz
  # 	and HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
  #   and it is set to output in HBR_Rep1.trimmed.R1.fastq.gz and HBR_Rep1.trimmed.R2.fastq.gz
  # You need to change the input for each run (directory stays the same)
  # You need to change the output file name for each run (use the output file names given in the assignment instructions)
  # (You need to change the output file directory from blah/blah/blah/project5/ to a folder of yours; make sure the folder exists)
  file=$1
  data=$2
  out=$3
  in1="$data/${file}1.fastq.gz"
  in2="$data/${file}2.fastq.gz"
  prefix="$(echo $file | awk '{print substr($0, 1, 8)}')"
  trimmed1="$out/${prefix}.trimmed.R1.fastq.gz"
  trimmed2="$out/${prefix}.trimmed.R2.fastq.gz"
  cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -o $trimmed1 \
  -p $trimmed2 \
  $in1 $in2
  
  
  # Aligning the reads
  # The inputs are the two files you made in the previous step
  # (so remember to update them  each time you run the script)
  # (and change the folder from /blah/blah/blah/project5/ of course)
  # For the --outFileNamePrefix argument, name the file HBR_Rep1, HBR_Rep2, UHR_Rep1, or UHR_Rep2 depending on which step you're on
  # 	also include the path to the folder you're using (i.e., replace /blah/blah/blah/project5)
  # Don't change the --genomeDir argument at all
  STAR \
  --runMode alignReads \
  --outSAMstrandField intronMotif \
  --twopassMode Basic \
  --readFilesIn $trimmed1 \
                $trimmed2 \
  --outFileNamePrefix "$out/$prefix" \
  --quantMode GeneCounts \
  --outStd Log \
  --outWigType bedGraph \
  --outWigReferencesPrefix ./bedgraph \
  --genomeDir /ix1/hugen2072-2026s/p5/STAR_reference/ \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate
}

for input in "${inputs[@]}"; do
  echo "---------------------------------------------"
  echo "Aligning: '$input'"
  align_rna_seq $input $DATA $OUT
  echo "---------------------------------------------"
done
