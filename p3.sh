#!/bin/bash
#SBATCH -A hugen2072-2026s
#SBATCH -N 1
#SBATCH -c 24
#SBATCH --output=Log_%x.out
#SBATCH --time=48:00:00
#SBATCH --mail-user=tsb31@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL,TIMEOUT

# load all modules
module load fastqc/0.11.9
module load cutadapt/5.1
module load gcc/8.2.0
module load bwa/0.7.17
module load samtools/1.22.1
module load gatk/4.5.0.0

set -euo pipefail
set -v

# initial setup of variables
scratch=$SLURM_SCRATCH
outdir=$SLURM_SUBMIT_DIR/p3
mkdir -p $outdir
p3=/ix1/hugen2072-2026s/p3

# test
# fastq1=$p3/toy_1.fastq.gz
# fastq2=$p3/toy_2.fastq.gz
# id=toy
# reference=$p3/human_g1k_v37.fasta

# indels=$p3/Mills_and_1000G_gold_standard.indels.b37.vcf
# snps=$p3/dbsnp_138.b37.vcf
# hapmap=$p3/hapmap_3.3.b37.vcf
# omni=$p3/1000G_omni2.5.b37.vcf
# phase=$p3/1000G_phase1.snps.high_confidence.b37.vcf

# prod
fastq1=$p3/p3_1.fastq.gz
fastq2=$p3/p3_2.fastq.gz
id=prod
reference=$p3/Homo_sapiens_assembly38.fasta

indels=$p3/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
snps=$p3/dbsnp_146.hg38.vcf.gz
hapmap=$p3/hapmap_3.3.hg38.vcf.gz
omni=$p3/1000G_omni2.5.hg38.vcf.gz
phase=$p3/1000G_phase1.snps.high_confidence.hg38.vcf.gz

outdir=$outdir/$id
mkdir -p $outdir
scratch=$scratch/$id
mkdir -p $scratch

# variables for checkpoint cram files
cram1cp=$outdir/cram1.cram
cram1=$scratch/cram1.cram
cram2cp=$outdir/cram2.cram
cram2=$scratch/cram2.cram

if [ -f "$cram1cp" ]; then
    echo "$cram1cp already exists ... jumping ahead."
    # move the cram file to scratch
    cp $cram1cp $cram1
else
    echo "Generating $cram1 ..."
    mkdir -p $outdir/pre
    mkdir -p $outdir/post
    # generating fastqc files before cutadapt
    fastqc $fastq1 -t 24 --outdir=$outdir/pre
    fastqc $fastq2 -t 24 --outdir=$outdir/pre
    # running cutadapt
    cutadapt -j 0 -m 10 -q 20 $fastq1 $fastq2 \
        -a AGATCGGAAGAG -A AGATCGGAAGAG \
        -o $scratch/1_trimmed.fastq.gz -p $scratch/2_trimmed.fastq.gz
    # generating fastqc files after cutadapt
    fastqc $scratch/1_trimmed.fastq.gz -t 24 --outdir=$outdir/post
    fastqc $scratch/2_trimmed.fastq.gz -t 24 --outdir=$outdir/post
    # run the alignment
    bwa mem -t 24 \
        $reference \
        $scratch/1_trimmed.fastq.gz $scratch/2_trimmed.fastq.gz \
        -R "@RG\tID:p3\tLB:P5\tSM:P5\tPL:ILLUMINA" \
        > $scratch/aligned.sam
    # convert to bam
    samtools view -b $scratch/aligned.sam | samtools sort -o $scratch/aligned.bam
    # index the bam file
    samtools index $scratch/aligned.bam
    # convert to a cram file
    samtools view -C -T $reference \
        --output-fmt-option version=3.0 \
        -o $cram1 $scratch/aligned.bam
    # generate the checkpoint
    cp $cram1 $cram1cp
    # clear the sam and bam files
    rm $scratch/aligned.sam $scratch/aligned.bam
fi
# either way, reindex the cram file here as we keep that in scratch
samtools index $cram1

if [ -f "$cram2cp" ]; then
    echo "$cram2cp already exists ... jumping ahead."
    # move the cram file to scratch
    cp $cram2cp $cram2
else
    echo "Generating $cram2 ..."
    # mark duplicates
    gatk MarkDuplicatesSpark -I $cram1 \
        -O $scratch/dupsmarked.bam \
	-R $reference
    # base recalibration calculation
    gatk BaseRecalibrator --java-options "-XX:ParallelGCThreads=24" \
        -I $scratch/dupsmarked.bam \
        -R $reference \
        -O $scratch/BQSR.table \
        --known-sites $indels \
        --known-sites $snps
    # apply the base recalibration calculations
    gatk ApplyBQSR --java-options "-XX:ParallelGCThreads=24" \
        -R $reference \
        -I $scratch/dupsmarked.bam \
        --bqsr-recal-file $scratch/BQSR.table \
        -O $scratch/dupsmarked_cleaned.bam
        # convert to a cram file
        samtools view -C -T $reference \
            --output-fmt-option version=3.0 \
            -o $cram2 $scratch/dupsmarked_cleaned.bam
        # generate the checkpoint
        cp $cram2 $cram2cp
        # clear the bam files
        rm $scratch/dupsmarked.bam $scratch/dupsmarked_cleaned.bam
        # keep the old cram file to enforce checkpoint, but clear the file of contents to save space
        echo "" > $cram1
fi
# either way, reindex the cram file here as we keep that in scratch
samtools index $cram2

# generate alignment statistics
samtools flagstat -@ 24 $cram2 \
    > $outdir/alignment_statistics.out
samtools depth $cram2 \
    | gzip > $outdir/depth_statistics.out.gz

# genotyping
gatk HaplotypeCaller -R $reference \
    -I $cram2 \
    -O $scratch/genotypes.g.vcf.gz \
    -ERC GVCF \
    -OVI \
    --native-pair-hmm-threads 24

# call genotypes
gatk GenotypeGVCFs --java-options "-XX:ParallelGCThreads=24" \
    -R $reference \
    -V $scratch/genotypes.g.vcf.gz \
    -O $scratch/genotypes.vcf.gz

# make the sites only file
gatk MakeSitesOnlyVcf -I $scratch/genotypes.vcf.gz -O $scratch/sites_only.vcf.gz

# calculate indel recalibration
gatk VariantRecalibrator --java-options "-XX:ParallelGCThreads=24" \
    -mode INDEL \
    -R $reference \
    -V $scratch/sites_only.vcf.gz \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
    -resource:mills,known=false,training=true,truth=true,prior=12 \
        $indels \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2 \
        $snps \
    -O $scratch/indels.recal \
    --tranches-file $scratch/indels.tranches

# calculate snp recalibration
gatk VariantRecalibrator --java-options "-XX:ParallelGCThreads=24" \
    -mode SNP \
    -R $reference \
    -V $scratch/sites_only.vcf.gz \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    -resource:hapmap,known=false,training=true,truth=true,prior=15 \
        $hapmap \
    -resource:omni,known=false,training=true,truth=true,prior=12 \
        $omni \
    -resource:1000G,known=false,training=true,truth=false,prior=10 \
        $phase \
    -resource:dbsnp,known=true,training=false,truth=false,prior=7 \
        $snps \
    -O $scratch/snps.recal \
    --tranches-file $scratch/snps.tranches

# apply indel recalibration
gatk ApplyVQSR --java-options "-XX:ParallelGCThreads=24" \
    -mode INDEL \
    -R $reference \
    -V $scratch/genotypes.vcf.gz \
    --recal-file $scratch/indels.recal \
    --tranches-file $scratch/indels.tranches \
    --truth-sensitivity-filter-level 99.0 \
    --create-output-variant-index true \
    -O $scratch/genotypes_indelqc.vcf.gz

# apply snp recalibration
gatk ApplyVQSR --java-options "-XX:ParallelGCThreads=24" \
    -mode SNP \
    -R $reference \
    -V $scratch/genotypes_indelqc.vcf.gz \
    --recal-file $scratch/snps.recal \
    --tranches-file $scratch/snps.tranches \
    --truth-sensitivity-filter-level 99.0 \
    --create-output-variant-index true \
    -O $outdir/final.vcf.gz

# calculate final statistics
gatk CollectVariantCallingMetrics -I $outdir/final.vcf.gz \
    --DBSNP $snps \
    -O $outdir/genotype_metrics
