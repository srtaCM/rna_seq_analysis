#!/bin/bash
#SBATCH --job-name=bowtie_GFP10_11
#SBATCH --account=rnaseq-pilotlab 
#SBATCH --output=bowtie_GFP10_11.%j.out
#SBATCH --error=bowtie_GFP10_11.%j.err
#SBATCH --mail-user=iliana@vt.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --partition=normal_q
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -e  # Exit on error

########################################
# Load Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

echo "Conda environment:"
conda list | grep -E "samtools|bowtie2"

########################################
# Directories

PROJDIR=/projects/rnaseq-pilotlab/icm/TRAP_CLE
INDATADIR=$PROJDIR/reports/fastp/TRAP+CLE_trimmed
OUTDIR=$PROJDIR/H_ratio/bowtie/GFP10_11

mkdir -p "$OUTDIR"

########################################
# Build Bowtie2 pseudo-reference index
########################################
cd $OUTDIR

if [ ! -f pseudo_ref_index.1.bt2 ]; then
  echo "Building bowtie2 index..."
  bowtie2-build --threads ${SLURM_CPUS_PER_TASK} \
      $PROJDIR/H_ratio/GFP10_11.fa \
      pseudo_ref_index
else
  echo "Index already exists. Skipping build."
fi

cd $INDATADIR

########################################
# Samples
########################################
samples=(
DIC1_S21_L002
DIC2_S22_L002
DIC3_S23_L002
DIC4_S24_L002
DIT1_S5_L001
DIT2_S6_L001
DIT3_S7_L001
DIT4_S8_L001
DIT5_S9_L001
DMC1_S18_L002
DMC2_S19_L002
DMC3_S20_L002
PIC1_S14_L002
PIC2_S15_L002
PIC3_S16_L002
PIC4_S17_L002
PIT1_S3_L001
PIT2_S1_L001
PIT3_S2_L001
PIT4_S3_L001
PIT5_S4_L001
PMC1_S10_L001
PMC2_S11_L001
PMC3_S12_L001
PMC4_S13_L002
PMT1_S10_L002
PMT2_S11_L002
PMT3_S12_L002
PMT4_S13_L002
PMT5_S14_L002
)

########################################
# Mapping
########################################
for sample in "${samples[@]}"; do

    FWD1=$INDATADIR/${sample}_R1.trimmed.fq
    REV2=$INDATADIR/${sample}_R2.trimmed.fq

    # Check if input files exist
    if [ ! -f "$FWD1" ] || [ ! -f "$REV2" ]; then
      echo "WARNING: Input files not found for $sample. Skipping."
      continue
    fi

    echo "Processing $sample at $(date)"

    bowtie2 \
        -x $OUTDIR/pseudo_ref_index \
        -1 $FWD1 \
        -2 $REV2 \
        --very-sensitive-local \
        -k 10 \
        -p ${SLURM_CPUS_PER_TASK} \
        -S $OUTDIR/${sample}.aligned.sam \
        2> $OUTDIR/${sample}.mapping_summary.log

    samtools view -bS $OUTDIR/${sample}.aligned.sam > $OUTDIR/${sample}.aligned.bam
    samtools sort $OUTDIR/${sample}.aligned.bam -o $OUTDIR/${sample}.aligned.sorted.bam
    samtools index $OUTDIR/${sample}.aligned.sorted.bam
    
    # Extract mapping statistics
    samtools flagstat $OUTDIR/${sample}.aligned.sorted.bam > $OUTDIR/${sample}.flagstat.txt
    
    # Count reads mapped to GFP10_11 specifically
    samtools view -c -F 4 $OUTDIR/${sample}.aligned.sorted.bam > $OUTDIR/${sample}.GFP10_11_mapped_count.txt

    # Clean up intermediate SAM
    rm $OUTDIR/${sample}.aligned.sam
    
    echo "$sample completed at $(date)"
done

########################################
# Generate mapping summary table
########################################
echo "All samples mapped. Creating summary..."
echo -e "sample_id\tGFP10_11_mapped_reads\ttotal_reads\tmapped_pct" > $OUTDIR/GFP10_11_mapping_summary.tsv

for sample in "${samples[@]}"; do
    count_file=$OUTDIR/${sample}.GFP10_11_mapped_count.txt
    
    if [ ! -f "$count_file" ]; then
      echo "WARNING: Count file not found for $sample"
      continue
    fi
    
    mapped=$(cat $count_file)
    total=$(samtools view -c $OUTDIR/${sample}.aligned.sorted.bam)
    
    if [ "$total" -gt 0 ]; then
      pct=$((mapped * 100 / total))
    else
      pct=0
    fi
    
    echo -e "${sample}\t${mapped}\t${total}\t${pct}" >> $OUTDIR/GFP10_11_mapping_summary.tsv
done

echo "Job finished at $(date)"

