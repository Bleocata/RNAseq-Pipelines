export WD="/home/scratch90/bleocata_20240912/MMRF-COMMPASS"
export reference="home/scratch90/bleocata_20240912/PRJEB37100/STAR_alignment/Reference"
export JOBS=4
export LOG_FILE="${WD}/Pipeline_RNAseq_$(date +%Y%m%d).log"

local bam_file="$1"
    local id
    id=$(basename "$bam_file" .rna_seq.genomic.gdc_realn.bam)
    
rnaseqc \
  "${reference}"/gencode.v44.primary_assembly.genes.gtf \
  "${WD}"/samples_BAM_BAI/${id}.rna_seq.genomic.gdc_realn.bam \
  ${WD}/output_rnaseqc/${id}. \
  --bed test_data/downsampled.bed \
  --fasta
  --coverage
  