#Pipeline atualizado em 28/12/2024
#1. Ajustes das variaves de ambiente para facilitar na reexecução quando os diretorios forem atualizado
#2. Alterado o diretorio onde estava o PON na scratch45 para /home/users/vlira/PanelOfNormals/PoN.100COVID.vcf.gz


############# Definicao de variaveis para as amostras originais #################

export wd="/home/scratch90/bleocata_20241205/MMRF-COMMPASS/Variant_Calling"
export path="/home/scratch90/vlira_18nov2024"
export BAM_LIST="/home/scratch90/bleocata_20241205/MMRF-COMMPASS/Variant_Calling/Cinco_samples_BAMs"

export MEM=100

export pre_processing_DIR="${wd}/preprocessing_result"

export REF_FASTA="$wd/TESTE/GRCh38.d1.vd1.fa"
export GATK="$path/tools/gatk-4.6.0.0/./gatk"
export PICARD="${path}/tools/picard-3.2.0/picard.jar"

export TARGET="$path/references/xgen-exome-research-panel-v2-targets-hg38.bed"
export GNOMAD="$path/references/af-only-gnomad.hg38.vcf.gz" #frequência de variantes genéticas em populações diversas
#export GNOMAD="$wd/references/af-only-gnomad.SABE1171.Abraom.hg38.vcf.gz"
export ANNOVAR="$path/tools/annovar/table_annovar.pl"
export ANNOVAR_DB="$path/humandb"
export CROSS_REFERENCE="$path/references/refGene_TARGET_COSMICv82CensusGene_F1.txt"

export INDEL_KNOWN="/home/projects2/LIDO/molPathol/oncoseek/nextseq/references/Mills_and_1000G_gold_standard.indels.hg38.vcf"
export DBSNP="$path/references/dbsnp_151.hg38.vcf.gz"

########################################## Variaveis para Amostras originais mutect2 #############################################################3

export wd="/home/scratch90/bleocata_20241205/MMRF-COMMPASS/Variant_Calling/"
export path="/home/scratch90/vlira_18nov2024"
export BAM_LIST="/home/scratch90/bleocata_20241205/MMRF-COMMPASS/Variant_Calling/Cinco_samples_BAMs"

export MEM=60 #limite maximo por tarefa (considerar paralelo)

export pre_processing_DIR="${wd}/preprocessing_result"
export OUTPUT_DIR="${wd}/Result_Mutect2.GATK4.6"
export LOG_DIR="${wd}/Arquivos_log"
export LOG_FILE="${LOG_DIR}/Pipeline_Mutect2_$(date +%Y%m%d).log"

export REF_FASTA="/home/projects2/LIDO/molPathol/oncoseek/nextseq/hg38/Homo_sapiens_assembly38.fasta"
export GATK="$path/tools/gatk-4.6.0.0/./gatk"
export PICARD="java -jar ${path}/tools/picard-3.2.0/picard.jar"

export TARGET="/home/scratch90/bleocata_20241205/MMRF-COMMPASS/Variant_Calling/exons_basic_hg38v47_chr.bed"
export PON="/home/scratch90/bleocata_20241205/MMRF-COMMPASS/Variant_Calling/somatic-hg38_1000g_pon.hg38.vcf.gz" #precisa ser gz
export GNOMAD="$path/references/af-only-gnomad.hg38.vcf.gz" #frequência de variantes genéticas em populações diversas
#export GNOMAD="$wd/references/af-only-gnomad.SABE1171.Abraom.hg38.vcf.gz"
export ANNOVAR="$path/tools/annovar/table_annovar.pl"
export ANNOVAR_DB="$path/humandb/"
export CROSS_REFERENCE="$path/references/refGene_TARGET_COSMICv82CensusGene_F1.txt"

export INDEL_KNOWN="/home/projects2/LIDO/molPathol/oncoseek/nextseq/references/Mills_and_1000G_gold_standard.indels.hg38.vcf"
export DBSNP="$path/references/dbsnp_151.hg38.vcf.gz"


#######################################
########### PRE-PROCECING #############
#######################################
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery

#1_STEP. Raw Unmapped Reads 
#2_STEP. Map to Reference: STARalignment & HISAT2
#3_STEP. Raw Mapped Reads: SAMtools view; SAMtools sort e SAMtools index


echo "            >>>>>> Starting Pipeline for Pre-processing Variant Calling  <<<<<< $(date) " >> "$pre_processing_DIR/preprocessing_$(date +%Y%m%d).log"
mkdir -p "$wd"
mkdir -p "$pre_processing_DIR"

######## 5_STEP. MarkDuplicates : Picard tool ############

stage_MarkDuplicates() {
  local bam_file="$1"
  local id
  id=$(basename "$bam_file" .rna_seq.genomic.gdc_realn.bam)

  # Mensagem de início no log
  echo ">>> Iniciando processamento com MarkDuplicates: $(date) <<<" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  echo ">>>>>> Processando arquivo: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

  # Execução do comando MarkDuplicates
  java -jar "$PICARD" MarkDuplicates \
      -I "$wd/Cinco_samples_BAMs/${id}.rna_seq.genomic.gdc_realn.bam" \
      -O "$pre_processing_DIR/marked_duplicates/${id}_marked_duplicates.bam" \
      -M "$pre_processing_DIR/marked_duplicates/${id}_marked_dup_metrics.txt" \
      2> "$pre_processing_DIR/marked_duplicates/${id}_marked_duplicates.bam.log"
      
      
if [[ $? -eq 0 ]]; then
        echo ">>>>>> MarkDuplicates concluído com sucesso para: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
    else
        echo ">>>>>> Erro ao executar MarkDuplicates para: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
    fi
}
export -f stage_MarkDuplicates

mkdir -p "$pre_processing_DIR/marked_duplicates"
find "$BAM_LIST" -name "*.bam" | parallel -j 4 stage_MarkDuplicates {}



######### 6_STEP. SplitNCigarRead: GATK Tools ############
#Adicionar ou recalcular os campos NM (número de diferenças), MD (mismatch descriptor), e UQ (qualidade única de alinhamento) nos alinhamentos de um arquivo BAM usando a ferramenta SetNmMdAndUqTags do Picard Tools.

echo ">>>>>> SplitNCigarRead: GATK Tools <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

SplitNCigarRead() {
  local sample="$1"
  local id
  id=$(basename "$sample" .rna_seq.genomic.gdc_realn.bam)

  echo ">>>>>> Iniciando SplitNCigarRead para Amostra: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

  $GATK --java-options "-Xmx${MEM}G" SplitNCigarReads \
      -R "${REF_FASTA}" \
      -I "${pre_processing_DIR}/marked_duplicates/${id}_marked_duplicates.bam" \
      -O "$pre_processing_DIR/SplitNCigarRead/${id}_marked_duplicates.tags.bam" \
      --tmp-dir "${pre_processing_DIR}/tmp_GATK" \
      2> "$pre_processing_DIR/SplitNCigarRead/${id}_marked_duplicates.tags.bam.log"
      
    if [[ $? -eq 0 ]]; then
        echo ">>>>>> SplitNCigarReads concluído com sucesso para: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
    else
        echo ">>>>>> Erro ao executar SplitNCigarReads para: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
    fi
}

export -f SplitNCigarRead

mkdir -p "$pre_processing_DIR/SplitNCigarRead"
mkdir -p "$pre_processing_DIR/tmp_GATK"
find "$BAM_LIST" -name "*.bam" | parallel -j 2 SplitNCigarRead {}


######### 7_STEP. Recalibrate Base Quality Score(BQSR): BaseRecalibrator##########

stage_base_recalibration(){
  local sample="$1"
  local id
  id=$(basename "$sample" .rna_seq.genomic.gdc_realn.bam)
  echo ">>>>>> Executando stage_base_recalibration para Amostra: $id <<<  $(date) " >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

  #Gera uma tabela de recalibracao com base nos padroes de erros detectados
  $GATK --java-options "-Xmx${MEM}G" BaseRecalibrator \
        -R "$REF_FASTA"\
        -I "${pre_processing_DIR}/SplitNCigarRead/${id}_marked_duplicates.tags.bam" \
        -known-sites ${INDEL_KNOWN} \
        -known-sites ${DBSNP} \
        -known-sites ${GNOMAD} \
        -L ${TARGET} \
        -O ${pre_processing_DIR}/BQSR_base_recalibration/"${id}_marked_duplicates.tags.recal.table" \
        > ${pre_processing_DIR}/BQSR_base_recalibration/"${id}_marked_duplicates.tags.recal.table.log" \
        2> ${pre_processing_DIR}/BQSR_base_recalibration/"${id}_marked_duplicates.tags.table.log2"
        
   if [[ $? -eq 0 ]]; then
      echo ">>>>>> BaseRecalibrator concluído com sucesso para: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  else
      echo ">>>>>> Erro ao executar BaseRecalibrator para: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
      return 1
  fi
  
  
  #Aplica as correcoes calculadas pelo BaseRecalibrator (table) ao arquivo .bam
  $GATK --java-options "-Xmx${MEM}G" ApplyBQSR \
        -R $REF_FASTA \
        -I ${pre_processing_DIR}/SplitNCigarRead/"${id}_marked_duplicates.tags.bam" \
        --bqsr-recal-file ${pre_processing_DIR}/BQSR_base_recalibration/"${id}_marked_duplicates.tags.recal.table" \
        -O ${pre_processing_DIR}/BQSR_base_recalibration/"${id}_marked_duplicates.tags.recal.bam" \
        > ${pre_processing_DIR}/BQSR_base_recalibration/"${id}_marked_duplicates.tags.recal.bam.log" \
        2> ${pre_processing_DIR}/BQSR_base_recalibration/"${id}_marked_duplicates.tags.recal.bam.log2"
        
  if [[ $? -eq 0 ]]; then
      echo ">>>>>> ApplyBQSR concluído com sucesso para: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  else
      echo ">>>>>> Erro ao executar ApplyBQSR para: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
      return 1
  fi
        
  #gerar um grafico para avaliacao da recalibracao 
  #$GATK --java-options "-Xmx${MEM}G" AnalyzeCovariates \
       #-bqsr ${pre_processing_DIR}/BQSR_base_recalibratio/"${id}_marked_duplicates.tags.recal.table"
       #-plots ${pre_processing_DIR}/BQSR_base_recalibratio/AnalyzeCovariates.pdf
}
export -f  stage_base_recalibration

mkdir -p "$pre_processing_DIR/BQSR_base_recalibration"
find "$BAM_LIST" -name "*.bam"  | parallel -j 4 stage_base_recalibration {}

