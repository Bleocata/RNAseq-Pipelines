#Pipeline atualizado em 14/06/2024
#1. Ajustes das variaves de ambiente para facilitar na reexecução quando os diretorios forem atualizado
#2. Alterado o diretorio onde estava o PON na scratch45 para /home/users/vlira/PanelOfNormals/PoN.100COVID.vcf.gz
#3. Agora o script recebe 2 agumentos: 1- lista de amostras; 2- diretorio SCRATCH
#4. Utualizado os databases: hg38_cosmic98_coding,hg38_avsnp150,hg38_clinvar_20220320, hg38_gnomad40_exome

#Melhorias pipeline: utilizar arquivos de eferencia sempre na mesa versao 
#Adcionar mais um chamador na pipeline
#utilizar as amostras PON de DNA d proprio projeto 
#realinhamento de resgate 


############# Definicao de variaveis #################

export wd="/home/scratch90/bleocata_20241205/MMRF-COMMPASS/Variant_Calling"
export path="/home/scratch90/vlira_18nov2024"
export BAM_LIST="/home/scratch90/bleocata_20241205/MMRF-COMMPASS/Variant_Calling/Cinco_samples_BAMs"

export MEM=200 #limite maximo por tarefa (considerar paralelo)

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
#Para processamento RNAseq, considerar: https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels

#1_STEP. Raw Unmapped Reads 
#2_STEP. Map to Reference: STARalignment & HISAT2
#3_STEP. Raw Mapped Reads: SAMtools view; SAMtools sort e SAMtools index

echo "            >>>>>> Starting Pipeline for Pre-processing Variant Calling  <<<<<< $(date) " >> "$pre_processing_DIR/preprocessing_$(date +%Y%m%d).log"

#1-verificar os cromossomos mapeados e presentes nos arquivos BAM baixados (se possuem contigs/cromossomos alternativos) com samtools idxstats
#2-remover cromossomos alternativos dos arquivos bam 
#3-atualizar o cabecalho dos arquivos bam filtrados (apenas com cromossomos padrao (1 a 22, X, Y e M))
#4-indexar os novos arquivos bam criados com samtools index 
#5-rodar a pipeline normalmente 

mkdir -p "$wd"


######### 4_Gerar o TOY : DownsampleSam Picard toolS #############

mkdir -p "$pre_processing_DIR"

stage_DownsampleSam() {
  local bam_file="$1"
  local id
  id=$(basename "$bam_file" .rna_seq.genomic.gdc_realn.bam)

  # Mensagem de início no log
  echo "" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  echo ">>> Iniciando processamento com DownsampleSam: $(date) <<<" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  echo ">>>>>> Processando arquivo: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

  # Execução do comando Downsamples
  if ${PICARD} DownsampleSam \
      I= ${wd}/TOYS_SAMPLES/${ID}.rna_seq.genomic.gdc_realn.bam\
      O="${wd}/TOYS_SAMPLES/${id}.downsampled.bam" \
      P=0.1 >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log" 2>&1; then
    # Mensagem de sucesso
    echo ">>>>>> Processamento concluído para: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  else
    # Mensagem de erro
    echo ">>>>>> Erro ao processar: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  fi
}
export -f stage_DownsampleSam

mkdir -p "${wd}/TOYS_SAMPLES"
find "$BAM_LIST" -name "*.bam" | parallel -j 4 stage_MarkDuplicates {}



######### 4_STEP. SortSam : Picard toolS #############
#caso nao esteja ordenado por coordenada(necessario mara markduplicate)
#stage_SortSam(){
#  local bam_file="$1"
#  local id
#  id=$(basename "$bam_file" .rna_seq.genomic.gdc_realn.bam)
#  
#  echo ">>> Iniciando processamento com SortSam: $(date) <<<" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
#  echo ">>>>>> Processando arquivo: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  
#  if "$PICARD" SortSam \
#     -I "$bam_file" \
#     -O "$pre_processing_DIR/sorted_bam/${id}_sorted.bam" \
#     -SORT_ORDER coordinate
#     2> "$pre_processing_DIR/sorted_bam/${id}_sorted.bam.log"; then
  
    # Mensagem de sucesso
#    echo ">>>>>> Processamento concluído para: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
#  else
    # Mensagem de erro
#    echo ">>>>>> Erro ao processar: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
#  fi
  
#}
  
#export -f stage_SortSam
#mkdir -p "$pre_processing_DIR/sorted_bam"
#find "${wd}/5samples_BAMs" -name "*.bam" | parallel -j 4 stage_SortSam {}


######## 5_STEP. MarkDuplicates : Picard tool ############

stage_MarkDuplicates() {
  local bam_file="$1"
  local id
  id=$(basename "$bam_file" .rna_seq.genomic.gdc_realn.bam)

  # Mensagem de início no log
  echo "" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  echo ">>> Iniciando processamento com MarkDuplicates: $(date) <<<" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  echo ">>>>>> Processando arquivo: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

  # Execução do comando MarkDuplicates
  if "$PICARD" MarkDuplicates \
      -I "$pre_processing_DIR/TOYS_SAMPLES/${id}.downsampled.bam" \
      -O "$pre_processing_DIR/marked_duplicates/${id}_marked_duplicates.bam" \
      -M "$pre_processing_DIR/marked_duplicates/${id}_marked_dup_metrics.txt" \
      2> "$pre_processing_DIR/marked_duplicates/${id}_marked_duplicates.bam.log"; then

    # Mensagem de sucesso
    echo ">>>>>> Processamento concluído para: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  else
    # Mensagem de erro
    echo ">>>>>> Erro ao processar: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
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

  echo "" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  echo ">>>>>> Iniciando SplitNCigarRead para Amostra: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

  $GATK --java-options "-Xmx${MEM}G" SplitNCigarReads \
      -R "${REF_FASTA}" \
      -I "${pre_processing_DIR}/marked_duplicates/${id}_marked_duplicates.bam" \
      -O "${pre_processing_DIR}/SplitNCigarRead/${id}_marked_duplicates.tags.bam" \
      --tmp-dir "${pre_processing_DIR}/tmp_GATK" \
      > "$pre_processing_DIR/SplitNCigarRead/${id}_marked_duplicates.tags.bam.log" \
      2> "$pre_processing_DIR/SplitNCigarRead/${id}_marked_duplicates.tags.bam.log2"
      
    if [[ $? -eq 0 ]]; then
        echo ">>>>>> SplitNCigarReads concluído com sucesso para: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
    else
        echo ">>>>>> Erro ao executar SplitNCigarReads para: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
    fi
}

export -f SplitNCigarRead

mkdir -p "$pre_processing_DIR/SplitNCigarRead"
mkdir -p "$pre_processing_DIR/tmp_GATK"
find "$BAM_LIST" -name "*.bam" | parallel -j 4 SplitNCigarRead {}


######### 7_STEP. Recalibrate Base Quality Score(BQSR): BaseRecalibrator##########

stage_base_recalibration(){
  local sample="$1"
  local id
  id=$(basename "$sample" .rna_seq.genomic.gdc_realn.bam)
  
  echo "" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
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




#######################################
############# GATK-MUTECT #############
#######################################

echo " >>>>>> Starting Pipeline to Run GATK-MUTECT2  <<<<<< $(date) " > "${LOG_FILE}"



################# 1_STEP. Running Mutect2 #####################

stage_Mutect2 (){
  local sample="$1"
  local id
  id=$(basename "$sample" .rna_seq.genomic.gdc_realn.bam)
  
  echo "" >> "$LOG_FILE"
  echo ">>>>>> Executando stage_Mutect2 para Amostra: "$id" <<<  $(date) " >> "$LOG_FILE"

  $GATK --java-options "-Xmx${MEM}G" Mutect2 \
        -R ${REF_FASTA} \
        -L ${TARGET} \
        -I ${pre_processing_DIR}/Filtered/"${id}_filtered_marked_duplicates.tags.recal.new.bam" \
        --germline-resource ${GNOMAD} \
        --panel-of-normals ${PON} \
        --f1r2-tar-gz ${OUTPUT_DIR}/Mutect2/"${id}.f1r2.tar.gz" \
        --bam-output ${OUTPUT_DIR}/Mutect2/"${id}.bamout.bam" \
        -O ${OUTPUT_DIR}/Mutect2/"${id}.unfiltered.vcf.gz" 2> $OUTPUT_DIR/Mutect2/"${id}.unfiltered.log"
        
   if [[ $? -eq 0 ]]; then
      echo ">>>>>> stage_Mutect2 concluído com sucesso para: $id <<<<<< $(date)" >> "$LOG_FILE"
  else
      echo ">>>>>> Erro ao executar stage_Mutect2 para: $id <<<<<< $(date)" >> "$LOG_FILE"
      return 1
  fi
}
export -f stage_Mutect2

mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}/Mutect2
mkdir -p ${LOG_DIR}
find "$BAM_LIST" -name "*.bam" | parallel -j 4 stage_Mutect2 {}




########### 2_STEP. LearnReadOrientationModel ###############
#gera uma tabela de prior de artefatos para cada amostra de tumor, a ser usada pelo FilterMutectCalls.

stage_LearnReadOrientationModel (){
  echo "" >> "$LOG_FILE"
  echo ">>>>>> STAGE_LearnReadOrientationModel<<<<<<  $(date) " >> "$LOG_FILE"

  local SAMPLE_READORIENTATION=$(find "$OUTPUT_DIR"/Mutect2/ -maxdepth 1 -mindepth 1  -name '*.f1r2.tar.gz')
  local SAMPLE_F1R2=$(echo $SAMPLE_READORIENTATION| sed 's/\s/ -I  /g')

  $GATK --java-options  "-Xmx${MEM}G"  LearnReadOrientationModel \
        -I $SAMPLE_F1R2 \
        -O "$OUTPUT_DIR"/LearnReadOrientationModel/read-orientation-model.tar.gz \
        2> $OUTPUT_DIR/LearnReadOrientationModel/read-orientation-model.log
  
  if [[ $? -ne 0 ]]; then
      echo "Erro ao executar LearnReadOrientationModel." >> "$LOG_FILE"
      exit 1
  fi
} 
export -f stage_LearnReadOrientationModel

mkdir -p "$OUTPUT_DIR/LearnReadOrientationModel"
stage_LearnReadOrientationModel


########### 3_STEP. GetPileupSummaries ###############
# Resume as contagens de leituras: Extrai informações sobre as frequências alélicas em regiões conhecidas de variantes germinativas 
# output utilizado como input para calculatecontamination 

stage_GetPileupSummaries (){
  local sample="$1"
  local id
  id=$(basename "$sample" .rna_seq.genomic.gdc_realn.bam)
  
  echo "" >> "$LOG_FILE"
  echo ">>>>>> Executando GetPileupSummaries para Amostra: "$id" <<<  $(date) " >> "$LOG_FILE"

  $GATK --java-options "-Xmx${MEM}G" GetPileupSummaries \
        -I ${pre_processing_DIR}/Filtered/"${id}_filtered_marked_duplicates.tags.recal.new.bam" \
        -V $GNOMAD \
        -L ${TARGET} \
        -O $OUTPUT_DIR/GetPileupSummaries/"${id}.getpileupsummaries.table" \
        2> $OUTPUT_DIR/GetPileupSummaries/"${id}".getpileupsummaries.log
      
  if [[ $? -eq 0 ]]; then
      echo ">>>>>> GetPileupSummaries concluído com sucesso para: $id <<<<<< $(date)" >> "$LOG_FILE"
  else
      echo ">>>>>> Erro ao executar GetPileupSummaries para: $id <<<<<< $(date)" >> "$LOG_FILE"
      return 1
  fi
}
export -f stage_GetPileupSummaries

mkdir -p $OUTPUT_DIR/GetPileupSummaries
find "$BAM_LIST" -name "*.bam" | parallel -j 2 stage_GetPileupSummaries {}
#se reduz o jobs do parallel, pode aumentar a memoria RAM 


########### 4_STEP. CalculateContamination ###############

stage_CalculateContamination (){
  local sample="$1"
  local id
  id=$(basename "$sample" .rna_seq.genomic.gdc_realn.bam)
  
  echo "" >> "$LOG_FILE"
  echo ">>>>>> Executando CalculateContamination para Amostra: "$id" <<<  $(date) " >> "$LOG_FILE"

  $GATK --java-options "-Xmx${MEM}G" CalculateContamination \
        -I $OUTPUT_DIR/GetPileupSummaries/"${id}.getpileupsummaries.table" \
        -tumor-segmentation $OUTPUT_DIR/CalculateContamination/"${id}.segments.table" \
        -O $OUTPUT_DIR/CalculateContamination/"${id}.calculatecontamination.table" \
        2> $OUTPUT_DIR/CalculateContamination/"${id}.calculatecontamination.log"

  if [[ $? -eq 0 ]]; then
      echo ">>>>>> CalculateContamination concluído com sucesso para: $id <<<<<< $(date)" >> "$LOG_FILE"
  else
      echo ">>>>>> Erro ao executar CalculateContamination para: $id <<<<<< $(date)" >> "$LOG_FILE"
      return 1
  fi
}
export -f stage_CalculateContamination

mkdir -p $OUTPUT_DIR/CalculateContamination/
find "$BAM_LIST" -name "*.bam" | parallel -j 4 stage_CalculateContamination {}



########### 5_STEP. FilterMutectCalls ###############

stage_FilterMutectCalls (){
  local sample="$1"
  local id
  id=$(basename "$sample" .rna_seq.genomic.gdc_realn.bam)
  
  echo "" >> "$LOG_FILE"
  echo ">>>>>> Executando FilterMutectCalls para Amostra: "$id" <<<  $(date) " >> "$LOG_FILE"
 
  $GATK --java-options "-Xmx${MEM}G" FilterMutectCalls \
        -R $REF_FASTA \
        -V $OUTPUT_DIR/Mutect2/"${id}.unfiltered.vcf.gz" \
        --tumor-segmentation $OUTPUT_DIR/CalculateContamination/"${id}.segments.table" \
        --contamination-table $OUTPUT_DIR/CalculateContamination/"${id}.calculatecontamination.table" \
        --stats $OUTPUT_DIR/Mutect2/"${id}.unfiltered.vcf.gz.stats" \
        --ob-priors "$OUTPUT_DIR"/LearnReadOrientationModel/read-orientation-model.tar.gz \
        -O $OUTPUT_DIR/FilterMutectCalls/"${id}.filtered.vcf.gz" \
        2> $OUTPUT_DIR/FilterMutectCalls/"${id}.filtered.vcf.log"

  if [[ $? -eq 0 ]]; then
      echo ">>>>>> FilterMutectCalls concluído com sucesso para: $id <<<<<< $(date)" >> "$LOG_FILE"
  else
      echo ">>>>>> Erro ao executar FilterMutectCalls para: $id <<<<<< $(date)" >> "$LOG_FILE"
      return 1
  fi
}
export -f stage_FilterMutectCalls

mkdir -p $OUTPUT_DIR/FilterMutectCalls/
find "$BAM_LIST" -name "*.bam" | parallel -j 4 stage_FilterMutectCalls {}


#######################################
######### FILTROS ADCIONAIS ###########
#######################################


########### 1_STEP. left normalization: BCFtools ###############

left_normalization () {
  local sample="$1"
  local id
  id=$(basename "$sample" .rna_seq.genomic.gdc_realn.bam)
  
  echo "" >> "$LOG_FILE"
  echo ">>>>>> Executando Normalization para Amostra: "$id" <<<" >> "$LOG_FILE"
 
  
  bcftools norm -m-both -O z -o $OUTPUT_DIR/left_normalization/"${id}.norm_Step1.vcf.gz" \
    $OUTPUT_DIR/FilterMutectCalls/"${id}.filtered.vcf.gz" \
    2> $OUTPUT_DIR/left_normalization/"${id}.norm_Step1.log"
    
  if [[ $? -eq 0 ]]; then
    echo ">>>>>> bcftools norm Step1 concluído com sucesso para Amostra: $id <<< $(date)" >> "$LOG_FILE"
  else
    echo ">>>>>> Erro ao executar bcftools norm Step1 para Amostra: $id <<< $(date)" >> "$LOG_FILE"
    return 1
  fi
    
  bcftools norm -O z -f ${REF_FASTA} -o $OUTPUT_DIR/left_normalization/"${id}.norm_Step2.vcf.gz" \
    $OUTPUT_DIR/left_normalization/"${id}.norm_Step1.vcf.gz" \
    2> $OUTPUT_DIR/left_normalization/"${id}.norm_Step2.log"
    
  if [[ $? -eq 0 ]]; then
    echo ">>>>>> bcftools norm Step2 concluído com sucesso para Amostra: $id <<< $(date)" >> "$LOG_FILE"
  else
    echo ">>>>>> Erro ao executar bcftools norm Step2 para Amostra: $id <<< $(date)" >> "$LOG_FILE"
    return 1
  fi
    
  bcftools index $OUTPUT_DIR/left_normalization/"${id}.norm_Step2.vcf.gz"
  2> $OUTPUT_DIR/left_normalization/"${id}.index.log"

  
  if [[ $? -eq 0 ]]; then
    echo ">>>>>> bcftools index concluído com sucesso para Amostra: $id <<< $(date)" >> "$LOG_FILE"
  else
    echo ">>>>>> Erro ao executar bcftools index para Amostra: $id <<< $(date)" >> "$LOG_FILE"
    return 1
  fi
}
export -f left_normalization

mkdir -p $OUTPUT_DIR/left_normalization
find "$BAM_LIST" -name "*.bam" | parallel -j 4 left_normalization {}


selectPASS(){
  local sample="$1"
  local id
  id=$(basename "$sample" .rna_seq.genomic.gdc_realn.bam)
  
  # Mensagem de início no log geral
  echo "" >> "$LOG_FILE"
  echo "> Executando selectPASS para Amostra: "$id" <<< $(date) " >> "$LOG_FILE"
  
  # Filtrar variantes com flag PASS
  bcftools view -f PASS -O z $OUTPUT_DIR/left_normalization/"${id}.norm_Step2.vcf.gz" > \
    $OUTPUT_DIR/FILTER_PASS/"${id}.pass.vcf.gz" \
    2> $OUTPUT_DIR/FILTER_PASS/"${id}.pass.log"

  # Indexar o VCF filtrado
  bcftools index $OUTPUT_DIR/FILTER_PASS/"${id}.pass.vcf.gz" \
    2>> $OUTPUT_DIR/FILTER_PASS/"${id}.pass.log"

  # Verificar sucesso e registrar no log geral
  if [[ $? -eq 0 ]]; then
    echo "> selectPASS concluído com sucesso para Amostra: $id <<< $(date)" >> "$LOG_FILE"
  else
    echo "> Erro ao executar selectPASS para Amostra: $id <<< $(date)" >> "$LOG_FILE"
  fi
}
export -f selectPASS


mkdir -p $OUTPUT_DIR/FILTER_PASS
find "$BAM_LIST" -name "*.bam" | parallel -j 4 selectPASS {}


selectVar(){
  local sample="$1"
  local id
  id=$(basename "$sample" .rna_seq.genomic.gdc_realn.bam)
  echo "> Executando selectVar para Amostra: "$id" <<< $(date) " >> "$LOG_FILE"
  zgrep "#" $OUTPUT_DIR/FILTER_PASS/"${id}.pass.vcf.gz" > $OUTPUT_DIR/FILTER_VAR/"${id}.pass.var.vcf"
  zgrep -v "#" $OUTPUT_DIR/FILTER_PASS/"${id}.pass.vcf.gz" |  grep -v -E "0/0:|0/\.|\.\/\." >> $OUTPUT_DIR/FILTER_VAR/"${id}.pass.var.vcf"
  bgzip $OUTPUT_DIR/FILTER_VAR/"${id}.pass.var.vcf" 
  bcftools index $OUTPUT_DIR/FILTER_VAR/"${id}.pass.var.vcf.gz" 
  
# Verificar sucesso e registrar no log geral
  if [[ $? -eq 0 ]]; then
    echo "> selectPASS concluído com sucesso para Amostra: $id <<< $(date)" >> "$LOG_FILE"
  else
    echo "> Erro ao executar selectPASS para Amostra: $id <<< $(date)" >> "$LOG_FILE"
  fi
}
export -f selectVar

mkdir -p $OUTPUT_DIR/FILTER_VAR
find "$BAM_LIST" -name "*.bam" | parallel -j 4 selectVar {}

################ EXCLUIR RNAeditind ##############
#database: REDIportal

Filter_rna() {
  local sample="$1"
  local id
  id=$(basename "$sample" .rna_seq.genomic.gdc_realn.bam)
  local BASE_PATH="/home/scratch90/bleocata_20241205/MMRF-COMMPASS/Variant_Calling"

  # Crie saída filtrada com base no arquivo de variantes
  zcat "$OUTPUT_DIR/FILTER_VAR/${id}.pass.var.vcf.gz" | \
  awk -v exclude_file="${BASE_PATH}/RNAediting_hg38_v2.txt" 'BEGIN {
      while (getline < exclude_file) exclude[$0]=1
  }
  /^#/ { print; next }
  {
      key = $1":"$2":"$4":"$5
      if (!(key in exclude)) print
  }' | gzip > "$OUTPUT_DIR/Filter_RNAediting/${id}.pass.var.rna.vcf.gz"
  
  
  if [[ $? -eq 0 ]]; then
    echo "> Filter RNA editing concluído com sucesso para Amostra: $id <<< $(date)" >> "$LOG_FILE"
  else
    echo "> Erro ao executar Filter RNA editing para Amostra: $id <<< $(date)" >> "$LOG_FILE"
  fi
}
export -f Filter_rna

mkdir -p "$OUTPUT_DIR/Filter_RNAediting"
find "$BAM_LIST" -name "*.bam" | parallel -j 4 Filter_rna {}




#######################################
############ ANNOTATION  ##############
#######################################

annotation (){
  local sample="$1"
  local id
  id=$(basename "$sample" .rna_seq.genomic.gdc_realn.bam)
  
  echo "" >> "$LOG_FILE"
  echo ">>>>>> Executando annovar para amostra $id: <<<<<< $(date)" >> "$LOG_FILE"
  
  $ANNOVAR  \
    --vcfinput $OUTPUT_DIR/FILTER_VAR/"${id}.pass.var.vcf.gz" \
    $ANNOVAR_DB -buildver hg38 --remove \
    --protocol refGene,avsnp150,gnomad41_exome_filt,abraom,cosmic99,icgc28,dbnsfp42a_filt,clinvar_20220320  \
    --operation gx,f,f,f,f,f,f,f --arg '-splicing 5',,,,,,, --polish \
    --xreffile $CROSS_REFERENCE --otherinfo --thread 10 \
    --outfile $OUTPUT_DIR/annotation/${id} \
    > $OUTPUT_DIR/annotation/${id}.log \
    2> $OUTPUT_DIR/annotation/${id}.log2
    
  if [[ $? -eq 0 ]]; then
    echo ">>>>>> ANNOVAR concluído com sucesso para Amostra: $id <<< $(date)" >> "$LOG_FILE"
  else
    echo ">>>>>> Erro ao executar ANNOVAR para Amostra: $id <<< $(date)" >> "$LOG_FILE"
    return 1
  fi
  
  # Corrigir caracteres no arquivo VCF anotadccco
  sed 's/\\x3b/;/g' $OUTPUT_DIR/annotation/${id}.hg38_multianno.vcf | sed 's/\\x3d/=/g' > $OUTPUT_DIR/annotation/${id}.hg38_multianno.correct.vcf  ## pelo jeirto nao precisava

 # Verificar se o arquivo corrigido foi gerado
  if [[ $? -eq 0 && -f $OUTPUT_DIR/annotation/${id}.hg38_multianno.correct.vcf ]]; then
    echo ">>>>>> Correção do VCF concluída para Amostra: $id <<< $(date)" >> "$LOG_FILE"
  else
    echo ">>>>>> Erro ao corrigir o VCF para Amostra: $id <<< $(date)" >> "$LOG_FILE"
    return 1
  fi
  
}
export -f annotation

mkdir -p $OUTPUT_DIR/annotation
find "$BAM_LIST" -name "*.bam" | parallel -j 4 annotation {}


################# Manipular as tabelas do annovar ###################

# Extrair o cabeçalho do primeiro arquivo
first_file=$(find $OUTPUT_DIR/annotation -name "*.hg38_multianno.txt" | head -n 1)
header_file="$OUTPUT_DIR/header.tmp"
if [ -n "$first_file" ]; then
  # Extrair o cabeçalho e adicionar "sample" como a primeira coluna
  head -n 1 "$first_file" | awk '{print "sample\t" $0}' > "$header_file"
else
  echo "Nenhum arquivo encontrado no diretório $OUTPUT_DIR/annotation."
  exit 1
fi


# Função para processar arquivos sem incluir o cabeçalho
concatenar() {
  local sample="$1"
  local id=$(basename "$sample" .hg38_multianno.txt)
  local temp_output="$OUTPUT_DIR/Final_GATK.4.6_annotated_$id.tmp"
  
  echo "" >> "$LOG_FILE"
  echo ">>>>>> Processando concatenamento para amostra $id <<<<<< $(date)" >> "$LOG_FILE"
  
  # Adicionar ID ao arquivo e remover o cabeçalho (linha 1)
  awk -OFS="\t" -v N="$id" 'NR > 1 {print N, $0}' "$sample" > "$temp_output"
  sed -i 's/\s\+/\t/g' "$temp_output"
}
export -f concatenar

# Processar os arquivos em paralelo
find $OUTPUT_DIR/annotation -name "*.hg38_multianno.txt" | parallel -j 4 concatenar {}

# Concatenar o cabeçalho e os arquivos processados
cat "$header_file" $OUTPUT_DIR/Final_GATK.4.6_annotated_*.tmp > $OUTPUT_DIR/Final_GATK.4.6_annotated.txt

# Limpar arquivos temporários
rm "$header_file"
rm $OUTPUT_DIR/Final_GATK.4.6_annotated_*.tmp
