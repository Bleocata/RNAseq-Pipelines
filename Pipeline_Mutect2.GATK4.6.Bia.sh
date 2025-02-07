#Pipeline atualizado em 14/06/2024
#1. Ajustes das variaves de ambiente para facilitar na reexecução quando os diretorios forem atualizado
#2. Alterado o diretorio onde estava o PON na scratch45 para /home/users/vlira/PanelOfNormals/PoN.100COVID.vcf.gz
#3. Agora o script recebe 2 agumentos: 1- lista de amostras; 2- diretorio SCRATCH
#4. Utualizado os databases: hg38_cosmic98_coding,hg38_avsnp150,hg38_clinvar_20220320, hg38_gnomad40_exome

############# Definicao de variaveis #################

export wd="/home/scratch90/bleocata_20241205/MMRF-COMMPASS/Variant_Calling"
export path="/home/scratch90/vlira_18nov2024"

export MEM=200
export JOBS=5

export pre_processing_DIR="${wd}/preprocessing_result"
export OUTPUT_DIR="${wd}/Result_Mutect2.GATK4.6"
export LOG_DIR="${wd}/Arquivos_log"
export LOG_FILE="${LOG_DIR}/Pipeline_Mutect2_$(date +%Y%m%d).log"

export REF_FASTA="/home/projects2/LIDO/molPathol/oncoseek/nextseq/hg38"
export GATK="$path/tools/gatk-4.6.0.0/./gatk"
export PICARD="java -jar ${path}tools/picard-3.2.0/picard.jar"

export TARGET="$path/references/xgen-exome-research-panel-v2-targets-hg38.bed"
export PON="${wd}/1000g_pon.hg38.vcf.gz"
export GNOMAD="$path/references/af-only-gnomad.hg38.vcf.gz" #frequência de variantes genéticas em populações diversas
#export GNOMAD="$wd/references/af-only-gnomad.SABE1171.Abraom.hg38.vcf.gz"
export ANNOVAR="$path/tools/annovar/table_annovar.pl"
export ANNOVAR_DB="$path/humandb/"
export CROSS_REFERENCE="$path/references/refGene_TARGET_COSMICv82CensusGene_F1.txt"

export INDEL_KNOWN="/home/projects2/LIDO/molPathol/oncoseek/nextseq/references/Mills_and_1000G_gold_standard.indels.hg38.vcf"
export DBSNP="$wd/references/dbsnp_151.hg38.vcf.gz"


#######################################
########### PRE-PROCECING #############
#######################################
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery

#1_STEP. Raw Unmapped Reads 
#2_STEP. Map to Reference: STARalignment & HISAT2
#3_STEP. Raw Mapped Reads: SAMtools view; SAMtools sort e SAMtools index


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
  echo ">>> Iniciando processamento com MarkDuplicates: $(date) <<<" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  echo ">>>>>> Processando arquivo: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

  # Execução do comando MarkDuplicates
  if $PICARD MarkDuplicates \
      -I "$bam_file" \
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



######### 6_STEP. SplitNCigarRead: GATK Tools ############
#Adicionar ou recalcular os campos NM (número de diferenças), MD (mismatch descriptor), e UQ (qualidade única de alinhamento) nos alinhamentos de um arquivo BAM usando a ferramenta SetNmMdAndUqTags do Picard Tools.

echo ">>>>>> SplitNCigarRead: GATK Tools <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

SplitNCigarRead() {
  local sample="$1"
  local id
  id=$(basename "$sample" _marked_duplicates.bam)

  echo ">>>>>> Iniciando SplitNCigarRead para Amostra: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

  "$GATK" --java-options "-Xmx${MEM}G" SplitNCigarReads \
      -R "${REF_FASTA}/Homo_sapiens_assembly38.fasta" \
      -I "$sample" \
      -O "$pre_processing_DIR/SplitNCigarRead/${id}_marked_duplicates.tags.bam" \
      2> "$pre_processing_DIR/SplitNCigarRead/${id}_marked_duplicates.tags.bam.log"
      
    if [[ $? -eq 0 ]]; then
        echo ">>>>>> SplitNCigarReads concluído com sucesso para: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
    else
        echo ">>>>>> Erro ao executar SplitNCigarReads para: $id <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
    fi
}

export -f SplitNCigarRead


######### 7_STEP. Recalibrate Base Quality Score(BQSR): BaseRecalibrator##########

stage_base_recalibration(){
  local sample="$1"
  local id
  id=$(basename "$sample" _marked_duplicates.tags.bam)
  echo ">>>>>> Executando stage_base_recalibration para Amostra: "$id" <<<  $(date) " >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

  $GATK --java-options "-Xmx${MEM}G" BaseRecalibrator \
        -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
        -I ${pre_processing_DIR}/"${id}._marked_duplicates.tags.bam" \
        -known-sites ${INDEL_KNOWN} \
        -known-sites ${DBSNP} \
        -known-sites ${GNOMAD} \
        -L ${TARGET} \
        -O ${pre_processing_DIR}/"${id}_marked_duplicates.tags.bqsr.table" \
        > ${pre_processing_DIR}/"${id}_marked_duplicates.tags.bqsr.table.log" \
        2> ${pre_processing_DIR}/"${id}_marked_duplicates.tags.table.log2"

  $GATK --java-options "-Xmx${MEM}G" ApplyBQSR \
        -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
        -I $wd/preprocessing_result/"${id}_marked_duplicates.tags.bam" \
        --bqsr-recal-file $wd/preprocessing_result/"${id}_marked_duplicates.tags.bqsr.table"\
        -O ${pre_processing_DIR}/"${id}_marked_duplicates.tags.bqsr.bam" \
        > ${pre_processing_DIR}/"${id}_marked_duplicates.tags.bqsr.bam.log" \
        2> ${pre_processing_DIR}/"${id}_marked_duplicates.tags.bqsr.bam.log2"
        
  $GATK --java-options "-Xmx${MEM}G" AnalyzeCovariates \
       -I my_reads.bam \
       -R reference.fasta \
       --known-sites sites_of_variation.vcf \
       --known-sites another/optional/setOfSitesToMask.vcf \
       -O recal_data.table
}
export -f  stage_base_recalibration




echo "            >>>>>> Starting Pipeline for Pre-processing Variant Calling  <<<<<< $(date) " >> "$pre_processing_DIR/preprocessing_$(date +%Y%m%d).log"

mkdir -p "$pre_processing_DIR"

mkdir -p "$pre_processing_DIR/marked_duplicates"
find "${wd}/5samples_BAMs" -name "*.bam" | parallel -j 4 stage_MarkDuplicates {}

mkdir -p "$pre_processing_DIR/SplitNCigarRead"
find "$pre_processing_DIR/marked_duplicates" -name "*_marked_duplicates.bam" | parallel -j 4 SplitNCigarRead {}

mkdir -p "$pre_processing_DIR/BQSR_base_recalibration"
find "$pre_processing_DIR/SplitNCigarRead" -name "*_marked_duplicates.tags.bam" | parallel -j 4 stage_base_recalibration {}



#######################################
############# GATK-MUTECT #############
#######################################



################# 1_STEP. Running Mutect2 #####################

stage_Mutect2 (){
  local SAMPLE=$1
  id="${SAMPLE##*/}"
  echo "" >> "$LOG_FILE"
  echo ">>>>>> Executando stage_Mutect2 para Amostra: "$id" <<<  $(date) " >> "$LOG_FILE"

  $GATK --java-options "-Xmx${MEM}G" Mutect2 \
        -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
        -L $TARGET \
        -I $wd/preprocessing_result/"${id%.dedup*}.dedup.tags.bqsr.bam" \
        --germline-resource $GNOMAD \
        --panel-of-normals $PON   \
        --f1r2-tar-gz $OUTPUT_DIR/Mutect2/"${id%.dedup*}.f1r2.tar.gz" \
        -bamout $OUTPUT_DIR/Mutect2/"${id%.dedup*}.bamout.bam" \
        -O $OUTPUT_DIR/Mutect2/"${id%.dedup*}.unfiltered.vcf.gz" 2> $OUTPUT_DIR/Mutect2/"${id%.dedup*}.unfiltered.log"
}
export -f stage_Mutect2



########### 2_STEP. LearnReadOrientationModel ###############
#gera uma tabela de prior de artefatos para cada amostra de tumor, a ser usada pelo FilterMutectCalls.

stage_LearnReadOrientationModel (){
  echo "" >> "$LOG_FILE"
  echo ">>>>>> STAGE_LearnReadOrientationModel<<<<<<  $(date) " >> "$LOG_FILE"

  local SAMPLE_READORIENTATION=$(find "$OUTPUT_DIR"/Mutect2/ -maxdepth 1 -mindepth 1  -name '*.f1r2.tar.gz')
  local SAMPLE_F1R2=$(echo $SAMPLE_READORIENTATION| sed 's/\s/ -I  /g')

  $GATK --java-options  "-Xmx${MEM}G"  LearnReadOrientationModel \
        -I $SAMPLE_F1R2 \
        -O "$OUTPUT_DIR"/LearnReadOrientationModel/read-orientation-model.tar.gz 2> $OUTPUT_DIR/LearnReadOrientationModel/read-orientation-model.tar.gz.log
} 
export -f stage_LearnReadOrientationModel


########### 3_STEP. GetPileupSummaries ###############

stage_GetPileupSummaries (){
  local SAMPLE=$1
  id="${SAMPLE##*/}"
  echo "" >> "$LOG_FILE"
  echo ">>>>>> Executando GetPileupSummaries para Amostra: "$id" <<<  $(date) " >> "$LOG_FILE"

  $GATK --java-options -Xmx300G GetPileupSummaries \
        -I $wd/preprocessing_result/"${id%.dedup*}.dedup.tags.bqsr.bam" \
        -V $GNOMAD \
        -L $GNOMAD \
        -O $OUTPUT_DIR/GetPileupSummaries/"${id%.dedup*}.getpileupsummaries.table" 2> $OUTPUT_DIR/GetPileupSummaries/"${id%.dedup*}".getpileupsummaries.log
}
export -f stage_GetPileupSummaries


########### 4_STEP. CalculateContaminationries ###############

stage_CalculateContamination (){
  local SAMPLE=$1
  id="${SAMPLE##*/}"
  echo "" >> "$LOG_FILE"
  echo ">>>>>> Executando CalculateContamination para Amostra: "$id" <<<  $(date) " >> "$LOG_FILE"

  $GATK --java-options "-Xmx${MEM}G" CalculateContamination \
        -I $OUTPUT_DIR/GetPileupSummaries/"${id%.dedup*}.getpileupsummaries.table" \
        -tumor-segmentation $OUTPUT_DIR/CalculateContamination/"${id%.dedup*}.segments.table" \
        -O $OUTPUT_DIR/CalculateContamination/"${id%.dedup*}.calculatecontamination.table" 2> $OUTPUT_DIR/CalculateContamination/"${id%.dedup*}.calculatecontamination.log"
}
export -f stage_CalculateContamination



########### 5_STEP. FilterMutectCalls ###############

stage_FilterMutectCalls (){
  local SAMPLE=$1
  id="${SAMPLE##*/}"
  echo "" >> "$LOG_FILE"
  echo ">>>>>> Executando FilterMutectCalls para Amostra: "$id" <<<  $(date) " >> "$LOG_FILE"
 
  $GATK --java-options "-Xmx${MEM}G" FilterMutectCalls \
        -R $REF_FASTA/Homo_sapiens_assembly38.fasta \
        -V $OUTPUT_DIR/Mutect2/"${id%.dedup*}.unfiltered.vcf.gz" \
        --tumor-segmentation $OUTPUT_DIR/CalculateContamination/"${id%.dedup*}.segments.table" \
        --contamination-table $OUTPUT_DIR/CalculateContamination/"${id%.dedup*}.calculatecontamination.table" \
        --stats $OUTPUT_DIR/Mutect2/"${id%.dedup*}.unfiltered.vcf.gz.stats" \
        --ob-priors "$OUTPUT_DIR"/LearnReadOrientationModel/read-orientation-model.tar.gz \
        -O $OUTPUT_DIR/FilterMutectCalls/"${id%.dedup*}.filtered.vcf.gz" 2> $OUTPUT_DIR/FilterMutectCalls/"${id%.dedup*}.filtered.vcf.log"
}
export -f stage_FilterMutectCalls


left_normalization () {
  local SAMPLE=$1
  id="${SAMPLE##*/}"
  echo "" >> "$LOG_FILE"
  echo ">>>>>> Executando Normalization para Amostra: "$id" <<<" >> "$LOG_FILE"
  date >> "$LOG_FILE"
  bcftools norm -m-both -O z -o $OUTPUT_DIR/left_normalization/"${id%.dedup*}.norm_Step1.vcf.gz" \
    $OUTPUT_DIR/FilterMutectCalls/"${id%.dedup*}.filtered.vcf.gz" 2> $OUTPUT_DIR/left_normalization/"${id%.dedup*}.norm_Step1.log"
  bcftools norm -O z -f $REF_FASTA/Homo_sapiens_assembly38.fasta -o $OUTPUT_DIR/left_normalization/"${id%.dedup*}.norm_Step2.vcf.gz" \
    $OUTPUT_DIR/left_normalization/"${id%.dedup*}.norm_Step1.vcf.gz" 2> $OUTPUT_DIR/left_normalization/"${id%.dedup*}.norm_Step2.log"
  bcftools index $OUTPUT_DIR/left_normalization/"${id%.dedup*}.norm_Step2.vcf.gz"
}
export -f left_normalization


selectPASS(){
  local SAMPLE=$1
  id="${SAMPLE##*/}"
  echo "> Executando selectPASS para Amostra: "$id" <<< $(date) " >> "$LOG_FILE"
  bcftools view -f PASS -O z $OUTPUT_DIR/left_normalization/"${id%.dedup*}.norm_Step2.vcf.gz" > $OUTPUT_DIR/FILTER_PASS/"${id%.dedup*}.pass.vcf.gz"
  bcftools index $OUTPUT_DIR/FILTER_PASS/"${id%.dedup*}.pass.vcf.gz" 
}
export -f selectPASS


selectVar(){
  local SAMPLE=$1
  id="${SAMPLE##*/}"
  echo "> Executando selectVar para Amostra: "$id" <<< $(date) " >> "$LOG_FILE"
  zgrep "#" $OUTPUT_DIR/FILTER_PASS/"${id%.dedup*}.pass.vcf.gz" > $OUTPUT_DIR/FILTER_VAR/"${id%.dedup*}.pass.var.vcf"
  zgrep -v "#" $OUTPUT_DIR/FILTER_PASS/"${id%.dedup*}.pass.vcf.gz" |  grep -v -E "0/0:|0/\.|\.\/\." >> $OUTPUT_DIR/FILTER_VAR/"${id%.dedup*}.pass.var.vcf"
  bgzip $OUTPUT_DIR/FILTER_VAR/"${id%.dedup*}.pass.var.vcf" 
  bcftools index $OUTPUT_DIR/FILTER_VAR/"${id%.dedup*}.pass.var.vcf.gz" 
}
export -f selectVar


annotation (){
  local SAMPLE=$1
  id="${SAMPLE##*/}"
  # # annovar
  echo "" >> "$LOG_FILE"
  echo ">>>>>> Executando annovar para amostra $id: <<<<<< $(date)" >> "$LOG_FILE"
  
  $ANNOVAR  \
    --vcfinput $OUTPUT_DIR/FILTER_VAR/"${id%.dedup*}.pass.var.vcf.gz"  $ANNOVAR_DB -buildver hg38 --remove \
    --protocol refGene,avsnp150,gnomad41_exome_filt,abraom,cosmic99,icgc28,dbnsfp42a_filt,clinvar_20220320  \
    --operation gx,f,f,f,f,f,f,f --arg '-splicing 5',,,,,,, --polish \
    --xreffile $CROSS_REFERENCE --otherinfo --thread 10 \
    --outfile $OUTPUT_DIR/annotation/${id%.dedup*} > $OUTPUT_DIR/annotation/${id%.dedup*}.log 2> $OUTPUT_DIR/annotation/${id%.dedup*}.log2
  sed 's/\\x3b/;/g' $OUTPUT_DIR/annotation/${id%.dedup*}.hg38_multianno.vcf | sed 's/\\x3d/=/g' > $OUTPUT_DIR/annotation/${id%.dedup*}.hg38_multianno.correct.vcf 
}
export -f annotation



echo "            >>>>>> Starting Pipeline to Run GATK-MUTECT2  <<<<<< $(date) " >> "$LOG_FILE"

mkdir -p $OUTPUT_DIR

mkdir -p $OUTPUT_DIR/Mutect2/
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P5 bash -c 'stage_Mutect2  "$@"' 'stage_Mutect2' &

mkdir $OUTPUT_DIR/LearnReadOrientationModel/
stage_LearnReadOrientationModel 

mkdir -p $OUTPUT_DIR/GetPileupSummaries/
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P1 bash -c 'stage_GetPileupSummaries "$@"' 'stage_GetPileupSummaries' &

mkdir -p $OUTPUT_DIR/CalculateContamination/
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P10 bash -c 'stage_CalculateContamination  "$@"' 'stage_CalculateContamination'

mkdir -p $OUTPUT_DIR/FilterMutectCalls/
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P10 bash -c 'stage_FilterMutectCalls  "$@"' 'stage_FilterMutectCalls'

mkdir -p $OUTPUT_DIR/left_normalization/
xargs -a ${SAMPLE_LIST_BAM}  -t -n1 -P10 bash -c 'left_normalization  "$@"' 'left_normalization'

mkdir -p $OUTPUT_DIR/FILTER_PASS
echo ">>>>>> selectPASS : <<< $(date) " >> "$LOG_FILE"
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P10 bash -c 'selectPASS  "$@"' 'selectPASS'

mkdir -p $OUTPUT_DIR/FILTER_VAR 
echo ">>>>>> selectVAR : <<< $(date) " >> "$LOG_FILE"
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P10 bash -c 'selectVar  "$@"' 'selectVar'

mkdir -p $OUTPUT_DIR/annotation/
xargs -a ${SAMPLE_LIST_BAM} -t -n1 -P${JOBS} bash -c 'annotation "$@"' 'annotation'

# USEI O LOOP FOR ABAIXO PARA ANNOTARA AS VARIANTES
for i in `cat ${SAMPLE_LIST_BAM}`; do
 annotation $i
done


# PASSO  - Concatenar tabelas do annovar

SAMPLES=$(find $OUTPUT_DIR/annotation -maxdepth 1 -mindepth 1  -name '*.hg38_multianno.txt')
for SAMPLE in $SAMPLES; do
  id="${SAMPLE##*/}"
  awk -OFS="\t" -v N=${id%.hg38*} '{print N,$_ }' $SAMPLE >> $OUTPUT_DIR/Final_GATK.4.6_annotated.txt
done

sed -i 's/\s/\t/' $OUTPUT_DIR/Final_GATK.4.6_annotated.txt

echo "" >> "$LOG_FILE"
echo "${TIME} >>>>>> End Pipeline <<< " >> "$LOG_FILE"
date >> "$LOG_FILE"
echo "" >> "$LOG_FILE"


mkdir $OUTPUT_DIR/vcfs_toMutalisk
SAMPLES=$(find $OUTPUT_DIR/mutations_to_filter_fromVCF/synonymous -maxdepth 1 -mindepth 1  -name '*.synonymous.bed')

for SAMPLE in ${SAMPLES}; do
  id="${SAMPLE##*/}"
  vcftools --gzvcf $OUTPUT_DIR/FILTER_VAR/"${id%.synonymous.bed*}.pass.var.vcf.gz" \
   --positions ${SAMPLE} \
   --out $OUTPUT_DIR/vcfs_toMutalisk/"${id%.synonymous.bed*}" --recode --keep-INFO-all \
   2> $OUTPUT_DIR/vcfs_toMutalisk/"${id%.synonymous.bed*}.log2"
done

   bedtools intersect -wa \
    -a $OUTPUT_DIR/FILTER_VAR/"${id%.synonymous.bed*}.pass.var.vcf.gz" \
    -b  ${SAMPLE} \
    > $OUTPUT_DIR/Mutation_toMutalisk/$id.bedtools
    
/home/scratch90/vlira_13may2024//Result_Mutect2.GATK4.6.2024-06-14/FILTER_VAR/ROP-97-ExC85-xgenV2_S65.pass.var.vcf.gz


http://mutalisk.org/result.php?rid=Snc2VuzCmJ

#synonymous
http://mutalisk.org/result.php?rid=YJqifIJ7TW

# PCAWG - SigProfiler 
http://mutalisk.org/result.php?rid=dn13TvlmTd