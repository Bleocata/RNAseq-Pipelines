#!/usr/bin/bash

#!/usr/bin/bash

# Pipeline para Análise de RNAseq em Pacientes com Mieloma Múltiplo
# Requerimentos: SRAToolkit,FastQC, MultiQC, Trimmomatic,STAR alignment,samtools
# Recomendação: Rodar este script em uma 'screen' para persistência da execução
# screen -S  TESTE_SCRIPT # para criar uma screen e nomeá-la

# Definindo Variáveis Globais
export WD="/home/scratch90/bleocata_20240912/TESTE_SCRIPT" # Diretório de trabalho
export MANIFEST="${WD}/manifest_list.txt" # Caminho do arquivo de manifest
export LOG_FILE="${WD}/Pipeline_RNAseq_$(date +%Y%m%d).log" #Monitorar e registrar as atividades e o status de execução do pipeline.
export ADAPTERS_FILE="/home/scratch90/bleocata_20240912/PRJEB37100/Adapters_sequences.fa" #Arquivo de adaptadores para o Trimmomatic
export THREADS=5
export JOBS=4
export ref_GRCh38="${WD}/STAR_alignment/Reference"
export BAM_DIR="${WD}/STAR_alignment/runSTAR"


# Variáveis para ferramentas 
export TRIMMOMATIC_JAR="/home/tools/manual/Trimmomatic-0.36/trimmomatic-0.36.jar"
export STAR="/home/tools/bin/STAR"

# Criação de Diretórios Necessários
mkdir -p "${WD}/sra/" "${WD}/fastq/" "${WD}/result_fastqc/" "${WD}/result_multiqc/" "${WD}/Results_trimming_MM/" "${WD}/result_fastqc_post_trimming/" "${WD}/result_multiqc_post_trimming/" "ref_GRCh38" "${WD}/STAR_alignment/genomeSTAR/" "${WD}/stats_reports2/" "${WD}/coverage_reports2/"

echo ">>> Starting Pipeline_RNAseq.sh <<<<<< $(date)" > "$LOG_FILE"

######## PRE-PROCESSAMENTO ##########
# Função para pré-processamento: download SRA, conversão Fastq e controle de qualidade inicial FastQC e MultiQC

step1_preProcessing (){
  local SAMPLE=$1 #define uma variável local chamada SAMPLE que armazena a amostra recebida como argumento quando a função é chamada.
  echo ">>>> Download SRA: ${SAMPLE} <<< $(date)" > "$LOG_FILE"  
  prefetch "${SAMPLE}" -O "${WD}/sra/"
  
  #converte SRA em FASTQ/divide os arquivos em F e R/comprime
  echo ">>>>>> Convert SRA to FASTQ: ${SAMPLE} <<< $(date)" > "$LOG_FILE"
  fastq-dump --split-files --gzip -O "${WD}/fastq/" "${WD}/sra/${SAMPLE}/${SAMPLE}.sra"
  
  echo ">>>>>> Quality Control FASTQC: ${SAMPLE} <<< $(date)" > "$LOG_FILE"
  fastqc -O "${WD}/result_fastqc/" "${WD}/fastq/${SAMPLE}_1.fastq.gz" "${WD}/fastq/${SAMPLE}_2.fastq.gz"
  
  # Exclui o arquivo SRA após conversão para economizar espaço
  rm "${WD}/sra/${SAMPLE}/${SAMPLE}.sra"

}

xargs -a "${MANIFEST}" -n1 -P"${JOBS}" bash -c 'step1_preProcessing "$@"' _

# MultiQC para análise agregada
echo ">>>>>> MultiQC <<< $(date)" > "$LOG_FILE"
multiqc "${WD}/result_fastqc/" -o "${WD}/result_multiqc/"

########### TRIMMOMATIC #############
# Etapa de Trimming com Trimmomatic para cada par de arquivos FASTQ

# Função para processar arquivos com o Trimmomatic
trimmomatic() {
  local forward_read="$1"
  local reverse_read="${forward_read/_1.fastq.gz/_2.fastq.gz}"
  
  # Verifica se o arquivo reverse correspondente existe
  if [ ! -f "$reverse_read" ]; then
  echo "Arquivo reverse não encontrado para: $forward_read" >> "$LOG_FILE"
  return 1
  fi
  
  #Extrai o nome base do arquivo para nomeação dos arquivos de saída
  local base_name=$(basename "$forward_read" "_1.fastq.gz")
  # Define os caminhos de saída para os arquivos paired e unpaired
  local output_forward_paired="${WD}/Results_trimming_MM/${base_name}_1_paired.fastq.gz"
  local output_forward_unpaired="${WD}/Results_trimming_MM/${base_name}_1_unpaired.fastq.gz"
  local output_reverse_paired="${WD}/Results_trimming_MM/${base_name}_2_paired.fastq.gz"
  local output_reverse_unpaired="${WD}/Results_trimming_MM/${base_name}_2_unpaired.fastq.gz"
  
  echo ">>>>>> Início do Trimming para: ${base_name} <<< $(date)" >> "$LOG_FILE"
  # Executa o Trimmomatic com os parâmetros específicos
  java -jar "$TRIMMOMATIC_JAR" PE -threads "$THREADS" \
  "$forward_read" "$reverse_read" \
  "$output_forward_paired" "$output_forward_unpaired" \
  "$output_reverse_paired" "$output_reverse_unpaired" \
  ILLUMINACLIP:"$ADAPTERS_FILE":2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  #obs:os arquivos fasta com as sequencias de adaptadores estao em home/tools/manual/Trimmomatic/adapters/
  
  if [ $? -eq 0 ]; then
  echo "Trimming finalizado com sucesso para: $base_name" >> "$LOG_FILE"
  else
    echo "Erro no trimming para: $base_name" >> "$LOG_FILE"
  fi
}
# Itera sobre os arquivos e chama a função trimmomatic usando xargs
find "${WD}/fastq/" -name "*_1.fastq.gz" | xargs -I {} -P "${JOBS}" bash -c 'trimmomatic "$@"' _ {}

######### Controle de qualidade Pos-trimagem ########
# FastQC pós-trimming e MultiQC
echo ">>>>>> Running FastQC on trimmed files <<< $(date)" > "$LOG_FILE"
find "${WD}/Results_trimming_MM" -name "*_paired.fastq.gz" |xargs -n 1 -P 10 -I {} fastqc -O "${WD}/result_fastqc_post_trimming/" {}

echo ">>>>>> Running MultiQC on FastQC results <<< $(date)" > "$LOG_FILE"
multiqc "${WD}/result_fastqc_post_trimming/" -o "${WD}/result_multiqc_post_trimming/"


########## ALINHAMENTO E MAPEAMENTO ###########
#1_Construir um genoma de referencia
buildIndexGenome (){
  wget -P "ref_GRCh38/" https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz 
  wget -P "ref_GRCh38/" https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz
  
  if ls "ref_GRCh38/"*.gz 1> /dev/null 2>&1; then
  gunzip "ref_GRCh38/"*.gz
  fi
  
  ${STAR} \
  --runThreadN ${JOBS} \
  --runMode genomeGenerate \
  --genomeDir "${WD}/STAR_alignment/genomeSTAR/" \
  --genomeFastaFiles "ref_GRCh38/GRCh38.primary_assembly.genome.fa" \
  --sjdbGTFfile "ref_GRCh38/gencode.v44.primary_assembly.annotation.gtf" \
  --sjdbOverhang 74 \
  > "${WD}/STAR_alignment/genomeSTAR/buildIndexGenome.out" 2> "${WD}/STAR_alignment/genomeSTAR/buildIndexGenome.log"
  
  # Explicação dos parâmetros:
  # --runThreadN: Número de threads para paralelizar o processo.
  # --runMode genomeGenerate: Define o modo para gerar o índice do genoma.
  # --genomeDir: Diretório onde o índice do genoma será salvo.
  # --genomeFastaFiles: Arquivo(s) FASTA contendo a sequência do genoma.
  # --sjdbGTFfile: Arquivo GTF de anotação do genoma (com informações de exon).
  # --sjdbOverhang: Tamanho da leitura - 1. Para leituras de 101 pares de bases, use 100.
  
  if [ $? -ne 0 ]; then
  echo "Erro durante a construção do índice do genoma." >> "$LOG_FILE"
  exit 1
  else
    echo "Índice do genoma criado com sucesso!" >> "$LOG_FILE"
  fi
}

echo ">>>>>> Construindo o índice do genoma <<< $(date)" >> "$LOG_FILE"
buildIndexGenome

#2_Mapeamento das reads com o genoma de referencia
# Gera a lista de IDs de amostras a partir dos arquivos .fastq.gz em fastq_MM
find "${WD}/Results_trimming_MM/" -maxdepth 1 -mindepth 1 -name '*_paired.fastq.gz' | \
xargs -I {} basename {} | \
cut -d'_' -f1 | \
sort -u > "${WD}/samples_ID.list"

runSTAR(){
  echo ">>>>>> Executando o mapeamento <<< $(date)" > "$LOG_FILE" 
  local ID=$1
  
  ${STAR} \
  --runThreadN ${JOBS} \
  --runMode alignReads \
  --genomeDir "${WD}/STAR_alignment/genomeSTAR" \
  --quantMode TranscriptomeSAM GeneCounts \
  --outSAMtype BAM SortedByCoordinate \
  --readFilesCommand zcat \
  --readFilesIn "${WD}/Results_trimming_MM/${ID}_1_paired.fastq.gz" "${WD}/Results_trimming_MM/${ID}_2_paired.fastq.gz" \
  --outFileNamePrefix "$BAM_DIR/${ID}." \
  > "$BAM_DIR/${ID}.mapped.log"
  
  ####PARA O STAR FUSION ###
  # --runThreadN ${JOBS} \
  # --runMode alignReads \
  # --genomeDir "${WD}/STAR_alignment/genomeSTAR" \
  #--quantMode TranscriptomeSAM GeneCounts \
  #--readFilesIn "${WD}/Results_trimming_MM/${ID}_1_paired.fastq.gz" "${WD}/Results_trimming_MM/${ID}_2_paired.fastq.gz" \
  #--readFilesCommand zcat \
  #--chimSegmentMin 15 \
  #--chimOutType Junctions \
  #--chimJunctionOverhangMin 20 \
  #--chimScoreMin 0 \
  #--chimScoreDropMax 20 \
  #--chimFilter banGenomicN \
  #--chimMainSegmentMultNmax 1 \
  #--twopassMode Basic \
  #--outSAMtype BAM SortedByCoordinate \
  #--outFileNamePrefix "$BAM_DIR/${ID}." \
  #> "$BAM_DIR/${ID}.mapped.log" 2>&1
  
  # Número de threads paralelas para acelerar o processo
  # Modo de execução para realizar o alinhamento
  # Diretório do genoma indexado
  # Gera contagens de genes e um SAM de transcritoma para análise de expressão
  # Arquivos FASTQ de entrada
  # Comando para descompactar arquivos .gz
  # Segmento mínimo de 15 pb para detecção de fusão
  # Gera o arquivo Chimeric.out.junction para STAR-Fusion
  # Tamanho mínimo de sobreposição para junções quiméricas
  # Pontuação mínima para considerar fusões (0 para máxima sensibilidade)
  # Máximo diferencial de pontuação para fusões
  # Filtra fusões que têm Ns na sequência genômica
  # Desativa multimapeamentos no segmento principal
  # Habilita o modo de duas passagens para melhor detecção de junções
  # Saída BAM ordenada por coordenadas
  # Prefixo para os arquivos de saída, incluindo logs e BAM
  # Redireciona a saída de log e erro para o mesmo arquivo
}

xargs -a "${WD}/samples_ID.list" -t -n1 -P${JOBS} bash -c 'runSTAR "$@"' _ 

############# Processamento dos arquivos BAM #############

process_bam() {
  echo ">>>>>> Criando uma lista de ID com os arquivos BAMs <<<<<< $(date)" > "$LOG_FILE"
  local bam_file="$1"
  local base_name
  base_name=$(basename "$bam_file" .Aligned.sortedByCoord.out.bam)
  local coverage_report="${WD}/coverage_reports2/${base_name}_coverage2.txt"
  local stats_report="${WD}/stats_reports2/${base_name}_stats2.txt"
  
  # Criar o índice do arquivo BAM ordenado
  echo ">>>>>> samtools index <<<<<< $(date)" > $LOG_FILE
  if ! samtools index "$bam_file"; then
  echo "Erro: falha ao indexar $bam_file" >&2
  return 1
  fi
  
  # Gerar o relatório de cobertura
  echo ">>>>>> samtools coverage <<<<<< $(date)" > "$LOG_FILE"
  if ! samtools coverage "$bam_file" > "$coverage_report"; then
  echo "Erro: falha ao gerar o relatório de cobertura para $bam_file" >&2
  return 1
  fi
  
  # Gerar o relatório de estatísticas
  echo ">>>>>> samtools stats <<<<<< $(date)" > "$LOG_FILE"
  if ! samtools stats "$bam_file" > "$stats_report"; then
  echo "Erro: falha ao gerar o relatório de estatísticas para $bam_file" >&2
  return 1
  fi
  
  echo "Processamento completo para: $base_name"
}

find "$BAM_DIR" -name "*.Aligned.sortedByCoord.out.bam" | xargs -I {} -P "${JOBS}" bash -c 'process_bam "$@"' _ 


echo ">>>>>> Fim do Pipeline_RNAseq.sh <<<<<< $(date)" > "$LOG_FILE"

############ Explorando os dados de RNAseq ##############
