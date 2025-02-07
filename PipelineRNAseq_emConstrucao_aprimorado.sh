#!/usr/bin/bash

# Pipeline para Análise de RNAseq em Pacientes com Mieloma Múltiplo
# Requerimentos: SRAToolkit,FastQC, MultiQC, Trimmomatic,STAR alignment,samtools
# Recomendação: Rodar este script em uma 'screen' para persistência da execução
# screen -S PRJEB37100 # para criar uma screen e nomeá-la

# Definindo Variáveis Globais
export wd="/home/scratch90/bleocata_20241205/PRJEB37100" # Diretório de trabalho
export MANIFEST="${wd}/manifest_list.txt" # Caminho do arquivo de manifest
export LOG_FILE="${wd}/Pipeline_RNAseq_$(date +%Y%m%d).log" #Monitorar e registrar as atividades e o status de execução do pipeline.
export ADAPTERS_FILE="${wd}/Adapters_sequences.fa" #CUIDADO: Arquivo de adaptadores para o Trimmomatic
export THREADS=15
export JOBS=4
export ref_GRCh38="${wd}/STAR_alignment/Reference"



# Variáveis para ferramentas 
export TRIMMOMATIC_JAR="/home/tools/manual/Trimmomatic-0.36/trimmomatic-0.36.jar"
export STAR="/home/tools/bin/STAR"

# Criação de Diretórios Necessários
mkdir -p "${wd}/sra/" "${wd}/fastq/" "${wd}/result_fastqc/" "${wd}/result_multiqc/" "${wd}/Results_trimming_MM/" "${wd}/result_fastqc_post_trimming/" "${wd}/result_multiqc_post_trimming/" "${ref_GRCh38}" "${wd}/STAR_alignment/genomeSTAR/" "${BAM_DIR}" "${wd}/stats_reports2/" "${wd}/coverage_reports2/"

echo ">>> Starting Pipeline_RNAseq.sh <<<<<< $(date)" > "$LOG_FILE"

######## PRE-PROCESSAMENTO ##########
# Função para pré-processamento: download SRA, conversão Fastq e controle de qualidade inicial FastQC e MultiQC

step1_preProcessing(){
  local SAMPLE=$1 #define uma variável local chamada SAMPLE que armazena a amostra recebida como argumento quando a função é chamada.
  echo ">>>> Download SRA: ${SAMPLE} <<< $(date)" >> "$LOG_FILE"  
  prefetch "${SAMPLE}" -O "${wd}/sra/"
  
  #converte SRA em FASTQ/divide os arquivos em F e R/comprime
  echo ">>>>>> Convert SRA to FASTQ: ${SAMPLE} <<< $(date)" >> "$LOG_FILE"
  fastq-dump --split-files --gzip -O "${wd}/fastq/" "${wd}/sra/${SAMPLE}/${SAMPLE}.sra"
  
  echo ">>>>>> Quality Control FASTQC: ${SAMPLE} <<< $(date)" >> "$LOG_FILE"
  fastqc -O "${wd}/result_fastqc/" "${wd}/fastq/${SAMPLE}_1.fastq.gz" "${wd}/fastq/${SAMPLE}_2.fastq.gz"
  
  # Exclui o arquivo SRA após conversão para economizar espaço
  rm "${wd}/sra/${SAMPLE}/${SAMPLE}.sra"

}
export -f step1_preProcessing
xargs -a "${MANIFEST}" -n1 -P"${JOBS}" bash -c 'step1_preProcessing "$@"' _

# MultiQC para análise agregada
echo ">>>>>> MultiQC <<< $(date)" >> "$LOG_FILE"
multiqc "${wd}/result_fastqc/" -o "${wd}/result_multiqc/"

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
  local output_forward_paired="${wd}/Results_trimming_MM/${base_name}_1_paired.fastq.gz"
  local output_forward_unpaired="${wd}/Results_trimming_MM/${base_name}_1_unpaired.fastq.gz"
  local output_reverse_paired="${wd}/Results_trimming_MM/${base_name}_2_paired.fastq.gz"
  local output_reverse_unpaired="${wd}/Results_trimming_MM/${base_name}_2_unpaired.fastq.gz"
  
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
export -f trimmomatic
# Itera sobre os arquivos e chama a função trimmomatic usando xargs
find "${wd}/fastq/" -name "*_1.fastq.gz" | xargs -I {} -P "${JOBS}" bash -c 'trimmomatic "$@"' _ {}

######### Controle de qualidade Pos-trimagem ########
# FastQC pós-trimming e MultiQC
echo ">>>>>> Running FastQC on trimmed files <<< $(date)" >> "$LOG_FILE"
find "${wd}/Results_trimming_MM" -name "*_paired.fastq.gz" |xargs -n 1 -P 10 -I {} fastqc -O "${wd}/result_fastqc_post_trimming/" {}

echo ">>>>>> Running MultiQC on FastQC results <<< $(date)" >> "$LOG_FILE"
multiqc "${wd}/result_fastqc_post_trimming/" -o "${wd}/result_multiqc_post_trimming/"


########## ALINHAMENTO E MAPEAMENTO ###########
#1_Construir um genoma de referencia
buildIndexGenome(){
  wget -P "${ref_GRCh38}/" https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz 
  wget -P "${ref_GRCh38}/" https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz
  
  gunzip "${ref_GRCh38}/"*.gz
  
  ${STAR} \
  --runThreadN ${THREADS} \
  --runMode genomeGenerate \
  --genomeDir "${wd}/STAR_alignment/genomeSTAR/" \
  --genomeFastaFiles "${ref_GRCh38}/GRCh38.primary_assembly.genome.fa" \
  --sjdbGTFfile "${ref_GRCh38}/gencode.v44.primary_assembly.annotation.gtf" \
  --sjdbOverhang 74 \
  > "${wd}/STAR_alignment/genomeSTAR/buildIndexGenome.out" 2> "${wd}/STAR_alignment/genomeSTAR/buildIndexGenome.log"
  
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
export -f buildIndexGenome
echo ">>>>>> Construindo o índice do genoma <<< $(date)" >> "$LOG_FILE"
buildIndexGenome

#2_Mapeamento das reads com o genoma de referencia
# Gera a lista de IDs de amostras a partir dos arquivos .fastq.gz em fastq_MM
find "${wd}/Results_trimming_MM/" -maxdepth 1 -mindepth 1 -name '*_paired.fastq.gz' | \
xargs -I {} basename {} | \
cut -d'_' -f1 | \
sort -u > "${wd}/samples_ID.list"

runSTAR(){
  echo ">>>>>> Executando o mapeamento <<< $(date)" >> "$LOG_FILE" 
  local ID=$1
  
  ${STAR} \
  --runThreadN ${THREADS} \
  --runMode alignReads \
  --genomeDir "${wd}/STAR_alignment/genomeSTAR" \
  --quantMode TranscriptomeSAM GeneCounts \
  --outSAMtype BAM SortedByCoordinate \
  --readFilesCommand zcat \
  --readFilesIn "${wd}/Results_trimming_MM/${ID}_1_paired.fastq.gz" "${wd}/Results_trimming_MM/${ID}_2_paired.fastq.gz" \
  --outFileNamePrefix "${BAM_DIR}/${ID}." > "${BAM_DIR}/${ID}.mapped.log"
  
  ####PARA O STAR FUSION ###
  # --runThreadN ${JOBS} \
  # --runMode alignReads \
  # --genomeDir "${wd}/STAR_alignment/genomeSTAR" \
  #--quantMode TranscriptomeSAM GeneCounts \
  #--readFilesIn "${wd}/Results_trimming_MM/${ID}_1_paired.fastq.gz" "${wd}/Results_trimming_MM/${ID}_2_paired.fastq.gz" \
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
export -f runSTAR
xargs -a "${wd}/samples_ID.list" -t -n1 -P${JOBS} bash -c 'runSTAR "$@"' _ 

############# Processamento dos arquivos BAM #############


mkdir -p "${wd}/coverage_reports2" "${wd}/stats_reports2"

process_bam() {
  echo ">>>>>> Criando uma lista de ID com os arquivos BAMs <<<<<< $(date)" >> "$LOG_FILE"
  local bam_file="$1"
  local base_name
  base_name=$(basename "$bam_file" .Aligned.sortedByCoord.out.bam)
  local coverage_report="${wd}/coverage_reports2/${base_name}_coverage2.txt"
  local stats_report="${wd}/stats_reports2/${base_name}_stats2.txt"
  
  # Criar o índice do arquivo BAM ordenado
  echo ">>>>>> samtools index <<<<<< $(date)" >> $LOG_FILE
  if ! samtools index "$bam_file"; then
  echo "Erro: falha ao indexar $bam_file" >&2
  return 1
  fi
  
  # Gerar o relatório de cobertura
  echo ">>>>>> samtools coverage <<<<<< $(date)" >> "$LOG_FILE"
  if ! samtools coverage "$bam_file" > "$coverage_report"; then
  echo "Erro: falha ao gerar o relatório de cobertura para $bam_file" >&2
  return 1
  fi
  
  # Gerar o relatório de estatísticas
  echo ">>>>>> samtools stats <<<<<< $(date)" >> "$LOG_FILE"
  if ! samtools stats "$bam_file" > "$stats_report"; then
  echo "Erro: falha ao gerar o relatório de estatísticas para $bam_file" >&2
  return 1
  fi
  
  echo "Processamento completo para: $base_name"
}
export -f process_bam
find "$BAM_DIR" -name "*.Aligned.sortedByCoord.out.bam" | \
parallel -j "${JOBS}" process_bam {}


######### Rodar o Multiqc para visualizacao do relatorio completo ##########

# Criar diretórios de saída se ainda não existirem
mkdir -p "${wd}/multiqc_coverage2/"
mkdir -p "${wd}/multiqc_stats2/"


multiqc_run() {
    local input_dir="$1"
    local output_dir="$2"

    # Verifica se os diretórios foram fornecidos
    if [[ -z "$input_dir" || -z "$output_dir" ]]; then
        printf "Erro: Diretórios de entrada ou saída não especificados.\n" >&2
        return 1
    fi

    # Verifica se o diretório de entrada existe
    if [[ ! -d "$input_dir" ]]; then
        printf "Erro: Diretório de entrada '%s' não encontrado.\n" "$input_dir" >&2
        return 1
    fi

    # Cria o diretório de saída se não existir
    mkdir -p "$output_dir"

    # Executa o MultiQC e captura possíveis errosls
    if ! multiqc "$input_dir" -o "$output_dir"; then
        printf "Erro ao executar MultiQC no diretório '%s'.\n" "$input_dir" >&2
        return 1
    fi
}

export -f multiqc_run

# Diretórios de entrada e saída para MultiQC
inputs=("${wd}/coverage_reports2/" "${wd}/stats_reports2/")
outputs=("${wd}/multiqc_coverage2/" "${wd}/multiqc_stats2/")

# Executar MultiQC em paralelo
parallel -j "${JOBS}" multiqc_run {1} {2} ::: "${inputs[@]}" ::: "${outputs[@]}"


############# Cicero: Chamada de variantes estruturais - Fusoes  ##################

export wd="/home/scratch90/bleocata_20241205/PRJEB37100"
export LOG_FILE="${wd}/Pipeline_RNAseq_$(date +%Y%m%d).log"
export JOBS=4
export bam_dir="${wd}/STAR_alignment/runSTAR"

echo ">>>>>> Baixando os arquivos de referência <<<" >> "$LOG_FILE"

# Baixa o arquivo de referência no diretório especificado
wget -P "${wd}" -O "${wd}/reference.tar.gz" "https://zenodo.org/records/5088371/files/reference.tar.gz?download=1"

# Extrai o arquivo tar.gz no diretório de trabalho
tar -xvzf "${wd}/reference.tar.gz" -C "${wd}"
mkdir -p "${wd}/output_RNApeg"

run_RNApeg() {
    local bam_file="$1"
    local id
    id=$(basename "$bam_file" .Aligned.sortedByCoord.out.bam)

    # Verifica se os arquivos de referência existem
    if [[ ! -f "${wd}/reference/Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa" ]] || \
       [[ ! -f "${wd}/reference/Homo_sapiens/GRCh38_no_alt/mRNA/RefSeq/refFlat.txt" ]]; then
        echo "Arquivos de referência faltando para a amostra $id" >> "$LOG_FILE"
        return 1
    fi

    echo ">>>>>> Executando RNApeg para amostra: $id <<<" >> "$LOG_FILE"

    # Executa o container Docker sem comentários in-line
    docker run --rm --privileged \
        -u "$(id -u):$(id -g)" \
        -v "${wd}/reference:/reference" \
        -v "${wd}/STAR_alignment:/data" \
        -v "${wd}/output_RNApeg:/output_RNApeg" \
        ghcr.io/stjude/rnapeg:latest \
        -b "/data/runSTAR/${id}.Aligned.sortedByCoord.out.bam" \
        -f "/reference/Homo_sapiens/GRCh38_no_alt/FASTA/GRCh38_no_alt.fa" \
        -r "/reference/Homo_sapiens/GRCh38_no_alt/mRNA/RefSeq/refFlat.txt" \
        -O "/output_RNApeg/output_${id}" \
        > "${wd}/output_RNApeg/RNApeg.${id}.log" \
        2> "${wd}/output_RNApeg/RNApeg.${id}.log2"

    # Verifica a saída do Docker
    if [[ $? -ne 0 ]]; then
        echo "Erro ao executar RNApeg para a amostra $id" >> "$LOG_FILE"
    else
        echo "RNApeg executado com sucesso para a amostra $id" >> "$LOG_FILE"
    fi
}

export -f run_RNApeg

# Executa a função RNApeg em paralelo para múltiplas amostras
find "${bam_dir}" -name "*.Aligned.sortedByCoord.out.bam" | parallel -j "${JOBS}" run_RNApeg {}



########## Rodando o Cicero ################

mkdir -p "${wd}/output_Cicero/"

run_Cicero (){
    local bam_files="$1"
    local id
    id=$(basename "$bam_files" .Aligned.sortedByCoord.out.bam)

    echo ">>>>>> Executando Cicero para amostra: $id <<<" >> "$LOG_FILE"
    date >> "$LOG_FILE"

    # Executa o container Docker para o Cicero
    docker run --rm --privileged \
        -u "$(id -u):$(id -g)" \
        -v "${wd}/reference:/reference" \
        -v "${wd}/STAR_alignment:/data" \
        -v "${wd}/output_Cicero:/output_Cicero" \
        -v "${wd}/output_RNApeg:/output_RNApeg" \
        ghcr.io/stjude/cicero:latest Cicero.sh \
        -b="/data/runSTAR/${id}.Aligned.sortedByCoord.out.bam" \
        -g="GRCh38_no_alt" \
        -r="/reference"  \
        -o="/output_Cicero/output_Cicero_${id}" \
        -j="/output_RNApeg/output_${id}/${id}.Aligned.sortedByCoord.out.bam.junctions.tab.shifted.tab"

    # Verifica o sucesso da execução
    if [[ $? -ne 0 ]]; then
        echo "Erro ao executar Cicero para a amostra $id" >> "$LOG_FILE"
    else
        echo "Cicero executado com sucesso para a amostra $id" >> "$LOG_FILE"
    fi
}

export -f run_Cicero

# Executa Cicero em paralelo para múltiplos BAMs
find "${bam_dir}" -name "*.Aligned.sortedByCoord.out.bam" | parallel -j "${JOBS}" run_Cicero {}


# CRIA UM LINK SIMBOLICO DE TODOS OS ARQUIVOS final_fusions.txt PARA O DIRETORIO output_Cicero_ALLsamples > Facilita visualizacao 

mkdir -p "${wd}/output_Cicero_ALLsamples/"

Cicero_Allsamples (){
    local id=$1
    echo "" >> $LOG_FILE
    echo ">>>>>> Executando Cicero_Allsamples para amostra: $ID <<<" >> $LOG_FILE
    date >> $LOG_FILE

    ln -s ${wd}/output_Cicero/output_Cicero_${id}/CICERO_DATADIR/${id}.Aligned.sortedByCoord.out/final_fusions.txt  ${wd}/output_Cicero_ALLsamples/${id}.final_fusions.txt
}
export -f Cicero_Allsamples
find "${bam_dir}" -name "*.Aligned.sortedByCoord.out.bam" | parallel -j "${JOBS}" Cicero_Allsamples {}

