#!/bin/bash
# Script para realizar realinhamento de arquivos BAM utilizando HISAT2.
# Este script inclui:
# 1. Download e preparação do índice do genoma GRCh38.
# 2. Realinhamento de arquivos BAM para um novo índice HISAT2.
# 
# Requisitos:
# - wget
# - gunzip
# - samtools
# - hisat2
# - GNU Parallel
#
# Uso:
# 1. Certifique-se de que as variáveis WD, BAM_LIST e THREADS estejam configuradas corretamente.
# 2. Execute o script diretamente: ./script.sh

# Configurações de variáveis
#!/bin/bash

# Configuração do ambiente
export WD="/home/scratch90/bleocata_20241205/MMRF-COMMPASS/Variant_Calling/TESTE"
export BAM_LIST="/home/scratch90/bleocata_20241205/MMRF-COMMPASS/Variant_Calling/Cinco_samples_BAMs"
export REALIGN_DIR="$WD/HISAT2_alignment"
export HISAT2_INDEX="${REALIGN_DIR}/reference/GRCh38.primary_assembly"
export THREADS=4
export LOG_FILE="${REALIGN_DIR}/Hisat2_realign_$(date +%Y%m%d).log"

# Configurações de segurança
set -euo pipefail

# Verificar pré-requisitos
command -v samtools >/dev/null 2>&1 || { echo "ERRO: samtools não encontrado"; exit 1; }
command -v hisat2 >/dev/null 2>&1 || { echo "ERRO: hisat2 não encontrado"; exit 1; }
command -v parallel >/dev/null 2>&1 || { echo "ERRO: parallel não encontrado"; exit 1; }

# Criar diretórios necessários
mkdir -p "${REALIGN_DIR}/reference" "${REALIGN_DIR}/realigned_reads" "${REALIGN_DIR}/realigned_bams"

# Função: Gerar o índice do genoma HG38
buildIndexGenome() {
  echo "            >>>>>> Iniciando download e preparação do índice do genoma HG38 <<<<<< $(date) " >> "$LOG_FILE"

  wget -P "${REALIGN_DIR}/reference" https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.primary_assembly.genome.fa.gz || {
    echo "            >>>>>> ERRO: Falha ao baixar GRCh38.primary_assembly.genome.fa.gz <<<<<< $(date) " >> "$LOG_FILE"
    exit 1
  }

  wget -P "${REALIGN_DIR}/reference" https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.primary_assembly.annotation.gtf.gz || {
    echo "            >>>>>> ERRO: Falha ao baixar gencode.v36.primary_assembly.annotation.gtf.gz <<<<<< $(date) " >> "$LOG_FILE"
    exit 1
  }

  if [ ! -f "${REALIGN_DIR}/reference/GRCh38.primary_assembly.genome.fa" ]; then
    echo "            >>>>>> Descompactando arquivos do genoma <<<<<< $(date) " >> "$LOG_FILE"
    gunzip "${REALIGN_DIR}/reference/"*.gz
  fi

  if [ ! -f "${REALIGN_DIR}/reference/GRCh38.primary_assembly.1.ht2" ]; then
    echo "            >>>>>> Gerando o índice do genoma com HISAT2 <<<<<< $(date) " >> "$LOG_FILE"
    hisat2-build -p "$THREADS" -f "${REALIGN_DIR}/reference/GRCh38.primary_assembly.genome.fa" "$HISAT2_INDEX"
    echo "            >>>>>> Índice do genoma gerado com sucesso <<<<<< $(date) " >> "$LOG_FILE"
  else
    echo "            >>>>>> Índice do genoma já existente. Pulando etapa <<<<<< $(date) " >> "$LOG_FILE"
  fi
  }

 # Check if the .fai index for the reference genome exists
if [ ! -f "${REALIGN_DIR}/reference/GRCh38.primary_assembly.genome.fa.fai" ]; then
    echo "            >>>>>> Gerando o índice do genoma de referência (.fai) <<<<<< $(date) " >> "$LOG_FILE"
    
    # Generate the .fai index
    if samtools faidx "${REALIGN_DIR}/reference/GRCh38.primary_assembly.genome.fa"; then
        echo "            >>>>>> Índice do genoma de referência (.fai) gerado com sucesso <<<<<< $(date) " >> "$LOG_FILE"
    else
        echo "ERRO: Falha ao gerar o índice do genoma de referência (.fai) <<<<<< $(date) " >> "$LOG_FILE"
        exit 1
    fi
else
    echo "            >>>>>> Índice do genoma de referência (.fai) já existente. Pulando etapa <<<<<< $(date) " >> "$LOG_FILE"
fi

}
export -f buildIndexGenome

buildIndexGenome

# Função: Realinhar arquivos BAM
realign_bam() {
  local bam_file=$1
  local id

  # Extrair o ID do arquivo BAM
  id=$(basename "$bam_file" .rna_seq.genomic.gdc_realn.bam)

  echo "" >> "$LOG_FILE"
  echo "            >>>>>> Iniciando processamento do arquivo BAM: $id <<<<<< $(date) " >> "$LOG_FILE"

  # Etapa 1: Selecionar as reads não alinhadas
  samtools view -b -f 4 "${WD}/TOYS_SAMPLES/${id}.downsampled.bam"> "${WD}/TOYS_SAMPLES/${id}.downsampled.unaligned.bam"

  # Etapa 2: Converter BAM para FASTQ
  if [ ! -f "${WD}/TOYS_SAMPLES/${id}_R1.fastq" ] || [ ! -f "${WD}/TOYS_SAMPLES/${id}_R2.fastq" ]; then
    echo "            >>>>>> Convertendo arquivo BAM para FASTQ: $id <<<<<< $(date) " >> "$LOG_FILE"
    samtools fastq "${WD}/TOYS_SAMPLES/${id}.downsampled.unaligned.bam" \
      -1 "${WD}/TOYS_SAMPLES/${id}_R1.fastq" \
      -2 "${WD}/TOYS_SAMPLES/${id}_R2.fastq" \
      -0 /dev/null -s /dev/null -n || {
      echo "ERRO: Falha ao converter BAM para FASTQ para $id" >> "$LOG_FILE"
      exit 1
    }
  fi

  # Etapa 3: Realinhamento com HISAT2
  echo "" >> "$LOG_FILE"
  echo "            >>>>>> Realinhando com HISAT2: $id <<<<<< $(date) " >> "$LOG_FILE"
  hisat2 -p "$THREADS" --phred33 -x "$HISAT2_INDEX" \
  -1 "${WD}/TOYS_SAMPLES/${id}_R1.fastq" \
  -2 "${WD}/TOYS_SAMPLES/${id}_R2.fastq" \
  -S "${REALIGN_DIR}/realigned_reads/${id}.sam" || {
    echo "ERRO: Falha no realinhamento com HISAT2 para $id" >> "$LOG_FILE"
    exit 1
  }
  echo "            >>>>>> Realinhamento com HISAT2 concluído para: $id <<<<<< $(date) " >> "$LOG_FILE"

  # Etapa 4: Converter SAM para BAM
  echo "            >>>>>> Convertendo SAM para BAM: $id <<<<<< $(date) " >> "$LOG_FILE"
  samtools view -Sb "${REALIGN_DIR}/realigned_reads/${id}.sam" > "${REALIGN_DIR}/realigned_bams/${id}_realigned.bam" || {
    echo "ERRO: Falha ao converter SAM para BAM para $id" >> "$LOG_FILE"
    exit 1
  }

  # Etapa 5: Ordenar e compactar BAM
  echo "            >>>>>> Ordenando BAM: $id <<<<<< $(date) " >> "$LOG_FILE"
  samtools sort -@ "$THREADS" \
    -o "${REALIGN_DIR}/realigned_bams/${id}_realigned.sorted.bam" \
    "${REALIGN_DIR}/realigned_bams/${id}_realigned.bam" 2> "${REALIGN_DIR}/realigned_bams/${id}_realigned.bam.log" || {
    echo "ERRO: Falha ao ordenar BAM para $id" >> "$LOG_FILE"
    exit 1
  }

  # Etapa 6: Indexar BAM
  echo "            >>>>>> Indexando BAM: $id <<<<<< $(date) " >> "$LOG_FILE"
  samtools index "${REALIGN_DIR}/realigned_bams/${id}_realigned.sorted.bam" || {
    echo "ERRO: Falha ao indexar BAM para $id" >> "$LOG_FILE"
    exit 1
  }

  # Limpeza de arquivos temporários
  rm -f "${WD}/TOYS_SAMPLES/${id}.unaligned.bam" \
        "${WD}/TOYS_SAMPLES/${id}_R1.fastq" \
        "${WD}/TOYS_SAMPLES/${id}_R2.fastq" \
        "${REALIGN_DIR}/realigned_reads/${id}.sam" \
        "${REALIGN_DIR}/realigned_bams/${id}_realigned.bam"

  echo "            >>>>>> Processamento concluído com sucesso para: $id <<<<<< $(date) " >> "$LOG_FILE"
}
export -f realign_bam

# Realinhar arquivos BAM em paralelo
echo "            >>>>>> Iniciando processamento de arquivos BAM <<<<<< $(date) " >> "$LOG_FILE"
echo "            >>>>>> Encontrados $(find "$BAM_LIST" -name "*.bam" | wc -l) arquivos BAM para processar <<<<<< $(date) " >> "$LOG_FILE"
find "$BAM_LIST" -name "*.bam" | parallel -j "$THREADS" realign_bam {}

echo "            >>>>>> Pipeline concluído <<<<<< $(date) " >> "$LOG_FILE"