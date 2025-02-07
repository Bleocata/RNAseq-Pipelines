#################1_step: filtrar arquivo bam, mantendo o cabecalho ####################

filter_bam () {
    local sample="$1"
    local id
    id=$(basename "$sample" .rna_seq.genomic.gdc_realn.bam)

    echo ">>>>> Iniciando Filtragem para Amostra: $id <<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

    samtools view -h "${pre_processing_DIR}/BQSR_base_recalibration/${id}_marked_duplicates.tags.recal.bam" | \
    awk '$1 ~ /^@/ || $3 ~ /chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/' | \
    samtools view -b -o "${pre_processing_DIR}/Filtered/${id}_filtered_marked_duplicates.tags.recal.bam" \
    2> "${pre_processing_DIR}/Filtered/${id}_filtered_marked_duplicates.tags.bam.log"
}
export -f filter_bam

mkdir -p "${pre_processing_DIR}/Filtered"
find "$BAM_LIST" -name "*.bam" | parallel -j 4 filter_bam {}

#Verificar correspondecnia das leituras mapeadas presentes no bam com o cabecalho: 
samtools idxstats <arquivo.bam>


#############2_step: reescrever o cabecalho ##############
for BAM_FILE in *_filtered_marked_duplicates.tags.recal.bam; do
    # Extrair e filtrar o cabeçalho, mantendo chrM
    samtools view -H "$BAM_FILE" | \
    awk '$0 ~ /^@SQ/ && $2 ~ /SN:chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ || $0 !~ /^@SQ/' > new_header.sam

    # Substituir o cabeçalho antigo pelo novo e criar o arquivo BAM atualizado
    {
        cat new_header.sam
        samtools view "$BAM_FILE"
    } | samtools view -b -o "${BAM_FILE%.bam}.new.bam"

    # Remover o arquivo temporário do cabeçalho
    rm new_header.sam

    # Mensagem de sucesso
    echo "Updated header for: $BAM_FILE -> ${BAM_FILE%.bam}.new.bam"
done

#Verificar correspondecnia das leituras mapeadas presentes no bam com o cabecalho: 
samtools idxstats <arquivo.bam>

################3_step: indexar os arquivos bam filtrados com o cabecalho atualizado####################

Index_bam() {
  local bam_file="$1"  # Caminho completo do arquivo BAM
  local base_name
  base_name=$(basename "$bam_file" .bam)  # Extrai o nome base do arquivo

  echo ">>>>>> Iniciando indexação do arquivo BAM: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

  # Executa o comando samtools index
  if samtools index "$bam_file"; then
    echo ">>>>>> Indexação concluída com sucesso para: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
  else
    echo ">>>>>> Erro: Falha ao indexar $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
    return 1
  fi
}
export -f Index_bam

# Localiza arquivos BAM e executa a função em paralelo
find . -name "*_filtered_marked_duplicates.tags.new.bam" | parallel -j 4 Index_bam {}


update_header() {
  local bam_file="$1"  # Caminho do arquivo BAM recebido
  local id=$(basename "$bam_file" _filtered_marked_duplicates.tags.recal.bam)  # Extrai o nome base do arquivo sem a extensão .bam

  echo ">>>>>> Atualizando cabeçalho para: $bam_file <<<<<< $(date)" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"

  # Extrair e filtrar o cabeçalho, mantendo chrM
  samtools view -H "$bam_file" | \
  awk '$0 ~ /^@SQ/ && $2 ~ /SN:chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ || $0 !~ /^@SQ/' > "${pre_processing_DIR}/Filtered/${id}_header.sam"

  # Substituir o cabeçalho antigo pelo novo e criar o arquivo BAM atualizado
  {
      cat "${pre_processing_DIR}/Filtered/${id}_header.sam"
      samtools view "$bam_file"
  } | samtools view -b -o "${bam_file%.bam}.new.bam"

  # Remover o arquivo temporário de cabeçalho
  rm "${id}_header.sam"

  echo "Updated header for: $bam_file -> ${bam_file%.bam}.new.bam" >> "${pre_processing_DIR}/preprocessing_$(date +%Y%m%d).log"
}
export -f update_header

# Usar find e parallel para processar os arquivos BAM
find "${pre_processing_DIR}/Filtered" -name "*_filtered_marked_duplicates.tags.recal.bam" | parallel -j 4 update_header {}
