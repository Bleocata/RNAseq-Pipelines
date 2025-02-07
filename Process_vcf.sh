#!/bin/bash

# Set the working directory and log directory
export wd="/home/scratch90/bleocata_20241205/MMRF-COMMPASS/Variant_Calling/SM_VCF/Cinco_vcf_gz"
export LOG_DIR="${wd}/Arquivos_log"
export LOG_FILE="${LOG_DIR}/Process_TSV_$(date +%Y%m%d).log"
output_dir="${wd}/Processed_TSVs"

# Create necessary directories
mkdir -p "$LOG_DIR"
mkdir -p "$output_dir"

# Define the function to transform VCF to TSV
Transform_tsv() {
  local sample="$1"
  local id=$(basename "$sample" .wxs.MuTect2.somatic_annotation.vcf.gz)
  
  # Process the VCF file
  bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t[%GT]\t[%AD]\t[%DP]\t[%INFO/CSQ]\n' "${sample}" | \
  Transform_tsv() {
    local sample="$1"       # Caminho do arquivo de entrada (VCF)
    local id=$(basename "$sample" .vcf.gz)  # Extrair o nome base do arquivo sem extensão
    local output_dir="/caminho/para/output" # Substitua pelo caminho do diretório de saída
    mkdir -p "$output_dir"  # Criar o diretório de saída, se não existir

    echo "Processando o arquivo $sample com ID $id"

    # Usar bcftools para extrair informações relevantes do VCF
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER\t[%GT\t%AD\t%DP]\t[%INFO/CSQ]\n' "${sample}" | \
    awk -v N="$id" -F'\t' '
        BEGIN {
            OFS="\t";  # Definir o separador de saída como tabulação
            # Cabeçalho do arquivo de saída
            print "Sample", "CHROM", "POS", "ID", "REF", "ALT", "FILTER", "GT", "AD", "DP", "SYMBOL", "Gene", "Consequence", "IMPACT";
        }
        {
            # Separar os campos do CSQ (assumindo "|" como delimitador)
            split($8, csq_fields, "|");
            
            # Verificar se há informações suficientes no campo CSQ
            if (length(csq_fields) >= 4) {
                SYMBOL=csq_fields[1];       # Nome do gene (ou símbolo)
                Gene=csq_fields[2];         # Gene
                Consequence=csq_fields[3];  # Consequência
                IMPACT=csq_fields[4];       # Impacto
            } else {
                SYMBOL="NA";
                Gene="NA";
                Consequence="NA";
                IMPACT="NA";
            }
            
            # Verificar se há informações genotípicas disponíveis
            if ($7 != "") {
                split($7, gt_fields, "\t");
                GT=gt_fields[1];   # Genótipo
                AD=gt_fields[2];   # Depth dos alelos
                DP=gt_fields[3];   # Profundidade
            } else {
                GT="NA";
                AD="NA";
                DP="NA";
            }
            
            # Imprimir linha formatada
            print N, $1, $2, $3, $4, $5, $6, GT, AD, DP, SYMBOL, Gene, Consequence, IMPACT;
        }
    ' > "${output_dir}/${id}.tsv"

    echo "Arquivo processado: ${output_dir}/${id}.tsv"
}
export -f Transform_tsv


# Find and process files in parallel
find "$wd/" -name "*.wxs.MuTect2.somatic_annotation.vcf.gz" -print0 | \
while IFS= read -r -d '' sample; do
Transform_tsv "$sample"
done


############ concatenar arquivos TSV ######################
# Define o arquivo de saída
output_file="merged_file.tsv"

# Seleciona o cabeçalho do primeiro arquivo encontrado e o salva no arquivo de saída
head -n 1 $(ls *.tsv | head -n 1) > "$output_file"

# Adiciona o conteúdo dos arquivos (sem os cabeçalhos) ao arquivo de saída
for file in *.tsv; do
tail -n +2 "$file" >> "$output_file"
done
