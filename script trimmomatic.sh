    ############ Script Trimmomatic-36 ############
    
#1_Definindo as Variaveis:

# Caminho para o diretório com arquivos FASTQ
WD="/home/scratch90/bleocata_20240912/PRJEB37100"
FASTQ_DIR="$WD/fastq_MM"

# Caminho para o Trimmomatic JAR
TRIMMOMATIC_JAR="/home/tools/manual/Trimmomatic-0.36/trimmomatic-0.36.jar"

# Caminho para a sequencia de adaptadores como TruSeq
ADAPTERS="$WD/Adapters_sequences.fa"

# Diretório de saída
OUTPUT_DIR="/$WD/Results_trimming_MM"

#2_Exportar variaveis 
export WD
export FASTQ_DIR
export TRIMMOMATIC_JAR
export OUTPUT_DIR 
export ADAPTERS

# Itera sobre os arquivos _1.fastq.gz ("forward reads")
for forward_read in "$FASTQ_DIR"/*_1.fastq.gz; do
# Cria o nome do arquivo reverse correspondente (_2.fastq.gz)
reverse_read="${forward_read/_1.fastq.gz/_2.fastq.gz}"

# Verificar se o arquivo reverse existe para evitar erros na execução
if [ ! -f "$reverse_read" ]; then
echo "Arquivo correspondente não encontrado para: $forward_read"
continue
fi #finaliza a estrutura condicional

# Extrai o nome base dos arquivos fastq (sem o _1.fastq.gz) para usar nos arquivos de saída
base_name=$(basename "$forward_read" "_1.fastq.gz")

# Definir as variaveis para os arquivos de saída: para os pares de reads "forward" e "reverse"
output_forward_paired="$OUTPUT_DIR/${base_name}_1_paired.fastq.gz"
output_forward_unpaired="$OUTPUT_DIR/${base_name}_1_unpaired.fastq.gz"
output_reverse_paired="$OUTPUT_DIR/${base_name}_2_paired.fastq.gz"
output_reverse_unpaired="$OUTPUT_DIR/${base_name}_2_unpaired.fastq.gz"

# Executa o Trimmomatic para cada par de arquivos
java -jar "$TRIMMOMATIC_JAR" PE -threads 4 \
"$forward_read" "$reverse_read" \
"$output_forward_paired" "$output_forward_unpaired" \
"$output_reverse_paired" "$output_reverse_unpaired" \
ILLUMINACLIP:"$ADAPTERS":2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#Remove sequencias especificas/ retira bases de baixa qualidade no inicio e final da read/

echo "Trimming finalizado para o par: $base_name"
done

echo "Todos os arquivos foram processados."