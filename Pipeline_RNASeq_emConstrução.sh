#!/usr/bin/bash

# PERSONALIZADO: CRIANDO VARIAVEIS E ATRIBUINDO VALORES PARA FACILITAR A ESCRITA DO SCRIPT
WD="/home/scratch90/bleocata_20240912/PRJEB37100/" #ATRIBUI O CAMINHO ONDE OS DADOS SRA ESTAO NO SERVIDOR

OUTPUT_DIR=$WD"/IDSamples_MM_RNASeq/" #A VARIÁVEL ARMAZENA O CAMINHO DE DIRETÓRIO ONDE OS DADOS SERÃO SALVOS 

JOBS=5 #NUMERO DE TAREFAS QUE DEVE SER EXECUTADA PARALELAMENTE 

STAR="/home/tools/bin/STAR"

mkdir -p "$WD/reference/" #Criação do diretório "reference"" no caminho de "wd"

OUTPUT_LOG="$OUTPUT_DIR.log" #Armazenar o caminho e o nome de um arquivo de log.

#Exportar as variaveis definidas, tornando visivel e disponivel para utilizacao em scripts no shell > fica visivel para subprocessos 
export wd
export OUTPUT_DIR
export OUTPUT_LOG
export STAR
export JOBS

      ########### BAIXANDO OS DADOS SRA COM SRATOOLKIT  ############
      
# 1_Crei uma pasta para os dados do meu projeto e uma subpasta 
mkdir PRJEB37100
mkdir data

# 2_criei um arquivo nano para colar os IDs so manifest,list baixado
nano manifest_list.txt

# 3_ baixar os arquivos SRA
screen -S DownloadPRJEB37100 #cria uma screen e nomeia 
cat manifest_list.txt| paralell -j 5 "prefetch {} -O data"
 



# OUTRA FORMA DE BAIXAR SRA FILES E TRANSFORMAR EM FASTQ COM SRA TOOLKIT 
	# Baixar as amostras NCBI ou GEO. Precisa conter os nomes dos arquivos separados por espaço " ".
	prefetch GSM6390349	GSM6390350	GSM6390351	GSM6390352	GSM6390353	GSM6390354	GSM6390355	GSM6390356 \
	GSM6390357	GSM6390358	GSM6390359	GSM6390360	GSM6390361	GSM6390362	GSM6390363	GSM6390364	GSM6390365 \
	GSM6390366	GSM6390367	GSM6390368	GSM6390369	GSM6390370	GSM6390371	GSM6390372	GSM6390373	GSM6390374 \
	GSM6390375	GSM6390376	GSM6390377	GSM6390378	GSM6390379	GSM6390380	GSM6390381	GSM6390382	GSM6390383 \
 	GSM6390384	GSM6390385	GSM6390386	GSM6390387	GSM6390388	GSM6390389	GSM6390390	GSM6390391	GSM6390392 \
	GSM6390393	GSM6390394	GSM6390395	GSM6390396	GSM6390397	GSM6390398	GSM6390399	GSM6390400	GSM6390401 \
	GSM6390402	GSM6390403	GSM6390404	GSM6390405	GSM6390406	GSM6390407	GSM6390408	GSM6390409	GSM6390410 \
	GSM6390411	GSM6390412	GSM6390413	GSM6390414	GSM6390415	GSM6390416	GSM6390417	GSM6390418	GSM6390419 \
	GSM6390420	GSM6390421	GSM6390422	GSM6390423	GSM6390424	GSM6390425	GSM6390426	GSM6390427	GSM6390428 \
	GSM6390429	GSM6390430	GSM6390431	GSM6390432	GSM6390433	GSM6390434	GSM6390435	GSM6390436	GSM6390437 \
	GSM6390438	GSM6390439	GSM6390440	GSM6390441	GSM6390442	GSM6390443	GSM6390444	GSM6390445	GSM6390446 \
	GSM6390447	GSM6390448	GSM6390449	GSM6390450	GSM6390451	GSM6390452	GSM6390453	GSM6390454	GSM6390455 \
	GSM6390456	GSM6390457	GSM6390458	GSM6390459	GSM6390460	GSM6390461	GSM6390462 -O $wd/data/
	
	    #################### CONTROLE DE QUALIDADE ########################
	    
	
#Criei a pasta fastq_MM antes, para poder enviar os dados fastq na proxima etapa. 
mkdir fastq_MM

	#1) Transformar os arquivos SRA em FastQ.
	#OBS: Em bash, na linha de comando, deve ser digitado o codigo em uma unica linha e sem "\", que indica uma quebra
	
	parallel -j 5 'fastq-dump --split-files --gzip \
	-O $wd/fastq_MM/' ::: $wd/data/*/*.sra #Executar em todos os arquivos .sra do diretorio e enviar para outro 
	
	#COMO HAVIA DADO ERRO AO RODAR O CODIGO, INSERI TODO O CAMINHO 
	#parallel -j 5 'fastq-dump --split-files --gzip -O /home/scratch90/bleocata_20240912/PRJEB37100/fastq_MM/' ::: /home/scratch90/bleocata_20240912/PRJEB37100/data/*/*.sra


# 2) CRIANDO UMA LISTA DE ID NAMES PARA CADA AMOSTRA

# listar todos os arquivos .fastq.gz no diretório fastq_MM sem incluir subdiretórios. 
	find "$wd/fastq_MM/" -maxdepth 1 -mindepth 1  -name '*.fastq.gz' | \ 
	xargs -I {} basename {} > $OUTPUT_DIR/samples_ID.list 	#remove o caminho completo de cada um dos arquivos e salva os nomes dos arquivos no arquivo samples_ID.list.

	SAMPLES=$(find "$wd/data/fastq_MM/" -maxdepth 1 -mindepth 1  -name '*.fastq.gz' | xargs -I {} basename {} ) #Atribui os arquivos (sem o path) a variavel SAMPLE
	echo "$SAMPLES" | cut -d'_' -f1| sort -u  > $OUTPUT_DIR/samples_ID.list #Imprime o nome dos arquivos (antes do primeiro _) armazenados na variavel SAMPLE, e organiza em ordem alfabetica excluindo duplicatas
	
	
	
	############### PRE-PROCESSAMENTO ###############


pre_processing (){

    # CONCAT LANES FASTQ
    # Percorra o array recém-criado de amostras e agrupe-as pelo número de leituras
    #for i in `ls */*/*.fastq.gz| cut -d '/' -f 3| cut -d '_' -f1 | sort -u`
    #do
    #    cat */*/${i}*_R1_*.fastq.gz > mergeFastq/"$i"_1.fastq.gz
    #    cat */*/${i}*_R2_*.fastq.gz > mergeFastq/"$i"_2.fastq.gz
        #rm fastq_files/${i}*
    #done
	
	
	   # 1º_FASTQC: Avaliação da qualidade das leituras
	  mkdir $wd/result_fastqc_MM
    parallel -j 10 'fastqc -O $wd/result_fastqc_MM/' ::: $wd/fastq_MM/*.fastq.gz
    
    # 2º_MULTIQC: Relatório resumido de qualidade 
	mkdir $wd/result_multiqc_MM
	multiqc $wd/result_fastqc_MM/ -O $wd/result_multiqc_MM/
    
  # 3º_TRIMMOMATIC:Remoção de leituras de baixa qualidade e adaptadores > caso nao passe pelo MULTIQC
  
  mkdir -p trimmed

  java -jar /home/tools/manual/Trimmomatic-0.33/trimmomatic-0.33.jar  PE  \  #Executa o trimmomatic para Pair-End
  $wd/fastq_MM/ERR3976_1.fastq.gz $wd/fastq_MM/ERR3976_2.fastq.gz \  #imput sequencia forward e reverse
  $wd/trimmed/ERR3976_1P.fq.gz $wd/trimmed/ERR3976_1U.fq.gz \  #output das sequencias forward
  $wd/trimmed/ERR3976_2P.fq.gz $wd/trimmed/ERR3976_2U.fq.gz \  #output das sequencias reverse
  ILLUMINACLIP:/home/tools/manual/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
  
  while read r1 r2; do 
  done < all.txt
  
  #Remove sequencias especificas/ retira bases de baixa qualidade no inicio e final da read/
    
  #3º_ Cutadapt:Remoção de leituras de baixa qualidade e adaptadores > caso nao passe pelo MULTIQC
	cutadapt -g Primer_F=AGCTTAGCTAGCTACG -G Primer_R=CGTAGCTAGCTAAGCT \
         -a "A{10}" -A "T{10}" \
         -o trimmed_1.fastq.gz -p trimmed_2.fastq.gz \
         input_1.fastq.gz input_2.fastq.gz
         
         
                      ####NOVAMENTE PRE-PROCESSAMENTO####

    # 2º_FASTQC:
    #mkdir pos_fastqc
    #fastqc TOY.mergeFastq/*.fastq.gz -o pos_fastqc
    # 2º_MULTIQC
    #mkdir pos_fastqc/multiqc_output
    #multiqc pos_fastqc/ -o pos_fastqc/multiqc_output
	
	# NÃO REALIZEI A ETAPA ACIMA, pois não rodei o Trimmomatic
}
export -f pre_processing #O codigo so vai rodar quando chamar a funcao 


# STAR BUILD INDEX GENOME
buildIndexGenome (){ 
  echo "" >> $OUTPUT_LOG
  echo ">>>>>> Contruindo o genoma de referencia <<<" >> $OUTPUT_LOG
  date >> $OUTPUT_LOG

    # mkdir reference
    # # download fasta and gtf
    # wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz -P reference/
    # wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz -P \
	# reference/
    # gunzip reference/*.gz

    # Generating genome indexes.
   ${STAR} \
 --runThreadN ${JOBS} \
 --runMode genomeGenerate \
 --genomeDir "${SCRATCH60}genomeSTAR/" \
 --genomeFastaFiles  "${SCRATCH60}reference/GRCh38.primary_assembly.genome.fa" \
 --sjdbGTFfile "${SCRATCH60}reference/gencode.v44.primary_assembly.annotation.gtf" \
 --sjdbOverhang 50  > "${SCRATCH60}genomeSTAR/buildIndexGenome.out" 2> "${SCRATCH60}genomeSTAR/buildIndexGenome.log"
}
export -f buildIndexGenome
	
runSTAR (){
        local ID=$1

    echo "" >> $OUTPUT_LOG
    echo ">>>>>> Executando runSTAR para Amostra: "$ID" <<<" >> $OUTPUT_LOG
    date >> $OUTPUT_LOG

    # Running mapping jobs
    ${STAR} \
    --runThreadN ${JOBS} \
        --runMode alignReads \
    --genomeDir "${SCRATCH60}/genomeSTAR/" \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat \
    --readFilesIn "${SCRATCH60}/fastq_chatila/${ID}_1.fastq.gz," "${SCRATCH60}/fastq_chatila/${ID}_2.fastq.gz" \
    --outFileNamePrefix "${OUTPUT_DIR}/runSTAR/${ID}." > "${OUTPUT_DIR}/runSTAR/${ID}.mapped.log"
        }
export -f runSTAR

xargs -a $OUTPUT_DIR/samples_ID.list -t -n1 -P${JOBS} bash -c 'runSTAR  "$@"' 'runSTAR'	

# PREPARANDO PARA PROCESSAR OS ARQUIVOS BAM E GERAR OS ARQUIVOS "coverage" e "stats"
	# Diretório contendo os arquivos BAM
	BAM_DIR="/home/scratch60/etores_ChatilaRNASeq_2024_07_08/Result_Chatila_RNASeq/runSTAR/"
	COVERAGE_DIR="/home/scratch60/etores_ChatilaRNASeq_2024_07_08/Result_Chatila_RNASeq/coverage_reports2"
	STATS_DIR="/home/scratch60/etores_ChatilaRNASeq_2024_07_08/Result_Chatila_RNASeq/stats_reports2"

	# Criação dos diretórios de saída, se não existirem
	mkdir -p $COVERAGE_DIR $STATS_DIR

	# Função para processar um arquivo BAM
process_bam() {
    local bam_file=$1
    local base_name=$(basename "$bam_file" .Aligned.sortedByCoord.out.bam)
    local coverage_report="$COVERAGE_DIR/${base_name}_coverage2.txt"
    local stats_report="$STATS_DIR/${base_name}_stats2.txt"

        echo "Processando $bam_file..."  # Adicionar log

    # Ordenar o arquivo BAM. Não precisa realizar, pois no passo anterior já foi gerado um arquivo sorted
	#".Aligned.sortedByCoord.out.bam" este aqui. Ficar atento que se usar "-n" ordena por nome
    #samtools sort -@ 4 -m4G -n "$bam_file" -o "$sorted_bam"

    # Criar o índice do arquivo BAM ordenado
    samtools index "$bam_file"
	if [ $? -ne 0 ]; then
		echo "Erro ao indexar $bam_file"
		return 1
	fi

    # Gerar o relatório de cobertura
    samtools coverage "$bam_file" > "$coverage_report"

    # Gerar o relatório de estatísticas
    samtools stats "$bam_file" > "$stats_report"

    echo "Processamento completo para: $base_name"

        }

export -f process_bam
export BAM_DIR COVERAGE_DIR STATS_DIR

# Encontrar todos os arquivos BAM e processá-los em paralelo
find "$BAM_DIR" -name "*.Aligned.sortedByCoord.out.bam" | parallel -j 5 process_bam {}

echo "Processamento paralelo completo!"



# RSEM - PREPARAR REFERENCE GENOME E NORMALIZAR
	mkdir RSEM_reference
	# RSEM PREPARAR PARA NORMALIZAR
	rsem-prepare-reference --gtf reference/gencode.v44.primary_assembly.annotation.gtf \
	--bowtie reference/GRCh38.primary_assembly.genome.fa RSEM_reference/human_gencode
		
	# NORMALIZAR USANDO RSEM PARA TPM

#!/bin/bash

# Diretório onde os arquivos BAM estão armazenados
BAM_DIR="/home/scratch60/etores_ChatilaRNASeq_2024_07_08/Result_Chatila_RNASeq/runSTAR/"
# Nome da referência ou prefixo da referência do RSEM
REFERENCE_NAME="/home/scratch60/etores_ChatilaRNASeq_2024_07_08/RSEM_reference/human_gencode"
# Diretório para armazenar os resultados do RSEM
OUTPUT_DIR="/home/scratch90/etorres_LARC/Chatila_Results/"

# Cria o diretório de saída, se não existir
mkdir -p "$OUTPUT_DIR"

# Função para processar um único arquivo BAM
process_rsem() {
    BAM_FILE="$1"
    SAMPLE_NAME=$(basename "$BAM_FILE" .Aligned.toTranscriptome.out.bam)
    SAMPLE_OUTPUT_DIR="$OUTPUT_DIR/${SAMPLE_NAME}_rsem"

    echo "Processando $SAMPLE_NAME..."

    # Cria o diretório de saída para a amostra, se não existir
    mkdir -p "$SAMPLE_OUTPUT_DIR"

    # Comando para rodar o RSEM
    rsem-calculate-expression --alignments --calc-ci --num-threads ${JOBS}*2 --paired-end --bam --no-bam-output \
    "$BAM_FILE" "$REFERENCE_NAME" "$SAMPLE_OUTPUT_DIR"

    if [ $? -eq 0 ]; then
        echo "Processamento de $SAMPLE_NAME concluído."
    else
        echo "Erro ao processar $SAMPLE_NAME" >&2
    fi
}

export -f process_rsem  # Exporta a função para que o GNU Parallel possa utilizá-la
export BAM_DIR REFERENCE_NAME OUTPUT_DIR   # Exportar as pastas locais

# Encontra todos os arquivos BAM e processa-os em paralelo
find "$BAM_DIR" -name "*.Aligned.toTranscriptome.out.bam" -print0 | parallel -0 -j 5 process_rsem {}

echo "Todos os samples foram processados!"

#> > > > > End Pipeline < < < < <