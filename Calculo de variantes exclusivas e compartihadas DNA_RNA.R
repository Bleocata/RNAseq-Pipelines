################# CHAMADA DE VARIANTES SOMATICAS EM DADOS DE RNASEQ #############

#1_Carregar bibliotecas:
library(dplyr)
library(readr)
library(tidyr)
library(maftools)
library(data.table)
library(VennDiagram)

#2_Carregar os arquivos: 
setwd("C:/Users/beatr/Sirio Libanes/Variant_calling")
MUT_DNA <- read_tsv("SM_DNA.tsv")
MUT_DNA2 <- read_delim("Final_GATK.4.6_annotated.txt", guess_max = 1000)

# Leitura com controle de delimitador
MUT_RNA <- read_delim("Final_GATK.4.6_annotated.txt", guess_max = 1000)

# Verificar número de colunas
ncol(data)


MUT_RNA <- fread("Final_GATK.4.6_annotated.tsv", header = T, sep = "\t")
var.annovar.maf = annovarToMaf(annovar = "Final_GATK.4.6_annotated.tsv", refBuild = 'hg38', 
                               tsbCol = 'sample', table = 'refGene')

#Criar a tabela DNA com Index para comparacao 

dna <- MUT_DNA %>%
  unite(col = "Index",Sample,CHROM,POS, sep = "_")

#Criar a tabela DNA com Index para comparacao 
rna <- MUT_RNA %>% select(-End)
rna <- rna %>% rename(POS = Start)
rna <- rna %>%
  unite(col = "Index",sample,Chr,POS, sep = "_")

# Identificar variantes presentes no DNA mas ausentes no RNA
present_in_dna_not_rna <- anti_join(dna, rna, by = "Index")

# Identificar variantes presentes no RNA mas ausentes no DNA
present_in_rna_not_dna <- anti_join(rna, dna, by = "Index")

# Identificar variantes presentes em ambas as tabelas
present_in_both <- inner_join(dna, rna, by = "Index")

# Resultados
print("Present in DNA but not RNA:")
print(present_in_dna_not_rna)

print("Present in RNA but not DNA:")
print(present_in_rna_not_dna)

print("Present in both DNA and RNA:")
print(present_in_both)

#calculos e porcentagem :
# Calcular números de variantes exclusivas e compartilhadas
#variante = linha(row)
total_dna <- nrow(dna)  # Total de variantes no DNA
total_rna <- nrow(rna)  # Total de variantes no RNA

exclusive_dna <- nrow(present_in_dna_not_rna)  # Variantes exclusivas do DNA
exclusive_rna <- nrow(present_in_rna_not_dna)  # Variantes exclusivas do RNA
shared_variants <- nrow(present_in_both)       # Variantes compartilhadas

# Cálculos de porcentagens
perc_exclusive_dna <- (exclusive_dna / total_dna) * 100
perc_exclusive_rna <- (exclusive_rna / total_rna) * 100
perc_dna_in_rna <- (shared_variants / total_dna) * 100  # % do DNA identificado no RNA
perc_rna_in_dna <- (shared_variants / total_rna) * 100  # % do RNA identificado no DNA

# Exibir os resultados
cat("Porcentagens calculadas:\n")
cat(sprintf("1. Porcentagem de variantes exclusivas do DNA: %.2f%%\n", perc_exclusive_dna))
cat(sprintf("2. Porcentagem exclusivas do RNA: %.2f%%\n", perc_exclusive_rna))
cat(sprintf("3. Porcentagem das variantes do DNA encontradas no RNA: %.2f%%\n", perc_dna_in_rna))
cat(sprintf("4. Porcentagem das variantes do RNA com suporte no DNA (verdadeiro positivo): %.2f%%\n", perc_rna_in_dna))

############Criar diagrama de venn para representacao grafica dos resultados obtidos ##################


# Criar o diagrama de Venn com contornos coloridos e preenchimento transparente
venn.plot <- draw.pairwise.venn(
  area1 = total_dna,                 # Total de variantes no DNA
  area2 = total_rna,                 # Total de variantes no RNA
  cross.area = shared_variants,      # Variantes compartilhadas
  category = c("DNA", "RNA"),        # Rótulos das categorias
  fill = c("transparent", "transparent"), # Preenchimento transparente
  col = c("deepskyblue3", "deeppink3"),          # Cor dos contornos
  lwd = 2,                           # Largura do contorno
  cex = 2,                           # Tamanho do texto interno
  cat.cex = 2,                       # Tamanho do texto das categorias
  cat.col = c("deepskyblue3", "deeppink3")       # Cor do texto das categorias
)

# Exibir o gráfico no RStudio
grid.draw(venn.plot)

# Salvar o gráfico
png("venn_diagram_adjusted.png", width = 800, height = 600)
grid.draw(venn.plot)
dev.off()
