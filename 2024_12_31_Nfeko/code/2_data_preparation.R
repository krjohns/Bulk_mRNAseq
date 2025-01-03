##setup----
# libraries
library(tidyverse) # all-purpose data analysis tools
library(fs) # filepath manipulation
library(Matrix) # for working with matrices (especially sparse ones)
library(cowplot) #ggplot add on (save_plot)
library(readxl) #for reading excel docs
library(scales) #for sci notation

# directories
wd <- getwd()
input_dir <- "/Volumes/nchevrier/projects/katej/projects/Bulk_mRNAseq/katej_12.24_bcbiogenome/STARSolo_Output/GeneFull/Filtered" # where bcbio output and meta data matrices live on project2 server
output_dir <- fs::path(wd, "2024_12_31_Nfeko/input/data/2_formatted_data/bcbio_genome")# where prepped data will live (metadata/counts csvs etc)

##read in sample info for de-multiplexing----
sample_info <- read_xlsx(fs::path(wd, "2024_12_31_Nfeko/input/other/20241218_nfeKO_rnaseq_sample_sheet.xlsx"))
#create meta data table
meta <- sample_info %>%
        separate(Sample_Name, into = c("Target", "Guide", "Replicate"), sep = "_", remove = F)
                 
#read in counts matrices for DCs and T cells (combined matrix)---
#read in tagcounts matrices
counts = readMM(fs::path(input_dir, 'matrix.mtx')) 
# add rownames and colnames
rownames(counts) = read_lines(fs::path(input_dir, 'features.tsv'))
colnames(counts) = read_lines(fs::path(input_dir, 'barcodes.tsv'))

#Shorten rownames to only gene name (remove ENSEMBL ID info)
rownames(counts) = gsub(".*\t","",gsub("\tGene.*","",rownames(counts)))

#Replace barcodes with sample names----
#convert to non-sparse format
NS_counts = as.matrix(counts)

# reorder to match metadata & pull out barcodes used in exp
NS_counts = NS_counts[, meta$Barcode]

# switch colnames to full name instead of barcode
colnames(NS_counts) = plyr::mapvalues(colnames(NS_counts),  from = meta$Barcode, to = meta$Sample_Name)


#reorder samples so replicates are clustered together for easier viz in morpheus
meta_order = meta %>% 
             mutate(Target = factor(Target, levels = c("sgCtrl", "sgNfe2l1", "sgNfe2l2", "doubleKO"))) %>%
             arrange(Target, Guide, Replicate)
  
NS_counts_order = NS_counts[,meta_order$Sample_Name]
  
#save separate count & meta data in case its useful later----
write.csv(as.matrix(NS_counts_order), fs::path(output_dir, "20241231_Nfekoseq_counts.csv"))
write.csv(as.matrix(meta_order), fs::path(output_dir, "20241231_Nfekoseq__metadata.csv"))


#----checking unfiltered matrix to see # of reads for unused barcodes
counts_unfilt = readMM(fs::path(input_dir, '../raw/matrix.mtx.gz')) 
# add rownames and colnames
rownames(counts_unfilt) = read_lines(fs::path(input_dir, '../raw/features.tsv.gz'))
colnames(counts_unfilt) = read_lines(fs::path(input_dir, '../raw/barcodes.tsv.gz'))

#Shorten rownames to only gene name (remove ENSEMBL ID info)
rownames(counts_unfilt) = gsub(".*\t","",gsub("\tGene.*","",rownames(counts_unfilt)))

NS_counts_unfilt = as.matrix(counts_unfilt) 
 
full_meta = data.frame(Barcode = colnames(NS_counts_unfilt)) %>%
            mutate(Barcode_used = ifelse(Barcode %in% meta$Barcode, T, F))

read_sum = data.frame(Reads = colSums(NS_counts_unfilt)) %>%
           rownames_to_column("Barcode") %>%
           full_join(full_meta, by = "Barcode") %>%
           mutate(Sample_Name = plyr::mapvalues(Barcode, from = meta$Barcode, to = meta$Sample_Name),
                  Sample_Name = factor(Sample_Name, levels = meta_order$Sample_Name))

#raw only has 24 barcodes in it so must be already filtered matrix...           
#save barcode QC plots
p<-ggplot(read_sum, aes(x = Sample_Name, y = Reads)) +
  scale_y_continuous(labels = scales::scientific) +
  xlab(NULL) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 45, hjust =1))

#reads might be too low to get meaninful DE data...
save_plot(fs::path(output_dir, "../../../../output/bcbio_genome/QC/read_by_sample.pdf"), p)
