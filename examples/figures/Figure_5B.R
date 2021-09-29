library("DESeq2")

cell_types <- c("Pericytes", "Endothelial Cells", "Microglia", "Myoc+ Astrocytes", "Myoc- Astrocytes", "Excitatory Neurons", "Inhibitory Neurons", "Perivascular Macrophages", "Vascular SMCs", "Committed Oligodendrocytes", "Mature Oligodendrocytes", "VLMCs", "OPCs")

conditions <- c("3DPI", "25DPI")

for (condition in conditions) {

        metadata <- read.csv(sprintf("sample_metadata_%s.csv", condition), header = TRUE, sep = ",")

        for (cell_type in cell_types) {
            
            
                if (!file.exists(sprintf("%s_gene_counts_%s.csv", cell_type, condition))) {
                    next
                }

                countData <- read.csv(sprintf("%s_gene_counts_%s.csv", cell_type, condition), header = TRUE, sep = ",")

                dds <- DESeqDataSetFromMatrix(countData=countData, colData=metadata, design= ~ Condition, tidy = TRUE)
                dds <- DESeq(dds)
                res <- results(dds, contrast = c("Condition", condition, "Control"), alpha = 0.05)
                summary(res)
                res <- res[order(res$padj),]
                head(res)

                write.csv(res, file=sprintf("deseq2_out/%s_%s_deseq2.csv", cell_type, condition), row.names = TRUE)
        }
}

warnings()