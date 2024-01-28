#Data download and DEA

GSE121212
gse <- GEOquery::getGEO("GSE121212", GSEMatrix = TRUE)
metadata <- gse$GSE121212_series_matrix.txt.gz |> Biobase::pData()
expr_path <- GEOquery::getGEOSuppFiles(GEO="GSE121212", fetch_files=TRUE)
expression_data <- data.table::fread(rownames(expr_path))

expression_data <- expression_data |>
  dplyr::distinct(V1, .keep_all = TRUE)  |>
  dplyr::rename(gene=V1)
 


GSE121212_obj <- list(expression_data=expression_data, metadata=metadata)

saveRDS(GSE121212_obj, "./GSE121212/GSE121212_object.rds")


GSE121212_data <- readRDS("./GSE121212/GSE121212_object.rds")
GSE121212_data_dds <- GSE121212_data$expression_data |>
  tibble::column_to_rownames(var="gene")
rownames(GSE121212_data$metadata ) <- GSE121212_data$metadata$title
GSE121212_data$metadata$`patient's condition:ch1`  <- factor(GSE121212_data$metadata$`patient's condition:ch1` )

test <- identical(GSE121212_data$metadata$title,colnames(GSE121212_data_dds))

coldata <- GSE121212_data$metadata[,c("patient's condition:ch1","skin type:ch1" )] 
coldata$`patient's condition:ch1` <- factor(coldata$`patient's condition:ch1`)
coldata$`skin type:ch1` <- factor(coldata$`skin type:ch1`)
coldata <- coldata |>
  dplyr::rename("tissue_type" =  "skin type:ch1",
                "condition" = "patient's condition:ch1")
coldata$condition <- relevel(coldata$condition, ref = "CTRL")
dds_obj <- DESeq2::DESeqDataSetFromMatrix(countData = GSE121212_data_dds, 
                                      colData = coldata, 
                                      design = ~condition)


keep <- rowSums(counts(dds_obj) >= 10) >= 37
dds_obj <- dds_obj[keep,]

dds <- DESeq2::DESeq(dds_obj)

res_AD_VS_CTRL <- DESeq2::results(dds, name="condition_AD_vs_CTRL", pAdjustMethod = "BH")
res_AD_VS_CTRL <- res_AD_VS_CTRL |>
  as.data.frame() |>
  dplyr::mutate(gene = rownames(res_AD_VS_CTRL))
res_PSO_VS_CTRL <- DESeq2::results(dds, name="condition_PSO_vs_CTRL", pAdjustMethod = "BH") |>
  as.data.frame() |>
  dplyr::mutate(gene = rownames(res_PSO_VS_CTRL))


saveRDS(res_AD_VS_CTRL, "./GSE121212/resDEA_AD_VS_CTRL.rds")
saveRDS(res_PSO_VS_CTRL, "./GSE121212/resDEA_PSO_VS_CTRL.rds")





