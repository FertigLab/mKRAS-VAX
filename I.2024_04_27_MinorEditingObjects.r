#Added bacth annotations to seurat objects. 
#Removed nonT cells that were still present in filtered seurat objects. 
#Previously objects had all cells belonging to 'T cell clusters' but some 
#inidvidual cells were later annotated as non T cells. 
#AAG 4 January 2024. 

setwd("/media/aag7319/WDRed/ZZZ_Homolig/data_exports/Executable-Code/DiffEx-TCR-Mapping/")
load( './seurat-objects/2024-04-18_Batch1TCells-KRASAnnot.rda') #seu1
load('./seurat-objects/2024-04-18_Batch2TCells-KRASAnnot.rda') #seu2
batch1 = read.csv('./cell-annotations/2024-01-04_batch1-tcell-annot.csv')
batch2 =  read.csv('./cell-annotations/2024-01-04_batch2-tcell-annot.csv')
batch_annot = rbind(batch1,batch2)

m1 = seu1@meta.data
idx = match(batch1$barcode, rownames(m1))
m1$batch_annotation = NA
m1$batch_annotation[idx] = batch1$annot
m2 = seu2@meta.data
idx = match(batch2$barcode, rownames(m2))
m2$batch_annotation = NA
m2$batch_annotation[idx] = batch2$annot
write.csv(m1, file = './2024_04_27_Batch1TCells-KRASAnnot-Metadata.csv', quote = FALSE, row.names = TRUE)
write.csv(m2, file = './2024_04_27_Batch2TCells-KRASAnnot-Metadata.csv', quote = FALSE, row.names = TRUE)

seu1@meta.data = m1
seu2@meta.data = m2
acceptable = c("CD4 Naive" ,"CD4 Proliferating" ,"CD4 TCM", "CD4 TEM","CD8 Naive", 
               "CD8 Proliferating" , "CD8 TCM"  ,"CD8 TEM"  ,"dnT", "MAIT" , "CD4 CTL")
seu1_final = seu1[,seu1@meta.data$batch_annotation %in% acceptable]
seu2_final = seu2[,seu2@meta.data$batch_annotation %in% acceptable]
save(seu1_final,seu2_final, file = './2024-04-27_TCells-BothBatches-final.rda')

m1 = seu1_final@meta.data
m2 = seu2_final@meta.data
write.csv(m1, file = './2024_04_27_Batch1TCells-KRASAnnot-MetadataFinal.csv', quote = FALSE, row.names = TRUE)
write.csv(m2, file = './2024_04_27_Batch2TCells-KRASAnnot-MetadataFinal.csv', quote = FALSE, row.names = TRUE)