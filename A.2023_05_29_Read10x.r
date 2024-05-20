#Working on multi-batch, aggregated expression matrix for J1994 single cell data. 
#Attempting to apply Batchelor to batch-correct and alter PCA reduction. 
#Then will cluster and define cell types in aggregate. 
#Prior approaches I attempted included (a) Azimuth alone applied in parallel to 
# each sequencing sample, (b) Azimuth + manual annotations on each sequencing file. 
#I sought to avoid analyzing aggregate matrix but will try now despite my attempt. 
# 29 May 2023 AAG. 

setwd("~/salvage_tmp/final/") #Parent directory of this code repository. 
library(Seurat)
library(parallel)
library(dplyr)
#Load in aggregate data matrix. ---------------------------------------------
aggr = Read10X_h5(filename = './raw_data/scRNA/aggr/outs/count/filtered_feature_bc_matrix.h5')
seu = CreateSeuratObject(counts = aggr)

#I need to compare to files 
meta = read.csv('./raw_data/scRNA/J1994_SC_run2_sampleIDs.csv')
meta$Sample.Name = gsub('-', '_', meta$Sample.Name)
sc_dirs_roots = c('AH_1GE_B9', 
                  'AH_2GE_B10', 
                  'AH_3GE_B11', 
                  'AH_4GE_B12', 
                  'AH_5GE_C1', 
                  'AH_6GE_C2', 
                  'AH_7GE_C3', 
                  'AH_8GE_C4', 
                  'C1D1GE', 
                  'C2D1GE',
                  'C3D1GE', 
                  'C6D1GE')


sc_dirs = paste0('./raw_data/scRNA//count/', sc_dirs_roots, '/filtered_feature_bc_matrix')

indiv_seurats = mclapply(sc_dirs, mc.cores = 5, FUN = function(x){
  d = Read10X(x)
  bcs = colnames(d)
  return(bcs)
})

names(indiv_seurats) <- sc_dirs_roots

indiv_seurats <- lapply(indiv_seurats, FUN = function(x){
  df = data.frame(Barcode = x)
  return(df)
})
bc_key = bind_rows(indiv_seurats, .id = 'Sample')

bc_key$Batch = 1
bc_key$Batch[grep('AH_',bc_key$Sample)] = 2
idx = match(bc_key$Sample, meta$Sample.Name)
bc_key$Patient = meta$Subject.ID[idx]
bc_key$Patient[is.na(bc_key$Patient)] = 'J1994.001'
bc_key$Timepoint = meta$Timepoint[idx]
bc_key$Timepoint[bc_key$Sample %in% c('C2D1GE', 'C3D1GE', 'C6D1GE')] = 'Post'
bc_key$Timepoint[bc_key$Sample %in% 'C1D1GE'] = 'Pre'

bc_key$Barcode = substr(bc_key$Barcode,1,nchar(bc_key$Barcode) - 2)
idx = match(bc_key$Sample, sc_dirs_roots)
bc_key$Barcode = paste0(bc_key$Barcode, '-', idx)
bc_key$VDJ_beta = FALSE
bc_key$VDJ_alpha = FALSE
#For reference, also load in collected VDJ data from all samples. -----------
vdj_dirs = c('AH_1T_B1', 
             'AH_2T_B2', 
             'AH_3T_B3', 
             'AH_4T_B4', 
             'AH_5T_B5', 
             'AH_6T_B6', 
             'AH_7T_B7', 
             'AH_8T_B8', 
             'C1D1T', 
             'C2D1T', 
             'C3D1T', 
             'C6D1T'
)
vdj_dirs = paste0('./raw_data/scVDJ/', vdj_dirs, '.tsv')

vdj_data = mclapply(vdj_dirs, mc.cores = 5, FUN = function(x){
  d = read.delim(x)
  return(d)
})

names(vdj_data) = sc_dirs_roots

for(f in 1:12){ #Annotate which barcodes in single cell had beta chains in VDJ data. 
  vdjs = vdj_data[[f]]
  beta_barcodes = unique(vdjs$cell_id[grep('TRBV', vdjs$v_call)])
  beta_barcodes = substr(beta_barcodes, 1, nchar(beta_barcodes) - 2)
  beta_barcodes = paste0(beta_barcodes, '-', f)
  idx = match(beta_barcodes, bc_key$Barcode)
  bc_key$VDJ_beta[idx] = TRUE
  
  alpha_barcodes = unique(vdjs$cell_id[grep('TRAV', vdjs$v_call)])
  alpha_barcodes = substr(alpha_barcodes, 1, nchar(alpha_barcodes) - 2)
  alpha_barcodes = paste0(alpha_barcodes, '-', f)
  idx = match(alpha_barcodes, bc_key$Barcode)
  bc_key$VDJ_alpha[idx] = TRUE
}
  
  
write.csv(bc_key, file = paste0('./code_outputs/', Sys.Date(), '_aggregate-meta.csv'))
rownames(bc_key) = bc_key$Barcode


#Relate metadata to seurat object.
m = seu@meta.data
m = cbind(m, bc_key)
seu@meta.data = m
save(seu, file = paste0('./code_outputs/',Sys.Date(), '_AggrSeurat.rda'))

