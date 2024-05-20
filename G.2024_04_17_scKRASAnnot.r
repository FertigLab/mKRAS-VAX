#Annotate in seurat objects anti-KRAS TCRs from expansion data. 
#AAG 3 January 2024 

library(Seurat)
library(ggplot2)
library(ggprism)
library(RColorBrewer)
library(ragg)
library(dplyr)
library(parallel)

setwd("~/salvage_tmp/final/")
source('./_rfunctions/ProcessAdaptiveFile.R')

#metadata for sc repertoire data
samples = c(paste0('AH_', c(1:8), 'T_B', c(1:8)), 
            paste0('C', c(1,2,3,6), 'D1T'))
meta = read.csv('./raw_data/scRNA/J1994_SC_run2_sampleIDs.csv')
meta1 = data.frame(
  Subject.ID = 'J1994.001', 
  Cycle = c('Cycle 1 - Week 0', 'Cycle 2 - Week 0', 'Cycle 3 - Week 0', 'Cycle 6 - Week 0'), 
  Sample.Name = samples[9:12], 
  Publication.ID.to.be.used.in.ms = NA, 
  Timepoint = c('Pre', 'Post', NA, NA)
)
meta = rbind(meta, meta1)
meta = meta[-c(1:8),] #This is a 12 row metadata frame corresponding to indices of scVDJ file
meta$Sample.Name = gsub('-', '_', meta$Sample.Name)
rm(meta1)
#Read in Tcell annotations--------
#These annotations are prepared via combination of manual curation and azimuth 
#on aggregated counts matrices by batch. 
batch1 = read.csv('./code_inputs/2024-01-04_batch1-tcell-annot.csv')
batch2 = read.csv('./code_inputs/2024-01-04_batch2-tcell-annot.csv')
batch_annot = rbind(batch1,batch2)
batch_annot = batch_annot[-which(batch_annot$annot == 'Eryth'),]

#Read in all scVDJ data-----------
samples = c(paste0('AH_', c(1:8), 'T_B', c(1:8)), 
            paste0('C', c(1,2,3,6), 'D1T'))
meta = read.csv('./raw_data/scRNA/J1994_SC_run2_sampleIDs.csv')
meta1 = data.frame(
  Subject.ID = 'J1994.001', 
  Cycle = c('Cycle 1 - Week 0', 'Cycle 2 - Week 0', 'Cycle 3 - Week 0', 'Cycle 6 - Week 0'), 
  Sample.Name = samples[9:12], 
  Publication.ID.to.be.used.in.ms = NA, 
  Timepoint = c('Pre', 'Post', NA, NA)
)
meta = rbind(meta, meta1)
meta = meta[-c(1:8),] #This is a 12 row metadata frame corresponding to indices of scVDJ file
meta$Sample.Name = gsub('-', '_', meta$Sample.Name)

bc_key = data.frame(sample = samples, idx = c(1:12))

files = paste0('./raw_data/scVDJ/',dir('./raw_data/scVDJ/'))
scVDJ = mclapply(files, FUN = read.delim, mc.cores = 10)

#Apply Batch Annotations to VDJ data----------------------------------------
bc_end = as.numeric(gsub("[^\\d]+", "", batch_annot$barcode, perl=TRUE))
annot_split = split(batch_annot, bc_end)
names(annot_split) = samples

process <- function(annot, vdj){
  vdj_strippedbc = substr(vdj$cell_id, 1, nchar(vdj$cell_id) - 2)
  annot$barcode = gsub('-.*', '', annot$barcode)
  idx = match(vdj_strippedbc, annot$barcode)
  vdj$batch_annotation = annot$annot[idx]
  return(list(vdj))
}

scVDJ_annot = mcmapply(FUN = process, annot = annot_split, vdj = scVDJ, mc.cores = 10)
names(scVDJ_annot) = samples

percent_mapped = unlist(lapply(scVDJ_annot, FUN = function(x){
  tot_ncells = length(unique(x$cell_id))
  tot_mapped = length(unique(x$cell_id[!is.na(x$batch_annotation)]))
  return(round(tot_mapped/tot_ncells, 3))
}))

meta$nBetaVDJ_cells = unlist(lapply(scVDJ_annot, FUN = function(x){
  x = x[grep('TRB', x$v_call),]
  tot_ncells = length(unique(x$cell_id))
  return(tot_ncells)
}))


patient_scVDJ_files = list(
  scVDJ_annot[9:12],
  scVDJ_annot[1:2], 
  scVDJ_annot[3:4], 
  scVDJ_annot[5:6], 
  scVDJ_annot[7:8]
)
patient_ids = paste0('J1994-', c('01', '02', '04', '08', '13'))
names(patient_scVDJ_files) = patient_ids

# readRDS('/home/aag7319/salvage_tmp/objects/2023-12-06_scVDJ-annot-from2023_06_08_CompleteFEST.rds')
scVDJ = lapply(scVDJ_annot, FUN = function(y){
  y = y[y$productive & y$v_call != '',]
  alpha = y[grep('TRAV', y$v_call),]
  beta = y[grep('TRBV', y$v_call),]
  beta$v_call = ProcessAdaptiveVgenes(beta$v_call)
  alpha$v_call = ProcessAdaptiveVgenes(alpha$v_call)
  beta$uid = paste0(beta$v_call, '_', beta$junction_aa)
  y = list(alpha = alpha, beta = beta)
  return(y)
})

betas = lapply(scVDJ, FUN = function(x){return(x$beta)})
betas = bind_rows(betas, .id = 'Sample')

#Load single cell objects##########################
tmp = load('./code_inputs/2024-01-03_batch2filt_TCellOnly.rda')
seu2 = filt

tmp = load('./code_inputs/2024-01-03_batch1filt_TCellOnly.rda')
seu1 = filt

m1 = seu1@meta.data
u1 = seu1@reductions$umap@cell.embeddings 
u1 = as.data.frame.matrix(u1)
m1 = cbind(m1,u1)
m2 = seu2@meta.data
u2 = seu2@reductions$umap@cell.embeddings 
u2 = as.data.frame.matrix(u2)
m2 = cbind(m2,u2)

m1$Barcode= substr(m1$Barcode, 1, 16)
m2$Barcode = substr(m2$Barcode, 1, (nchar(m2$Barcode) - 2))
is_annot = betas[!is.na(betas$batch_annotation) & betas$is_cell, ]
yields = table(is_annot$Sample) #The number of annotated T cells per sample. Used for Tot Rep Freq
# yield1 = table(m1$Sample)
# yield2 = table(m2$Sample)
# yields = c(yield2, yield1); rm(yield1, yield2)
meta$scSampleName = c(unique(m2$Sample), unique(m1$Sample))
meta$scTCellCount = as.numeric(yields)


#load diffex data, with NO frequency cutoff. ######
MIN_OR = 5; MIN_FREQ = 0.00
pts = c('01', '02', '04', '08', '13')
diffEx_all = lapply(
  paste0('./code_inputs/diffex-data/2023-06-07_J1994-', pts,'_diffEx.tsv'), 
  read.delim, 
  sep = ' '
)
names(diffEx_all) = paste0('J1994.0', pts)
diffEx_all = lapply(diffEx_all, FUN = function(diffEx_all){
  diffEx_all$xr_or = NULL; diffEx_all$xr_baseline = NULL; diffEx_all$xr_all = NULL; diffEx_all$xr = NULL
  diffEx_all$baseline_thresh = NULL; diffEx_all$or_thresh =NULL; diffEx_all$count_thresh = NULL
  return(diffEx_all)
})


agg = bind_rows(diffEx_all, .id = 'Subject')
agg_filt = agg[(agg$or >= MIN_OR & agg$f1 >= MIN_FREQ), ]
tab = table(agg_filt$id, agg_filt$Subject) > 0
public_clones = rownames(tab)[which(rowSums(tab) > 1)]
public_individuals = apply(tab[which(rowSums(tab) > 1),], 1, FUN = function(x){
  paste0(names(x[which(x > 0)]), collapse = ';')
})
diffEx_all = lapply(diffEx_all, FUN = function(de){
  de = de[de$or >= MIN_OR & de$f1 >= MIN_FREQ,]
  tab = table(de$id[de$or >= MIN_OR & de$f1 >= MIN_FREQ], de$Antigen[de$or >= MIN_OR & de$f1 >= MIN_FREQ])
  xr_ids = rownames(tab)[rowSums(tab) > 1]
  xr_ags = apply(tab[rowSums(tab) > 1, ], 1, FUN = function(x){
    ags = paste0(names(x)[which(x > 0)], collapse = ';')
  })
  xr_ags = unlist(xr_ags)
  de$public = de$id %in% public_clones
  idx = match(de$id[de$public], public_clones)
  de$public_individuals = NA
  de$public_individuals[de$public] = public_individuals[idx]
  de$xr = de$id %in% xr_ids
  idx = match(de$id[de$xr], xr_ids)
  de$xrAgs = NA
  de$xrAgs[de$xr] = xr_ags[idx]
  #de = de[match(unique(de$uid), de$uid), ]
  
  return(de)
})

diffex = lapply(diffEx_all, FUN = function(de){
  de = de[match(unique(de$id), de$id),]
  x = data.frame(
    Antigen = de$Antigen, 
    XR = de$xr, 
    TRBV = gsub('_.*', '', de$id),
    CDR3.beta.aa =  gsub('.*_', '', de$id),
    uid = de$id,
    x = de$Antigen
  )
  x$adaptive.tcrb.id = x$uid
  x$uid = paste0(ProcessAdaptiveVgenes(x$TRBV), '_', x$CDR3.beta.aa)
  x$Antigen[x$XR] =de$xrAgs[de$xr]
  x$x[x$XR] = 'XR'
  return(x)
})

diffex = bind_rows(diffex, .id = 'Patient')
# diffex = read.csv('./homolig_inputs/2023_07_25_DiffExTCRGroups/Collected_All.csv')
# diffex = diffex[diffex$TRBV != '',]
# diffex$TRBV = ProcessAdaptiveVgenes(diffex$TRBV)
# diffex$uid = paste0(diffex$TRBV, '_', diffex$CDR3.beta.aa)
# diffex$x = diffex$Antigen; diffex$x[diffex$XR] = 'XR'

###Key Helper functions#####
digestTRBV <- function(uid = NA, trbv = NA){
  #break down uid into TRBV gene. 
  #split into V family, V gene, and V allele components. 
  #return data frame with each field. 
  #used to match unclear V congruence between 10x and Adaptive data. 
  if(all(is.na(uid))) uid = trbv
  N = length(uid)
  if(all(is.na(trbv))) trbv = gsub('_.*', '', uid)
  
  has_allele =  c(1:N) %in%  grep('\\*', trbv)
  has_gene = c(1:N) %in%  grep('\\-',trbv)
  
  allele = vector(mode = 'integer', length = N)
  gene = vector(mode = 'integer', length = N)
  tmp = gsub( '\\*', 'Z', trbv)
  
  allele[has_allele] = gsub('.*Z','',tmp[has_allele])
  gene[has_gene] = gsub( 'Z.*' , '', gsub('.*\\-','',tmp[has_gene]))
  
  
  family =gsub( '\\-.*' , '', gsub('.*BV', '',trbv)) 
  family[!has_gene & has_allele] = gsub( '\\*.*' , '', gsub('.*BV', '',trbv[!has_gene & has_allele])) 
  
  df = data.frame(
    uid = uid, 
    trbv = trbv, 
    family = as.numeric(family), 
    gene = as.numeric(gene), 
    allele = as.numeric(allele)
  )
  return(df)
}
CheckTRBVEquivalence <- function(input1, input2){
  if(length(input1) == length(input2)) N = length(input1)
  v1 = digestTRBV(trbv =input1)
  v2 = digestTRBV(trbv = input2)
  
  good = vector(mode = 'logical', length = N)
  good = (v1$family == v2$family & v1$gene == v2$gene) | 
    (v1$family == v2$family & (v1$gene == 0 | v2$gene == 0) ) | 
    (v1$family == v2$family & v1$gene == v2$gene )
  
  return(good)
}
getKRASBarcodes <- function(subject){
  samples = meta$Sample.Name[meta$Subject.ID == subject]
  tcrs = betas[betas$Sample %in% samples,]
  mapped = which(tcrs$junction_aa %in% diffex$CDR3.beta.aa) 
  if(length(mapped) == 0) return(NA)
  barcodes = data.frame(Sample = tcrs$Sample[mapped], Barcode = tcrs$cell_id[mapped], 
                        uid = tcrs$uid[mapped], 
                        cdr3 = tcrs$junction_aa[mapped],
                        trbv.10x = tcrs$v_call[mapped],
                        batch_annotation = tcrs$batch_annotation[mapped], 
                        is_cell = tcrs$is_cell[mapped])
 # idx = match(barcodes$uid, diffex$uid) #revise. 
  idx = match(barcodes$cdr3, diffex$CDR3.beta.aa)
  barcodes$reactivity = diffex$x[idx]
  barcodes$all.antigens = diffex$Antigen[idx]
  barcodes$adaptive.tcrb.id = diffex$adaptive.tcrb.id[idx]
  barcodes$adaptive.proc.id = diffex$uid[idx]
  barcodes$adaptive.trbv = diffex$TRBV[idx]
  barcodes = barcodes[order(barcodes$uid),]
  
  barcodes$PASS = barcodes$uid == barcodes$adaptive.proc.id
  barcodes$PASS[!(barcodes$PASS) & CheckTRBVEquivalence(barcodes$trbv.10x, barcodes$adaptive.trbv)] = TRUE
  
  #Manually need to check an edge case: 
  #If a diffEx CDR3 has multiple TRBV genes, I may have not found it in single cell 
  #because I only checked scTCR equivalence to ONE of those V-CDR3 pairs. 
  tmp = table(diffex$CDR3.beta.aa, diffex$TRBV) > 0
  tmp2 = rowSums(tmp)
  to.check = names(tmp2)[tmp2 > 1] #CDR3s with multiple TRBV genes in diffex data. 
  barcodes$CHECK = barcodes$cdr3 %in% to.check
  
  barcodes = barcodes[barcodes$PASS | barcodes$CHECK,]
  return(barcodes)
}
SummarizeMapped <- function(subject, b = NA){
  if(all(is.na(b)))  b = getKRASBarcodes(subject)
  idx = match(b$Sample, meta$Sample.Name)
  b$percent = 1/meta$scTCellCount[idx] * 100
  b$Timepoint = meta$Timepoint[idx]
  b$batch_annotation = factor(b$batch_annotation, levels = celltypes )
  b$reactivity = factor(b$reactivity, levels = ags)
  b$Timepoint = factor(b$Timepoint, levels = c('Pre', 'Post'))
  r =  tapply(b, x ~ b$batch_annotation + b$Sample + b$reactivity + b$Timepoint, FUN = function(x){sum(x$percent)})
  df = reshape2::melt(r)
  return(df)
}
###SCRIPT######
subjects = unique(meta$Subject.ID)[c(5,1:4)]
barcodes = lapply(subjects, getKRASBarcodes); names(barcodes) = subjects
barcodes = bind_rows(barcodes, .id = 'Subject')
barcodes$Barcode = substr(barcodes$Barcode, 1, (nchar(barcodes$Barcode) - 2))

#Now manually check the possibility that some failed, flagged barcodes may have a good pairing.
to.check =which(barcodes$CHECK & !barcodes$PASS)
for(i in to.check){
  tmp = barcodes[i, ]
  de = diffex[diffex$CDR3.beta.aa == tmp$cdr3,]
  print(paste0('INDEX: ', i))
  print(tmp)
  print(de)
  cat(c('\n \n \n'))
  
}
#One single TCR has met this rare edge case criteria. 
idx =grep('CASSQDGVATDTQYF', barcodes$cdr3)
print(barcodes[idx,])
print( diffex[diffex$CDR3.beta.aa == tmp$cdr3,])
barcodes$adaptive.tcrb.id[idx] = 'TCRBV03-01/03-02*01_CASSQDGVATDTQYF'
barcodes$adaptive.trbv[idx] = 'TCRBV03-01/03-02*01'
barcodes$adaptive.proc.id[idx] = 'TRBV3-1*01_CASSQDGVATDTQYF'
barcodes$CHECK[idx] = FALSE
barcodes = barcodes[!(barcodes$CHECK & !barcodes$PASS),]
barcodes_spl = split(barcodes, barcodes$Subject)


barcodes$PASS = NULL;barcodes$CHECK = NULL; barcodes$adaptive.proc.id = NULL; 


celltypes = unique(barcodes$batch_annotation)
celltypes = unique(barcodes$batch_annotation)[!is.na(celltypes)]
ags = unique(barcodes$reactivity)
ags = ags[order(ags)]

subjects = subjects

summary = vector(mode = 'list', length = length(subjects))

for( s in 1:length(subjects)){
  summary[[s]] = SummarizeMapped(subjects[s], b = barcodes_spl[[s]])
}
summary = bind_rows(summary)
colnames(summary) = c('celltype', 'Sample', 'Reactivity' , 'Timepoint', 'Freq')
idx = match(summary$Sample, meta$Sample.Name)
summary$Subject = meta$Subject.ID[idx]

#Append KRAS Reactivity to Seurat metadata
#I know I could obviously use match() somehow but its giving me trouble. 
tmp2 = paste0(m2$Sample, '_', m2$Barcode)
tmp1 = paste0(m1$Sample, '_', m1$Barcode)
m2$antiKRAS = FALSE
m2$tcr_uid = NA
m2$reactivity = NA
m1$antiKRAS = FALSE
m1$tcr_uid = NA 
m1$reactivity = NA
for(i in 1:nrow(meta)){
  gesamp = meta$scSampleName[i]
  barinfo = barcodes[barcodes$Sample == meta$Sample.Name[i],]
  barinfo = barinfo[!is.na(barinfo$batch_annotation),]
  bars = barinfo$Barcode
  queries = paste0(gesamp, '_', bars)
  if(length(bars) == 0) next 
  if(meta$Subject.ID[i] == 'J1994.001'){
    idx = match(queries, tmp1)
    m1$tcr_uid[idx] = barinfo$uid
    m1$reactivity[idx] = barinfo$reactivity
    m1$antiKRAS[idx] = TRUE
  }else{
    idx = match(queries, tmp2)
    m2$tcr_uid[idx] = barinfo$uid
    m2$reactivity[idx] = barinfo$reactivity
    m2$antiKRAS[idx] = TRUE
  }
}

write.csv(barcodes, paste0('./code_outputs/',Sys.Date(), '_MappedKRASBarcodes-XR.csv'), row.names = FALSE, quote = FALSE)





#Load TCellCluster Identities#######
b1_cid = read.csv('./code_inputs/2024-01-03_batch1-TCluster-Assignments.csv')
b1_cid$Barcode = substr(b1_cid$Barcode, 1, 16)
b2_cid = read.csv('./code_inputs/2024-01-03_batch2-TCluster-Assignments.csv')
b2_cid$Barcode = substr(b2_cid$Barcode, 1, 16)
#Re-append to Seurat objects 
m1$TCluster = b1_cid$TCluster
m2$TCluster = b2_cid$TCluster

idx = match(rownames(m1), batch_annot$barcode)
m1$batch_annotation = batch_annot$annot[idx]
idx = match(rownames(m2), batch_annot$barcode)
m2$batch_annotation = batch_annot$annot[idx]
seu1@meta.data = m1
seu2@meta.data = m2

#save new seurat objects 
save(seu1, file = paste0('./code_outputs/',Sys.Date(), '_Batch1TCells-KRASAnnot.rda'))
save(seu2, file = paste0('./code_outputs/',Sys.Date(), '_Batch2TCells-KRASAnnot.rda'))
write.csv(m1, file =paste0('./code_outputs/',Sys.Date(),  '_Batch1TCells-KRASAnnot-Metadata.csv'), quote = FALSE, row.names = TRUE)
write.csv(m2, file = paste0('./code_outputs/',Sys.Date(), '_Batch2TCells-KRASAnnot-Metadata.csv'), quote = FALSE, row.names = TRUE)
#Plot Summary Mapped Data------------------

spl = split(summary, summary$Subject)

cols = colorRampPalette(brewer.pal(8, 'Set2'))(12)
names(cols) = c('CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 
                'CD4 TCM', 'CD4 TEM', 'CD8 Naive', 'CD8 Proliferating', 'CD8 TCM', 
                'CD8 TEM', 'dnT', 'MAIT', 'Treg')
plots = lapply(spl, FUN = function(x){
  x = x[!is.na(x$Timepoint),]
  gg = ggplot(x, aes(x = Reactivity, y = Freq, fill = celltype)) + 
    geom_col() + 
    facet_wrap(~Timepoint, scale = 'free_y') + 
    theme_prism() + 
    coord_flip() + 
    scale_fill_manual(values = cols)
    
})

for(p in 1:5){
  fname = paste0('./code_outputs/', Sys.Date(), '_',subjects[p],
                 '_mappedExpansionFreqs-XR.pdf')
  pdf(file = fname, width = 8, height = 3.1)
  plot(plots[[p]])
  dev.off()
  
  fname = paste0('./code_outputs/', Sys.Date(), '_',subjects[p],
                 '_mappedExpansionFreqs-XR.png')
  
  agg_png(filename = fname, width = 4000, height = 1550, scaling = 6)
  plot(plots[[p]])
  dev.off()
  
}
write.csv(summary, file = './code_outputs/MappedExpansionFreqs-XR.csv', row.names = FALSE)
#Now distribute XR Sequences across respective antigen categories. ######
barcodes_spl_xrdist = lapply(barcodes_spl, FUN = function(x){
  xr =x[x$reactivity == 'XR', ]
  sr =x[x$reactivity != 'XR', ]
  
  xr.dup = vector(mode= 'list', length = nrow(xr))
  for( i in 1:nrow(xr)){
    tmp = xr[i,]
    ags = strsplit(xr$all.antigens[i], split = ';')
    df = data.frame(
      Subject = tmp$Subject, 
      Sample = tmp$Sample, 
      Barcode = tmp$Barcode, 
      uid = tmp$uid, 
      cdr3 = tmp$cdr3, 
      trbv.10x = tmp$trbv.10x, 
      batch_annotation = tmp$batch_annotation, 
      is_cell = tmp$is_cell, 
      reactivity = ags, 
      all.antigens = tmp$all.antigens, 
      adaptive.tcrb.id = tmp$adaptive.tcrb.id, 
      adaptive.proc.id = tmp$adaptive.proc.id, 
      adaptive.trbv = tmp$adaptive.trbv, 
      PASS = tmp$PASS, 
      CHECK = tmp$CHECK
    )
    colnames(df)[9] = 'reactivity'
    rownames(df) = NULL
    xr.dup[[i]] = df
  }
  
  xr.dup = bind_rows(xr.dup)
  x.new = rbind(xr.dup, sr)
  return(x.new)
})


ags = unique(bind_rows(barcodes_spl_xrdist)$reactivity)
ags = ags[order(ags)]

summary = vector(mode = 'list', length = length(subjects))

for( s in 1:length(subjects)){
  summary[[s]] = SummarizeMapped(subjects[s], b = barcodes_spl_xrdist[[s]])
}
summary = bind_rows(summary)
colnames(summary) = c('celltype', 'Sample', 'Reactivity' , 'Timepoint', 'Freq')
idx = match(summary$Sample, meta$Sample.Name)
summary$Subject = meta$Subject.ID[idx]



spl = split(summary, summary$Subject)

cols = colorRampPalette(brewer.pal(8, 'Set2'))(12)
names(cols) = c('CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 
                'CD4 TCM', 'CD4 TEM', 'CD8 Naive', 'CD8 Proliferating', 'CD8 TCM', 
                'CD8 TEM', 'dnT', 'MAIT', 'Treg')
plots = lapply(spl, FUN = function(x){
  x = x[!is.na(x$Timepoint),]
  gg = ggplot(x, aes(x = Reactivity, y = Freq, fill = celltype)) + 
    geom_col() + 
    facet_wrap(~Timepoint, scale = 'free_y') + 
    theme_prism() + 
    coord_flip() + 
    scale_fill_manual(values = cols)
  
})

for(p in 1:5){
  fname = paste0('./code_outputs/', Sys.Date(), '_',subjects[p],
                 '_mappedExpansionFreqs-noXR.pdf')
  pdf(file = fname, width = 8, height = 3.1)
  plot(plots[[p]])
  dev.off()
  
  fname = paste0('./code_outputs/', Sys.Date(), '_',subjects[p],
                 '_mappedExpansionFreqs-noXR.png')
  
  agg_png(filename = fname, width = 4000, height = 1550, scaling = 6)
  plot(plots[[p]])
  dev.off()
  
}
write.csv(summary, file = './code_outputs/MappedExpansionFreqs-noXR.csv', row.names = FALSE)
#Now filter by Pre timepoint. Include XR ##################
barcodes_spl_preFilt = lapply(barcodes_spl, FUN = function(x){
  idx = match(x$Sample, meta$Sample.Name)
  x$Timepoint = meta$Timepoint[idx]
  in_pre = unique(x$uid[x$Timepoint == 'Pre'])
  x.filt = x[!(x$uid %in% in_pre),]
  return(x.filt)
})



ags = unique(bind_rows(barcodes_spl_preFilt)$reactivity)
ags = ags[order(ags)]

summary = vector(mode = 'list', length = length(subjects))

for( s in 1:length(subjects)){
  summary[[s]] = SummarizeMapped(subjects[s], b = barcodes_spl_preFilt[[s]])
}
summary = bind_rows(summary)
colnames(summary) = c('celltype', 'Sample', 'Reactivity' , 'Timepoint', 'Freq')
idx = match(summary$Sample, meta$Sample.Name)
summary$Subject = meta$Subject.ID[idx]



spl = split(summary, summary$Subject)

cols = colorRampPalette(brewer.pal(8, 'Set2'))(12)
names(cols) = c('CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 
                'CD4 TCM', 'CD4 TEM', 'CD8 Naive', 'CD8 Proliferating', 'CD8 TCM', 
                'CD8 TEM', 'dnT', 'MAIT', 'Treg')
plots = lapply(spl, FUN = function(x){
  x = x[!is.na(x$Timepoint),]
  gg = ggplot(x, aes(x = Reactivity, y = Freq, fill = celltype)) + 
    geom_col() + 
    #facet_wrap(~Timepoint, scale = 'free_y') + 
    theme_prism() + 
    coord_flip() + 
    scale_fill_manual(values = cols)
  
})

for(p in 1:5){
  fname = paste0('./code_outputs/', Sys.Date(), '_',subjects[p],
                 '_mappedExpansionFreqs-XR-preFilt.pdf')
  pdf(file = fname, width = 4.5, height = 3.1)
  plot(plots[[p]])
  dev.off()
  
  fname = paste0('./code_outputs/', Sys.Date(), '_',subjects[p],
                 '_mappedExpansionFreqs-XR-preFilt.png')
  
  agg_png(filename = fname, width = 2500, height = 1550, scaling = 6)
  plot(plots[[p]])
  dev.off()
  
}
write.csv(summary, file = './code_outputs/MappedExpansionFreqs-XR-preFilt.csv', row.names = FALSE)
#Now filter by Pre timepoint. Distribute XR ##################
barcodes_spl_xrdist_preFilt = lapply(barcodes_spl_xrdist, FUN = function(x){
  idx = match(x$Sample, meta$Sample.Name)
  x$Timepoint = meta$Timepoint[idx]
  in_pre = unique(x$uid[x$Timepoint == 'Pre'])
  x.filt = x[!(x$uid %in% in_pre),]
  return(x.filt)
})



ags = unique(bind_rows(barcodes_spl_xrdist_preFilt)$reactivity)
ags = ags[order(ags)]

summary = vector(mode = 'list', length = length(subjects))

for( s in 1:length(subjects)){
  summary[[s]] = SummarizeMapped(subjects[s], b = barcodes_spl_xrdist_preFilt[[s]])
}
summary = bind_rows(summary)
colnames(summary) = c('celltype', 'Sample', 'Reactivity' , 'Timepoint', 'Freq')
idx = match(summary$Sample, meta$Sample.Name)
summary$Subject = meta$Subject.ID[idx]



spl = split(summary, summary$Subject)

cols = colorRampPalette(brewer.pal(8, 'Set2'))(12)
names(cols) = c('CD4 CTL', 'CD4 Naive', 'CD4 Proliferating', 
                'CD4 TCM', 'CD4 TEM', 'CD8 Naive', 'CD8 Proliferating', 'CD8 TCM', 
                'CD8 TEM', 'dnT', 'MAIT', 'Treg')
plots = lapply(spl, FUN = function(x){
  x = x[!is.na(x$Timepoint),]
  gg = ggplot(x, aes(x = Reactivity, y = Freq, fill = celltype)) + 
    geom_col() + 
    #facet_wrap(~Timepoint, scale = 'free_y') + 
    theme_prism() + 
    coord_flip() + 
    scale_fill_manual(values = cols)
  
})

for(p in 1:5){
  fname = paste0('./code_outputs/', Sys.Date(), '_',subjects[p],
                 '_mappedExpansionFreqs-noXR-preFilt.pdf')
  pdf(file = fname, width = 4.5, height = 3.1)
  plot(plots[[p]])
  dev.off()
  
  fname = paste0('./code_outputs/', Sys.Date(), '_',subjects[p],
                 '_mappedExpansionFreqs-noXR-preFilt.png')
  
  agg_png(filename = fname, width = 2500, height = 1550, scaling = 6)
  plot(plots[[p]])
  dev.off()
  
}

write.csv(summary, file = './code_outputs/MappedExpansionFreqs-noXR-preFilt.csv', row.names = FALSE)