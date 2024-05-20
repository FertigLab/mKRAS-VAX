#produce table of filtered diffEx data meeting specified criteria. 
#Append raw expansion data for reference. 
#AAG 18 April 2024 

library(dplyr)
library(xlsx)

setwd("~/salvage_tmp/final/")
source('./_rfunctions/ProcessAdaptiveFile.R')

MIN_OR = 5; MIN_FREQ = 0.00
pts = c('01', '02', '04', '08', '13')
diffEx_all = lapply(
  paste0('./raw_data/diffex-data/2023-06-07_J1994-', pts,'_diffEx.tsv'), 
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

diffEx_all_complete = bind_rows(diffEx_all, .id = 'Subject')
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
  
  x$uid = paste0(ProcessAdaptiveVgenes(x$TRBV), '_', x$CDR3.beta.aa)
  x$adaptive.tcrb.id = paste0(x$TRBV, '_', x$CDR3.beta.aa)
  x$Antigen[x$XR] =de$xrAgs[de$xr]
  x$x[x$XR] = 'XR'
  return(x)
})


diffex= bind_rows(diffex, .id = 'Subject')

m = matrix(data = NA, ncol = 12, nrow = nrow(diffex))
colnames(m) = paste0(rep(c('G12A', 'G12C', 'G12D', 'G12R', 'G12V', 'G13D'),each = 2), c('_OR', '_F1'))
tmp = gsub('\\*', '9', diffEx_all_complete$id)
for(r in 1:nrow(diffex)){
  x = diffex[r,]
  cdr3beta = gsub('\\*', '9', x$adaptive.tcrb.id)
  df = diffEx_all_complete[c(1:nrow(diffEx_all_complete) %in% grep(cdr3beta, tmp)),]
  
  #if df contains expansion info for multiple patients, filter to that patient. 
  df = df[df$Subject == x$Subject,]
  
  if(length(grep('G12A', df$Antigen)) > 0){
    y = df[which(df$Antigen == 'G12A'),]
    if(y$f2 == 0) y$or = Inf
    m[r, 1] = y$or
    m[r,2] = y$f1
  }
  if(length(grep('G12C', df$Antigen)) > 0){
    y = df[which(df$Antigen == 'G12C'),]
    if(y$f2 == 0) y$or = Inf
    m[r, 3] = y$or
    m[r,4] = y$f1
  }
  if(length(grep('G12D', df$Antigen)) > 0){
    y = df[which(df$Antigen == 'G12D'),]
    if(y$f2 == 0) y$or = Inf
    m[r, 5] = y$or
    m[r,6] = y$f1
  }
  if(length(grep('G12R', df$Antigen)) > 0){
    y = df[which(df$Antigen == 'G12R'),]
    if(y$f2 == 0) y$or = Inf
    m[r, 7] = y$or
    m[r,8] = y$f1
  }
  if(length(grep('G12V', df$Antigen)) > 0){
    y = df[which(df$Antigen == 'G12V'),]
    if(y$f2 == 0) y$or = Inf
    m[r, 9] = y$or
    m[r,10] = y$f1
  }
  if(length(grep('G13D', df$Antigen)) > 0){
    y = df[which(df$Antigen == 'G13D'),]
    if(y$f2 == 0) y$or = Inf
    m[r, 11] = y$or
    m[r,12] = y$f1
  }
}

m = as.data.frame.matrix(m)
final_table = cbind(diffex, m)
pub.key = data.frame(
  orig = c('J1994.001', 'J1994.002', 'J1994.004', 'J1994.008', 'J1994.013'), 
  new = c('J1994_12', 'J1994_5', 'J1994_8', 'J1994_2', 'J1994_6')
)
idx = match(final_table$Subject, pub.key$orig)
final_table$Subject = pub.key$new[idx]
write.table(final_table, file = paste0('./code_outputs/',Sys.Date(), '_all-filtered-diffex.tsv'), sep = '\t')
spl = split(final_table, final_table$Subject)
for(i in 1:length(spl)){
  write.xlsx(spl[[i]], file = paste0('./code_outputs/',Sys.Date(), '_All-Expansion-TCRB.xlsx'), row.names = FALSE,
             sheetName = names(spl)[i], append = TRUE)
}


