#Run diffEx testing on AH Expansion data under a variety of parameters. 
#Determine optimal cutoffs to call clones 'expanded'. 
#AAG 5 May 2023 

#Order of operations: 
#For each patient, 
# For each expansion condition, 
#   Calculate diffex TCRs against Baseline (a) and against No Peptide (b) 
#   Display Odds Ratio vs. P value 
library(parallel)
setwd("~/salvage_tmp/final/")
diffExTest <- function(name1,name2, alpha = 0.05, filter = TRUE){
  require(mutoss)
  
  d1 = read.delim(name1)
  d2 = read.delim(name2)
  
  s1 = sum(d1$count)
  s2 = sum(d2$count)
  
  ids = unique(c(d1$id, d2$id)) 
  
  df = data.frame(id = ids, c1 = 0, c2 = 0, pval = NA, or = NA)
  
  df$c1[match(d1$id, df$id)] = d1$count
  df$c2[match(d2$id, df$id)] = d2$count
  df$f1 = df$c1/s1
  df$f2 = df$c2/s2
  
  #Filtering TCRs which appear no more than once in both datasets. 
  torm = which(df$c1 < 2 & df$c2 < 2)
  if(length(torm) > 0)   df = df[-torm, ]

  
  for(i  in 1:nrow(df)){
    ctab = matrix(nrow = 2, ncol = 2)
    ctab[1,1] = df$c1[i]
    ctab[1,2] = df$c2[i]
    ctab[2,1] = s1 - df$c1[i]
    ctab[2,2] = s2 - df$c2[i]
    tmp = fisher.test(ctab, alternative = 'greater')
    df$pval[i] = tmp$p.value
    df$or[i] = tmp$estimate[[1]][1]
  }
  
  require(mutoss)
  df$padj = BY(df$pval, alpha = 0.05, silent = TRUE)$adjPValues
  df = df[order(df$pval),]
  row.names(df) = NULL
  
  if(filter){
    df = df[which(df$padj <= alpha),]
  }
  return(df)
}


getFreqThreshold = function(n, p)
{
  #Copied from Luda source code 
  #Determine minimum count threshold according to 
  # calculate frequency threshold using the number of cells and probability
  # (1-((1-p)^(1/n)))
  # where n=number of cells per well and p=selected probability
  return (1-((1-p)^(1/n)))
}

diffExProcessor <- function(patient_id, data_dir = './code_inputs/expansions_agg/' ){
  files = dir(data_dir)[grep(patient_id, dir(data_dir))]
  files = paste0(data_dir, files)
  KRASCond = files[grep('G12|G13', files)] #Experimental Stimulation Conditions.
  nopep = files[grep('peptide', files)]
  baseline = files[grep('Baseline|baseline', files)]
  
  #Ignore baseline and compile results against no peptide condition. 
 noPep_res = lapply(KRASCond, FUN = diffExTest, name2 = nopep)
 names(noPep_res) = c('G12A', 'G12C', 'G12D', 'G12R', 'G12V', 'G13D') #A bit risky! Assumes static filenames
 res = bind_rows(noPep_res, .id = 'Antigen')

 #Ignore no peptide and compile results against baseline. 
  baseline_res = lapply(KRASCond, FUN = diffExTest, name2 = baseline)
  names(baseline_res)= c('G12A', 'G12C', 'G12D', 'G12R', 'G12V', 'G13D') #A bit risky! Assumes static filenames
  baseline_res = bind_rows(baseline_res, .id = 'Antigen')
  
 #Determine count threshold
 baseline_dat = read.delim(baseline)
 min_freq = getFreqThreshold(n = 7e5, p = .01) #Most of the time, this numeber is < 1 count
 baseline_dat = baseline_dat[which(baseline_dat$freq > min_freq), ]

 res$baseline_thresh = res$id %in% baseline_dat$id
 res$or_thresh = res$or >= 5
 res$count_thresh = res$c1 >= 10

 xr  =  names(table(res$id)[which(table(res$id) > 1)])
 xr_or = names(table(res$id[which(res$or_thresh)])[which(table(res$id[which(res$or_thresh)]) > 1)])
 xr_baseline =  names(table(res$id[which(res$baseline_thresh)])[which(table(res$id[which(res$baseline_thresh)]) > 1)])
 xr_all = xr_baseline[which(xr_baseline %in% xr_or)]

 res$xr = res$id %in% xr
 res$xr_or = res$id %in% xr_or
 res$xr_baseline = res$id %in% xr_baseline
 res$xr_all = res$id %in% xr_all
 
 write.table(res, row.names = FALSE, file = paste0(Sys.Date(), '_', 
                                patient_id, '_diffEx.tsv'))
 
 df = as.data.frame.table(table(res$Antigen))
 
 
 df2 = data.frame(
   Var1 = c('nTot', 'nXR', 'nXR_baseline', 'nXR_OR', 'nXR_all'), 
   Freq = c(length(unique(res$id)), length(xr),
            length(xr_baseline),
            length(xr_or),
            length(xr_all))
 )

 
   
 df = rbind(df, df2)
 return(df)
}

patient_ids = paste0('J1994-', c('01', '02', '04', '08', '13'))

system.time(
out  <- mclapply(patient_ids, diffExProcessor, mc.cores = 5)
)
names(out) = patient_ids
out = bind_rows(out, .id = 'Patient')


