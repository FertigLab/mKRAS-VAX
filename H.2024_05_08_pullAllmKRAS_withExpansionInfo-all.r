#Generate holistic list meeting indicated OR and Freq criteria 
#displaying all mapped anti-KRAS TCRs found in single cell data, 
#along with their compiled sc annotation information and 
#differential expression statistics. 
#AAG 12 April 2024

######Set Up Workspace########
library(parallel)
library(dplyr)
library(xlsx)
#setwd("/media/aag7319/WDRed/ZZZ_Homolig/data_exports/Executable-Code/DiffEx-TCR-Mapping/")
setwd("~/salvage_tmp/final/")
patient_ids = paste0('J1994-', c('01', '02', '04', '08', '13'))
MIN_OR = 5; MIN_FREQ = 0.00; MIN.OR=5;MIN.FREQ=0.00 #unfortunately both variable names are used through script. Define both here.
####STAGE ONE: Get all single cell annotations for mKRAS TCR beta chains found in scVDJ data################
#List of annotated scVDJ files for each patient. 
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

#files = paste0('./test_data/NZ01JHU506_000_analysis/cellranger/vdj/', samples, '/outs/airr_rearrangement.tsv')
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
names(patient_scVDJ_files) = patient_ids
#Load All DiffEx Data----------------------------------------------------------
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


#Clean up workspace------
agg = bind_rows(diffEx_all, .id = 'Subject')
agg_filt = agg[(agg$or >= MIN_OR & agg$f1 >= MIN_FREQ), ]
tab = table(agg_filt$id, agg_filt$Subject) > 0
public_clones = rownames(tab)[which(rowSums(tab) > 1)]
public_individuals = apply(tab[which(rowSums(tab) > 1),], 1, FUN = function(x){
  paste0(names(x[which(x > 0)]), collapse = ';')
})

diffEx_all = lapply(diffEx_all, FUN = function(de){
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
  return(de)
})

#For each patient, assign consensus cell phenotype-------
idx = match(names(scVDJ_annot), meta$Sample.Name)
sc_by_subject = split(scVDJ_annot, meta$Subject.ID[idx]);
sc_by_subject = lapply(sc_by_subject, FUN = bind_rows, .id = 'Sample')
process = function(de, scdat){
  tab = table(de$id[de$or >= MIN_OR & de$f1 >= MIN_FREQ], de$Antigen[de$or >= MIN_OR & de$f1 >= MIN_FREQ])
  xr_ids = rownames(tab)[rowSums(tab) > 1]
  de$cdr3 = gsub('.*_','',de$id)
  de$IN_SC = de$cdr3 %in% scdat$junction_aa
  
  
  alpha_chains = scdat[grep('TRAV', scdat$v_call),]
  alpha_uids = paste0(alpha_chains$Sample, alpha_chains$cell_id, alpha_chains$clone_id)
  insc = de[de$IN_SC,]
  #insc = de #Do NOT filter whether you are in sc data or not. 
  out = lapply(insc$cdr3, FUN = function(x){
    matched = scdat[which(scdat$junction_aa == x),]
    nSamp = length(unique(matched$Sample))
    Samps = paste0(unique(matched$Sample), collapse = ';')
    annot = table(matched$batch_annotation)
    if(length(which(annot == max(annot))) == 1){
      consensus_annot = names(annot)[which(annot == max(annot))]
    }else{
      consensus_annot = NA
    }
    if(length(annot) >= 1){
      full_annot = paste0(names(annot), '(', as.character(annot), ')', collapse = ';')
    }else{
      full_annot = NA
    }
    tmp = table(matched$cell_id, matched$Sample) > 0
    nCells = colSums(tmp)
    totCells= length(unique(matched$cell_id))
    nCells = paste0(names(nCells), '(', as.character(nCells), ')', collapse = ';')
    cellids = paste0(matched$Sample, matched$cell_id, matched$clone_id)
    has_alpha = any(cellids %in%  alpha_uids)
    
    df = data.frame(
      consensus_annot = consensus_annot, 
      nSamp = nSamp, 
      Samples = Samps,
      nCells = totCells,
      nCells_by_sample = nCells, 
      has_alpha = has_alpha, 
      annot = full_annot
    )
    return(df)
  })
  out = bind_rows(out)
  
  de_annot = cbind(insc, out)
  
  
  notsc = de[!de$IN_SC,]
  df = data.frame(
    consensus_annot = rep(NA, nrow(notsc)),
    nSamp = rep(NA, nrow(notsc)), 
    Samples = rep(NA, nrow(notsc)),
    nCells = rep(NA, nrow(notsc)),
    nCells_by_sample = rep(NA, nrow(notsc)), 
    has_alpha = rep(NA, nrow(notsc)),
    annot = rep(NA, nrow(notsc))
  )
  
  notannot = cbind(notsc, df)
  
  
  de_annot = rbind(de_annot, notannot)
  
  nDiffEx = length(unique(de$id))
  nxr = length(xr_ids)
  nMapped = length(unique(insc$id))
  nxr_mapped = sum(unique(insc$id) %in% xr_ids)
  # de_annot$has_alpha = NULL; de_annot$IN_SC = NULL; de_annot$cdr3 = NULL
  hasPaired = length(unique(de_annot$id))
  nxr_paired = sum(unique(de_annot$id %in% xr_ids))
  
  summary = data.frame(
    nDiffEx = nDiffEx, 
    nMapped = nMapped, 
    nPaired = hasPaired, 
    nxr = nxr, 
    nxr_mapped = nxr_mapped, 
    nxr_paired = nxr_paired
  )
  de_annot = de_annot[order(de_annot$Antigen, de_annot$or, de_annot$f1, decreasing = TRUE),]
  de_annot = split(de_annot, de_annot$xr)
  sr = de_annot[[1]]
  sr$xr = NULL
  sr$xrAgs = NULL
  single_reactive = split(sr, sr$Antigen)
  cross_reactive = de_annot[[2]]
  cross_reactive$xr = NULL
  
  return(list(single_reactive = single_reactive, cross_reactive = cross_reactive, summary = summary))
}

results = lapply(diffEx_all, FUN = process, scdat = bind_rows(sc_by_subject))
results = c(results[[1]], results[[2]], results[[3]], results[[4]], results[[5]])
#results = mapply(process, de = diffEx_all, scdat = sc_by_subject)

summaries = results[c(3,6,9,12,15)]
names(summaries) = names(diffEx_all)
summaries = bind_rows(summaries, .id = 'Subject')
single_reactives= results[c(1,4,7,10,13)]; names(single_reactives) = names(diffEx_all)
cross_reactives = results[c(2,5,8,11,14)]; names(cross_reactives) = names(diffEx_all)


refineXR <- function(x){
  #x = x[x$IN_SC & x$or >= MIN_OR & x$f1 >= MIN_FREQ & x$has_alpha,]
  x = x[x$IN_SC & x$or >= MIN_OR & x$f1 >= MIN_FREQ,]
  x = x[match(unique(x$id), x$id),]
  x = x[,colnames(x) %in% c('id', 'xrAgs', 'cdr3', 'consensus_annot', 
                            'nSamp', 'nCells', 'Samples', 'nCells_by_sample', 
                            'annot','has_alpha')]
  x$XR = TRUE
  colnames(x)[colnames(x) == 'xrAgs'] = 'Antigen'
  return(x)
}

refineSR <- function(x){
  x = bind_rows(x)
 # x = x[x$IN_SC & x$or >= MIN_OR & x$f1 >= MIN_FREQ & x$has_alpha,]
  x = x[x$IN_SC & x$or >= MIN_OR & x$f1 >= MIN_FREQ,]
  x = x[,colnames(x) %in% c('id', 'Antigen', 'cdr3', 'consensus_annot', 
                            'nSamp', 'nCells', 'Samples', 'nCells_by_sample', 
                            'annot', 'has_alpha')]
  x$XR = FALSE
  return(x)
}
xrfilt = lapply(cross_reactives, refineXR)
srfilt = lapply(single_reactives, refineSR)

srfilt = bind_rows(srfilt, .id = 'subject')
xrfilt = bind_rows(xrfilt, .id = 'subject')

allfilt1 = rbind(srfilt, xrfilt)


##Do again but only permit single cell from that patient to be matched to expansion
results = mapply(process, de = diffEx_all, scdat = sc_by_subject)

summaries = results[c(3,6,9,12,15)]
names(summaries) = names(diffEx_all)
summaries = bind_rows(summaries, .id = 'Subject')
single_reactives= results[c(1,4,7,10,13)]; names(single_reactives) = names(diffEx_all)
cross_reactives = results[c(2,5,8,11,14)]; names(cross_reactives) = names(diffEx_all)


refineXR <- function(x){
  #x = x[x$IN_SC & x$or >= MIN_OR & x$f1 >= MIN_FREQ & x$has_alpha,]
  x = x[x$IN_SC & x$or >= MIN_OR & x$f1 >= MIN_FREQ,]
  x = x[match(unique(x$id), x$id),]
  x = x[,colnames(x) %in% c('id', 'xrAgs', 'cdr3', 'consensus_annot', 
                            'nSamp', 'nCells', 'Samples', 'nCells_by_sample', 
                            'annot','has_alpha')]
  x$XR = TRUE
  colnames(x)[colnames(x) == 'xrAgs'] = 'Antigen'
  return(x)
}

refineSR <- function(x){
  x = bind_rows(x)
  #x = x[x$IN_SC & x$or >= MIN_OR & x$f1 >= MIN_FREQ & x$has_alpha,]
  x = x[x$IN_SC & x$or >= MIN_OR & x$f1 >= MIN_FREQ,]
  x = x[,colnames(x) %in% c('id', 'Antigen', 'cdr3', 'consensus_annot', 
                            'nSamp', 'nCells', 'Samples', 'nCells_by_sample', 
                            'annot','has_alpha')]
  x$XR = FALSE
  return(x)
}
xrfilt = lapply(cross_reactives, refineXR)
srfilt = lapply(single_reactives, refineSR)

srfilt = bind_rows(srfilt, .id = 'subject') 
xrfilt = bind_rows(xrfilt, .id = 'subject')
allfilt2 = rbind(srfilt, xrfilt)


#Now identify ones in allfilt1 which are not ided in allfilt2
allfilt1$Endogeneous.Mapping = allfilt1$id %in% allfilt2$id
allfilt = allfilt1
meta1 = meta
####STAGE TWO: Obtain Paired-Chain info and diffEx statistics########

meta = data.frame(
  patient = paste0('J1994.0', c('01', '01', '01', '01', '02', '02', 
                                '04', '04', '08', '08', '13', '13')), 
  sample = c('C1D1T', 'C2D1T', 'C3D1T', 'C6D1T', 
             paste0('AH_', c(1:8), 'T_B', c(1:8)))
)

scVDJ = bind_rows(scVDJ_annot, .id = 'Sample')
scVDJ$Patient = meta$patient[match(scVDJ$Sample, meta$sample)]

#Make list of all patient TCRs to pull. 


toGet = diffEx_all

filtered = diffEx_all
filtered = lapply(filtered, FUN = function(x){
  x = x[x$or >= MIN.OR & x$f1 >= MIN.FREQ,]
  return(x)
})
toGet = filtered

#For each patient set of TCRs, collapse to unique V/beta chains (there may be duplicates)
toGet = lapply(toGet, FUN = function(x){
  colnames(x)[1] = 'Antigen'
  x = x[match(unique(x$id), x$id),]
  x$cdr3_aa = gsub('.*_', '', x$id)
  x$v_call = gsub('_.*', '', x$id)
  return(x)
})
#For each patient unique CDR3, go through w process function. 
process <- function(dat, patient = 'J1994.001'){
  #First, confirm all beta chains can be matched within that patient. 
  #If not,  there has been some technical error prior to this. 
  ALL_PRESENT = all(dat$cdr3_aa %in% dat$cdr3_aa)
  if(!ALL_PRESENT){print('Error: Not all beta chains present.'); return(1)}
  
  #Get unique cell ids and pull alpha info. 
  #Make consensus alpha info if multiple TCRalphas present with matching Beta aa. 
  if(is.na(patient)){
    scVDJ = scVDJ[scVDJ$Patient != 'J1994.004',]
  }else{
    scVDJ = scVDJ[scVDJ$Patient == patient, ]  
  }
  
  alpha_chains = scVDJ[grep('TRAV', scVDJ$v_call), ] #all alpha chains from this patient. 
  beta_chains = scVDJ[grep('TRBV', scVDJ$v_call), ] #all beta chains from this patient. 
  matches = beta_chains[beta_chains$junction_aa %in% dat$cdr3_aa,]
  
  spl = split(matches, matches$junction_aa)
  cell_ids = lapply(spl, FUN = function(x){
    return(unique(x$cell_id))}) #All cell ids containing aa beta chain match 
  
  alpha_info  = lapply(cell_ids, FUN = function(cell_ids){
    alphas = alpha_chains[ alpha_chains$cell_id %in% cell_ids, ]
    
    if(nrow(alphas) == 0){
      df = data.frame(
        TRAV = NA, 
        CDR3.alpha.aa = NA,
        TRAJ = NA, 
        multiple_alpha_calls = FALSE
      )
    }else if(nrow(alphas) == 1){
      df = data.frame(
        TRAV = alphas$v_call, 
        CDR3.alpha.aa = alphas$junction_aa,
        TRAJ = alphas$j_call, 
        multiple_alpha_calls = FALSE
      )
    }else{
      vtab = table(alphas$v_call)
      consensus_v = names(vtab)[which(vtab == max(vtab))[1]]
      CDR3.alphatab = table(alphas$junction_aa)
      consensus_cdr3 = names(CDR3.alphatab)[which(CDR3.alphatab == max(CDR3.alphatab))[1]]
      jtab = table(alphas$j_call)
      consensus_j = names(jtab)[which(jtab == max(jtab))[1]]
      
      df = data.frame(
        TRAV = consensus_v,  
        CDR3.alpha.aa = consensus_cdr3, 
        TRAJ = consensus_j, 
        multiple_alpha_calls = TRUE
      )
      
      if(length(unique(alphas$junction_aa)) == 1 & 
         length(unique(alphas$v_call)) == 1 & 
         length(unique(alphas$j_call)) == 1){
        df$multiple_alpha_calls = FALSE
      }
    }
    
    return(df)
    
  })
  
  alpha_info = bind_rows(alpha_info, .id = 'CDR3.beta.aa')
  
  beta_info = lapply(cell_ids, FUN = function(cell_ids){
    betas = beta_chains[ (beta_chains$cell_id %in% cell_ids & beta_chains$junction_aa %in% dat$cdr3_aa), ]
    
    if(nrow(betas) == 1){
      df = data.frame(
        TRBV = betas$v_call, 
        CDR3.beta.aa = betas$junction_aa,
        TRBJ = betas$j_call, 
        multiple_beta_calls = FALSE
      )
    }else{
      vtab = table(betas$v_call)
      consensus_v = names(vtab)[which(vtab == max(vtab))[1]]
      CDR3.betatab = table(betas$junction_aa)
      consensus_cdr3 = names(CDR3.betatab)[which(CDR3.betatab == max(CDR3.betatab))[1]]
      jtab = table(betas$j_call)
      consensus_j = names(jtab)[which(jtab == max(jtab))[1]]
      
      df = data.frame(
        TRBV = consensus_v,  
        CDR3.beta.aa = consensus_cdr3, 
        TRBJ = consensus_j, 
        multiple_beta_calls = TRUE
      )
    }
    
    return(df)
  })
  
  beta_info = bind_rows(beta_info)
  idx = match(beta_info$CDR3.beta.aa, alpha_info$CDR3.beta.aa)
  alpha_info$CDR3.beta.aa = NULL
  combo_info = cbind(beta_info, alpha_info)
  
 #combo_info = combo_info[!is.na(combo_info$CDR3.alpha.aa),] #Filter to paired chain 
  
  return(combo_info)
}


all_tcrs = list(
  J1994.001 =  process(toGet$J1994.001, patient = NA),
  J1994.002 = process(toGet$J1994.002, patient = NA),
  J1994.004 = process(toGet$J1994.004, patient = NA),
  J1994.008 = process(toGet$J1994.008, patient = NA),
  J1994.013 = process(toGet$J1994.013, patient = NA)
)

#all_tcrs = all_tcrs[names(all_tcrs) != 'J1994.004']

all_tcrs = bind_rows(all_tcrs, .id = 'Patient')
#rm(list = setdiff(ls(), c('all_tcrs', 'allfilt', 'meta')))

#Discrepancies in membership between all_tcrs and allfilt are due to the 
#inclusion/exclusion of J1994.004 at various stages of processing. 
#If one wants to totally omit subject 4 from output, only consider 
#the intersection of these lists. 
all_cdr3 = unique(c(all_tcrs$CDR3.beta.aa, allfilt$cdr3))
table(all_cdr3 %in% all_tcrs$CDR3.beta.aa, all_cdr3 %in% allfilt$cdr3)


case2 = all_cdr3[!(all_cdr3 %in% all_tcrs$CDR3.beta.aa) & (all_cdr3 %in% allfilt$cdr3)]

toGet2 = bind_rows(toGet, .id ='subject' )
toGet2 = toGet2[toGet2$cdr3_aa %in% case2,]
toGet2[toGet2$subject != 'J1994.004',]
supp = process(toGet2, patient = NA)

both = all_tcrs[all_tcrs$CDR3.beta.aa %in% allfilt$cdr3,]
idx = match(both$CDR3.beta.aa, allfilt$cdr3)
both = cbind(both, allfilt[idx, colnames(allfilt) %in% c('subject', 'Antigen', 'id', 
                                                         'consensus_annot', 'nSamp', 'Samples', 'nCells', 'nCells_by_sample', 'annot', 'XR', 'Endogenous.Mapping')])

pre_samples = paste0(meta1$Sample.Name[meta1$Timepoint == 'Pre'][1:5], collapse = '|')

both$Pre.Tx.Timepoint = FALSE
both$Pre.Tx.Timepoint[grep(pre_samples, both$Samples)] = TRUE
both$subject = NULL

colnames(both)[colnames(both) == 'id'] = 'adaptive.tcrb.id'
final_table = both[both$Patient != 'J1994.004',]

#Obtain exact diffex data for all TCRs. 
diffex= bind_rows(diffEx_all, .id = 'Subject')

m = matrix(data = NA, ncol = 12, nrow = nrow(final_table))
colnames(m) = paste0(rep(c('G12A', 'G12C', 'G12D', 'G12R', 'G12V', 'G13D'),each = 2), c('_OR', '_F1'))
tmp = gsub('\\*', '9', diffex$id)
for(r in 1:nrow(final_table)){
  x = final_table[r,]
  cdr3beta = gsub('\\*', '9', x$adaptive.tcrb.id)
  df = diffex[c(1:nrow(diffex) %in% grep(cdr3beta, tmp)),]
  
  #if df contains expansion info for multiple patients, filter to that patient. 
  df = df[df$Subject == x$Patient,]
  
  if(any(df$Antigen == 'G12A')){
    y = df[which(df$Antigen == 'G12A'),]
    if(y$f2 == 0) y$or = Inf
    m[r, 1] = y$or
    m[r,2] = y$f1
  }
  if(any(df$Antigen == 'G12C')){
    y = df[which(df$Antigen == 'G12C'),]
    if(y$f2 == 0) y$or = Inf
    m[r, 3] = y$or
    m[r,4] = y$f1
  }
  if(any(df$Antigen == 'G12D')){
    y = df[which(df$Antigen == 'G12D'),]
    if(y$f2 == 0) y$or = Inf
    m[r, 5] = y$or
    m[r,6] = y$f1
  }
  if(any(df$Antigen == 'G12R')){
    y = df[which(df$Antigen == 'G12R'),]
    if(y$f2 == 0) y$or = Inf
    m[r, 7] = y$or
    m[r,8] = y$f1
  }
  if(any(df$Antigen == 'G12V')){
    y = df[which(df$Antigen == 'G12V'),]
    if(y$f2 == 0) y$or = Inf
    m[r, 9] = y$or
    m[r,10] = y$f1
  }
  if(any(df$Antigen == 'G13D')){
    y = df[which(df$Antigen == 'G13D'),]
    if(y$f2 == 0) y$or = Inf
    m[r, 11] = y$or
    m[r,12] = y$f1
  }
}

m = as.data.frame.matrix(m)
final_table = cbind(final_table, m)

write.table(final_table, file = paste0('./code_outputs/',Sys.Date(), '_mKRAS-TCRs-NoV.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)

#Swap patient ids with 'publication ids' 
pub.key = data.frame(
  orig = c('J1994.001', 'J1994.002', 'J1994.004', 'J1994.008', 'J1994.013'), 
  new = c('J1994_12', 'J1994_5', 'J1994_8', 'J1994_2', 'J1994_6')
)
#Split into big sheet and write by subject 
idx = match(final_table$Patient, pub.key$orig)
final_table$Patient = pub.key$new[idx]

spl = split(final_table, final_table$Patient)

for(i in 1:length(spl)){
  write.xlsx(spl[[i]], file = paste0('./code_outputs/',Sys.Date(), '_mKRAS-TCRs-NoV.xlsx'), row.names = FALSE,
                                     sheetName = names(spl)[i], append = TRUE)
}

           
#Optional: Consider V Match to call Single-cell mapping.  ####
source('./_rfunctions/ProcessAdaptiveFile.R')
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

vdf = data.frame(
  trbv.10x = final_table$TRBV, 
  trbv.adaptive = gsub('_.*', '', final_table$adaptive.tcrb.id)

)
vdf$failed = (ProcessAdaptiveVgenes(vdf$trbv.adaptive) != vdf$trbv.10x & 
                    !CheckTRBVEquivalence(vdf$trbv.adaptive, vdf$trbv.10x) & 
                    ProcessAdaptiveVgenes(vdf$trbv.adaptive) != ProcessAdaptiveVgenes(vdf$trbv.10x))
final_table = final_table[!vdf$failed, ]
write.table(final_table, file = paste0('./code_outputs/',Sys.Date(), '_mKRAS-TCRs-VMatch.tsv'), sep = '\t', row.names = FALSE, quote = FALSE)

spl = split(final_table, final_table$Patient)
for(i in 1:length(spl)){
  write.xlsx(spl[[i]], file = paste0('./code_outputs/',Sys.Date(), '_mKRAS-TCRs-VMatch.xlsx'), row.names = FALSE,
                                     sheetName = names(spl)[i], append = TRUE)
}




#Lastly, compute the total frequency of CD4 vs CD8 in the expansions. ####
#Use V matched table. 

exp.res = lapply(spl, FUN = function(s){
  CD4.types = c(unique(final_table$consensus_annot)[grep('CD4', unique(final_table$consensus_annot))], 'Treg')
  CD8.types = unique(final_table$consensus_annot)[grep('CD8', unique(final_table$consensus_annot))]
  
  tmp = s[s$consensus_annot %in% CD4.types,]
  cd4.df = data.frame(
    type = 'CD4', 
    mKRAS = c('G12A', 'G12C', 'G12D', 'G12R', 'G12V', 'G13D'), 
    freq = c(sum(tmp$G12A_F1, na.rm = TRUE), 
             sum(tmp$G12C_F1, na.rm = TRUE), 
             sum(tmp$G12D_F1, na.rm = TRUE), 
             sum(tmp$G12R_F1, na.rm = TRUE), 
             sum(tmp$G12V_F1, na.rm = TRUE), 
             sum(tmp$G13D_F1, na.rm = TRUE) )
  )
  
  tmp = s[s$consensus_annot %in% CD8.types,]
  cd8.df = data.frame(
    type = 'CD8', 
    mKRAS = c('G12A', 'G12C', 'G12D', 'G12R', 'G12V', 'G13D'), 
    freq = c(sum(tmp$G12A_F1, na.rm = TRUE), 
             sum(tmp$G12C_F1, na.rm = TRUE), 
             sum(tmp$G12D_F1, na.rm = TRUE), 
             sum(tmp$G12R_F1, na.rm = TRUE), 
             sum(tmp$G12V_F1, na.rm = TRUE), 
             sum(tmp$G13D_F1, na.rm = TRUE) )
  )
  
  df = rbind(cd4.df, cd8.df)
  return(df)
})
exp.res = bind_rows(exp.res, .id = 'Subject')
exp.res2 = lapply(spl, FUN = function(s){
  CD4.types = c(unique(final_table$consensus_annot)[grep('CD4', unique(final_table$consensus_annot))], 'Treg')
  CD8.types = unique(final_table$consensus_annot)[grep('CD8', unique(final_table$consensus_annot))]
  
  tmp = s[s$consensus_annot %in% CD4.types,]
  cd4.df = data.frame(
    type = 'CD4', 
    mKRAS = rep(c('G12A', 'G12C', 'G12D', 'G12R', 'G12V', 'G13D'), each = nrow(tmp)), 
    freq = c(tmp$G12A_F1, 
             tmp$G12C_F1,
             tmp$G12D_F1,  
             tmp$G12R_F1, 
            tmp$G12V_F1, 
             tmp$G13D_F1 ), 
    tcrb.id.adaptive =rep(tmp$adaptive.tcrb.id, 6),
    tcrb.id.10x = rep(paste0(tmp$TRBV, '_', tmp$CDR3.beta.aa), 6)
  )
  
  tmp = s[s$consensus_annot %in% CD8.types,]
  cd8.df = data.frame(
    type = 'CD8', 
    mKRAS = rep(c('G12A', 'G12C', 'G12D', 'G12R', 'G12V', 'G13D'), each = nrow(tmp)), 
    freq = c(tmp$G12A_F1, 
             tmp$G12C_F1,
             tmp$G12D_F1,  
             tmp$G12R_F1, 
             tmp$G12V_F1, 
             tmp$G13D_F1 ), 
    tcrb.id.adaptive =rep(tmp$adaptive.tcrb.id, 6),
    tcrb.id.10x = rep(paste0(tmp$TRBV, '_', tmp$CDR3.beta.aa), 6)
  )
  
  df = rbind(cd4.df, cd8.df)
  return(df)
})
exp.res2 = bind_rows(exp.res2, .id = 'Subject')
exp.res$Subject = factor(exp.res$Subject, levels = c('J1994_12', 'J1994_5', 'J1994_2', 'J1994_6'))
exp.res2$Subject = factor(exp.res2$Subject, levels = c('J1994_12', 'J1994_5', 'J1994_2', 'J1994_6'))
ggplot(exp.res, aes(x = Subject, y = freq * 100, fill = type)) + 
  geom_col(color = 'black', position = 'dodge')+
  scale_fill_manual(values = c('white', 'black') ) + 
  facet_wrap(~mKRAS, scales = 'free') + 
  theme_prism() + 
  ylab('Percent Expansion') + 
  theme(panel.grid.major = element_line(color = 'grey90')) + 
   coord_flip()

ggplot(exp.res2, aes(x = Subject, y = freq * 100, fill = type)) + 
  geom_boxplot(color = 'black', outlier.shape = NA)+
  scale_fill_manual(values = c('white', 'black') ) + 
  facet_wrap(~mKRAS, scales = 'free') + 
  theme_prism() + 
  ylab('Percent Expansion') + 
  theme(panel.grid.major = element_line(color = 'grey90')) + 
  coord_flip()
ggplot(exp.res2, aes(x = Subject, y = freq * 100, fill = type)) + 
  geom_violin(color = 'black', outlier.shape = NA)+
  scale_fill_manual(values = c('white', 'black') ) + 
  facet_wrap(~mKRAS, scales = 'free') + 
  theme_prism() + 
  ylab('Percent Expansion') + 
  theme(panel.grid.major = element_line(color = 'grey90')) + 
  coord_flip()

exp.res2 = exp.res2[!is.na(exp.res2$freq),]
write.csv(exp.res2, file = paste0('./code_outputs/',Sys.Date(), 'CD4.CD8.Exp.Freqs.csv'),row.names = FALSE )