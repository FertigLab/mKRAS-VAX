#ProcessAdaptiveVgenes Function modified from ProcessAdaptiveFile.R written as part of GIANA. 
#Other functions are original. 
#Accessed from GIANA Github 17 September 2023 
# https://github.com/s175573/GIANA/blob/master/ProcessAdaptiveFile.R
ProcessAdaptiveVgenes <- function(input){

  #AAG 17 September 2023 
  Vgene_o=c(1,2,9,13:19,26:28,30)
  Vgene_o=as.character(Vgene_o)

  gsub('TCRBV[0]{0,1}','TRBV', input)->tmpV
  ## Multiple calls
  vv.m=grep('/',tmpV)
  if(length(vv.m)>0){
    TMP=strsplit(tmpV[vv.m],'/')
    tmpV[vv.m]=unlist(sapply(TMP,function(x)x[1]))
  }
  vv.m=grep(',',tmpV)
  if(length(vv.m)>0){
    TMP=strsplit(tmpV[vv.m],',')
    tmpV[vv.m]=unlist(sapply(TMP,function(x)x[1]))
  }
  tmpV=gsub('-0','-',tmpV)
  v_digit=grep('\\*',tmpV)
  if(length(v_digit)>0)tmpV[- v_digit] = paste(tmpV[- v_digit],'*01',sep='') else tmpV=paste(tmpV,'*01',sep='')
  ## 1. Orphan V genes do not have "-1", need to remove
  Vnumbers=gsub('TRBV','',tmpV)
  Vnumbers=gsub('\\*.+','',Vnumbers)
  Vnumbers=gsub('-.+','',Vnumbers)
  vv.o=which(Vnumbers %in% Vgene_o)
  tmpV1=tmpV
  tmpV1[vv.o]=gsub('-1','',tmpV1[vv.o])
  ## 2. Non-orphan V genes but without "-1", need to add
  vv.no=which(! Vnumbers %in% Vgene_o)
  vv.non=grep('-',tmpV1[vv.no])
  if(length(vv.non)>0)tmpV1[vv.no][-vv.non]=gsub('\\*01','-1*01',tmpV1[vv.no][-vv.non])
  output=tmpV1
  return(output)
}
CheckTRBVEquivalence <- function(input1, input2){
  if(length(input1) == length(input2)) N = length(input1)
  v1 = digestTRBV(trbv =input1)
  v2 = digestTRBV(trbv = input2)
  
  good = vector(mode = 'logical', length = N)
  good = (v1$family == v2$family & v1$gene == v2$gene) | 
    (v1$family == v2$family & (v1$gene == 0 | v2$gene == 0) ) | 
    (v1$family == v2$family & v1$gene == v2$gene )
  
  good2 = ((input1 == ProcessAdaptiveVgenes(input2)) | 
             (ProcessAdaptiveVgenes(input1) == input2) | 
             (ProcessAdaptiveVgenes(input1) == ProcessAdaptiveVgenes(input2)))
  good3 = good + good2  > 0
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

