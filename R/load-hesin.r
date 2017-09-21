##########################################################################################
#### -----#### -----#### -----#### -----#### -----#### -----#### -----#### -----#### -----
#' LoadHesinTable 
#'
#' Combining Hesin tables for analyses by CreateUKBiobankPhentoypes() 
#' 
#' @name LoadHesinTable
#' @param dfUKbioDataset Dataframe:The main UKBiobank dataset, currently only STATA files are allowed that are processed using the ukb_conv tool. Load with as.data.frame(read.dta13(file,convert.dates = TRUE))
#' @param hesin_file containing tsv file primary diagnoses for in hospital data, diagnoses + operation
#' @param hesin_diagicd10_file tsv file containing data on secondary ICD10 codes  for in hospital data
#' @param hesin_diagicd9_file tsv file containing data on secondary ICD9 codes for in hospital data
#' @param hesin_oper_file tsv file containing data on secondary operation codes for in hospital data
#' @keywords LoadHesinTable
#' @export
#' @examples
#' 
#'
LoadHesinTable = function(UKbioDataset,hesin_file,hesin_diagicd10_file,hesin_diagicd9_file,hesin_oper_file) {
  ### what are the _nb//addendum variables?
  #### HESIN 
  dfhesin<-data.table::fread(hesin_file,header=T,sep="\t", stringsAsFactors=FALSE, na.strings="")

  dfhesin[is.na(dfhesin$epistart),]$epistart<-dfhesin[is.na(dfhesin$epistart),]$admidate
  dfhesin[is.na(dfhesin$opdate),]$opdate<-dfhesin[is.na(dfhesin$opdate),]$epistart
  VctHesinvars<-c("admidate","cause_icd10","diag_icd10","diag_icd9","disdate","epiend","epistart","opdate","oper4","operstat","posopdur","preopdur")
  names(dfhesin)[ names(dfhesin) %in% VctHesinvars]<-paste(VctHesinvars,"_1",sep="")
  dfhesin<-dfhesin[,names(dfhesin) %in% c("eid","record_id",paste(VctHesinvars,"_1",sep="")), with = FALSE]
  #### convert to dates:
  DateVars=c("epistart_1","opdate_1")
  dfhesin[,(DateVars):= lapply(.SD, as.Date), .SDcols = DateVars]
  ### merge with visit dates:

  UKbioDataset_subset<-UKbioDataset[c("n_eid",names(UKbioDataset[grepl("ts_53_",names(UKbioDataset))]))]
  dfhesin<-merge(dfhesin,UKbioDataset_subset,by.x = "eid",by.y = "n_eid")
  #### 
  dfhesin_diagicd10_file<-loadtable(hesin_diagicd10_file,c("diag_icd10"),"_2")
  dfhesin_oper_file<-loadtable(hesin_oper_file,c("opdate","oper4"),"_2")
  dfhesin_diagicd9_file<-loadtable(hesin_diagicd9_file,c("diag_icd9"),"_2")
  
  LstTbls<-list(dfhesin,dfhesin_diagicd10_file,dfhesin_oper_file,dfhesin_diagicd9_file)
  dfMerged<-Reduce(function(x,y) {merge(x,y,by="record_id",all=TRUE,allow.cartesian=T)}, LstTbls)
  names(dfMerged)[2]<-"n_eid"
  dfMerged<-as.data.frame(dfMerged[ !is.na(dfMerged$epistart_1) ])### THERE ARE SOME THAT HAVE NO EPISTART DATE, BUT THOSE ARE MOSTLY ACCIDENTS ICD10, NO EXTRA INFO, LEAVING THEM OUT. 
  
  return(dfMerged)
  ##"opdate_1" 
  #dfMerged$opdate_1!=dfMerged$opdate_2
}
##########################################################################################
loadtable = function(file,vars,suffix){
  df<-data.table::fread(file,header=T,sep="\t", stringsAsFactors=FALSE, na.strings="")
  names(df)[names(df) %in% vars]<-paste(vars,suffix,sep="")
  df<-df[,names(df) %in% c("record_id",paste(vars,suffix,sep="")) , with = FALSE]
  return(df)
}
#CreateUKBiobankPhentoypes(3,2 ,UKbioDataset,dfmaster_SQL_merge=dfmaster_SQL_merge,dfSQL_TraitTable,Outputdir )




