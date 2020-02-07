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
#' @export
loadtable = function(file,vars,suffix){
  df<-data.table::fread(file,header=T,sep="\t", stringsAsFactors=FALSE, na.strings="")
  names(df)[names(df) %in% vars]<-paste(vars,suffix,sep="")
  df<-df[,names(df) %in% c("record_id",paste(vars,suffix,sep="")) , with = FALSE]
  return(df)
}
#CreateUKBiobankPhentoypes(3,2 ,UKbioDataset,dfmaster_SQL_merge=dfmaster_SQL_merge,dfSQL_TraitTable,Outputdir )


#' @export
loadGPTable <- function(UKbioDataset,
                        gp_file,
                        cols_tokeep=c("eid","event_dt","read_2","read_3"),
                        cols_rename=c("n_eid","event_dt","read_2","read_3"),
                        mindate = "1930-01-01",
                        maxdate = "2021-01-01"
){
  dfgp <- fread(gp_file,select=cols_tokeep)
  names(dfgp) <- cols_rename
  dfgp$event_dt = as.Date(dfgp$event_dt, "%d/%m/%Y")
  dfgp$epiend_1 = dfgp$event_dt
  UKbioDataset_subset <- UKbioDataset[c("n_eid", names(UKbioDataset[grepl("ts_53_", names(UKbioDataset))]))]
  dfgp <- merge(dfgp, UKbioDataset_subset, by.x = "n_eid",  by.y = "n_eid")
  dfgp<-as.data.frame(dfgp)

  print(paste("missing dates for ", sum(is.na(dfgp$event_dt)),"of",nrow(dfgp),"entries = ",100*sum(is.na(dfgp$event_dt))/nrow(dfgp),"%, excluding these." ))
  dfgp <- dfgp[!is.na(dfgp$event_dt),]

  print(paste("QC dates.. "))
  dfgp <- subset(dfgp, event_dt > as.Date(mindate) ) # removing before 1930-ish.. some 1900, 1902, 1903 observations..
  dfgp <- subset(dfgp, event_dt < as.Date(maxdate) ) # removing after today-ish.. some 2037 observations.

  print(paste("#individuals",length(unique(dfgp$n_eid))))
  # hist( unique(dfgp$event_dt), "years", freq = TRUE)
  # hist( dfgp$event_dt, "years", freq = TRUE,bars=100)
  return(dfgp)
}

#' LoadHesinTable_v2, v2 of uk biobank in hospital tables
#'
#' Combining Hesin tables for analyses by CreateUKBiobankPhentoypes()
#'
#' @name LoadHesinTable_v2
#' @param dfUKbioDataset Dataframe:The main UKBiobank dataset, currently only STATA files are allowed that are processed using the ukb_conv tool. Load with as.data.frame(read.dta13(file,convert.dates = TRUE))
#' @param fhesin containing tsv file hes records
#' @param fhesin_diag tsv file containing data on secondary ICD10/ICD9 codes
#' @param fhesin_oper tsv file containing data on secondary OPCS codes
#' @keywords LoadHesinTable_v2
#' @export
#' @examples
#'
#'
LoadHesinTable_v2 <- function(UKbioDataset,fhesin, fhesin_diag,fhesin_oper){

  # read.
  dfhesin <- data.frame(fread(fhesin,header=T,sep="\t", stringsAsFactors=FALSE, na.strings=""))
  dfhesin_diag <- data.frame(fread(fhesin_diag,header=T,sep="\t", stringsAsFactors=FALSE, na.strings=""))
  dfhesin_oper <- data.frame(fread(fhesin_oper,header=T,sep="\t", stringsAsFactors=FALSE, na.strings=""))

  # c("eid", "ins_index", "dsource", "source", "epistart", "epiend",
  #   "epidur", "bedyear", "epistat", "epitype", "epiorder", "spell_index",
  #   "spell_seq", "spelbgin", "spelend", "speldur", "pctcode", "gpprpct",
  #   "category", "elecdate", "elecdur", "admidate", "admimeth_uni",
  #   "admimeth", "admisorc_uni", "admisorc", "firstreg", "classpat_uni",
  #   "classpat", "intmanag_uni", "intmanag", "mainspef_uni", "mainspef",
  #   "tretspef_uni", "tretspef", "operstat", "disdate", "dismeth_uni",
  #   "dismeth", "disdest_uni", "disdest", "carersi")
  #dfhesin[dfhesin$source==10,]

  dfhesin <- dfhesin[,c("eid", "ins_index", "dsource", "source", "epistart", "epiend", "admidate",  "disdate")]
  dfhesin_diag<- dfhesin_diag[,c("eid", "ins_index", "arr_index", "level", "diag_icd9","diag_icd10")] # c("eid", "ins_index", "arr_index", "level", "diag_icd9", "diag_icd9_nb","diag_icd10", "diag_icd10_nb")
  dfhesin_oper<- dfhesin_oper[,c("eid", "ins_index", "arr_index", "level", "opdate", "oper3", "oper4")] # c("eid", "ins_index", "arr_index", "level", "opdate", "oper3", "oper3_nb", "oper4", "oper4_nb", "posopdur", "preopdur")
  print("merging hesin + diagnosis, please wait..")
  dfhes <- merge(dfhesin,dfhesin_diag,by = c("eid","ins_index"))
  print("merging hesin + operation, please wait..")
  dfhes <- merge(dfhes,dfhesin_oper,by = c("eid","ins_index","arr_index"),all=T)

  # convert dates
  dfhes$epistart <- as.Date(as.character(dfhes$epistart),format="%Y%m%d")
  dfhes$epiend <- as.Date(as.character(dfhes$epiend),format="%Y%m%d")
  dfhes$admidate <- as.Date(as.character(dfhes$admidate),format="%Y%m%d")
  dfhes$disdate <- as.Date(as.character(dfhes$disdate),format="%Y%m%d")
  dfhes$opdate <- as.Date(as.character(dfhes$opdate),format="%d/%m/%Y")

  dfhes[is.na(dfhes$epistart),"epistart"] <- dfhes[is.na(dfhes$epistart),"admidate"]
  dfhes[is.na(dfhes$epistart),"epiend"] <- dfhes[is.na(dfhes$epistart),"disdate"]
  dfhes[is.na(dfhes$epistart),"epistart"] <- dfhes[is.na(dfhes$epistart),"opdate"]
  dfhes[is.na(dfhes$opdate),"opdate"] <- dfhes[is.na(dfhes$opdate),"epistart"]
  dfhes<- dfhes[!is.na(dfhes$epistart),]

  dfhes$level <- rowMins(as.matrix(dfhes[,c("level.x","level.y")]),na.rm = T)
  print("merging ukb dates + hesin, please wait..")
  UKbioDataset_subset<-UKbioDataset[c("n_eid",names(UKbioDataset[grepl("ts_53_",names(UKbioDataset))]))]
  dfhes<-merge(dfhes,UKbioDataset_subset,by.x = "eid",by.y = "n_eid")
  dfhes <- data.frame(dfhes)

  names(dfhes)[which(names(dfhes) %in% c("eid","epistart","epiend","diag_icd9","diag_icd10","opdate","oper3" ,"oper4"    ))] <- c("n_eid",'epistart_1',"epiend_1","diag_icd9_1","diag_icd10_1","opdate_1","oper3_1" ,"oper4_1" )

  dfhes$diag_icd9_2 <- dfhes$diag_icd9_1
  dfhes[dfhes$level==2,"diag_icd9_1"] <- NA # consistent with Version 1 of the data, not very elegant...but too lazy to change everything, and since I  should recode everything there is no point..
  dfhes[dfhes$level==1,"diag_icd9_2"] <- NA

  dfhes$diag_icd10_2 <- dfhes$diag_icd10_1
  dfhes[dfhes$level==2,"diag_icd10_1"] <- NA
  dfhes[dfhes$level==1,"diag_icd10_2"] <- NA

  dfhes$oper3_2 <- dfhes$oper3_1
  dfhes[dfhes$level==2,"oper3_1"] <- NA
  dfhes[dfhes$level==1,"oper3_2"] <- NA

  dfhes$oper4_2 <- dfhes$oper4_1
  dfhes[dfhes$level==2,"oper4_1"] <- NA
  dfhes[dfhes$level==1,"oper4_2"] <- NA

  return(dfhes)

}

