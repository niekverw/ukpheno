
#### -----#### -----#### -----#### -----#### -----#### -----#### -----#### -----#### -----
#' CreateUKBiobankPhentoypes - Function
#'
#' CreateUKBiobankPhentoypes creates dichotomous phenotypes for UK Biobank using definitions that are based on the available variables (e.g. HES-data, mortality, toucschreen,  SR data). The current output includes variables on history, study visit, future, time-to-first-event; WORK IN PROGRESS
#' @name CreateUKBiobankPhentoypes
#' @param Nvisits #Number of follow-up visits in your dataset (a.t.m 3, baseline,1 and 2) .
#' @param visitreference #Number visit that should be used as reference; should be 0 (baseline),1 or 2
#' @param dfUKbioDataset Dataframe:The main UKBiobank dataset, currently only STATA files are allowed that are processed using the ukb_conv tool. Load with as.data.frame(read.dta13(file,convert.dates = TRUE,convert.factors=FALSE))
#' @param dfmaster_SQL_merge Dataframe: The Hesin table of UKBiobank that can be obtained by SQL; load via (WORK IN PROGRESS)
#' @param dfTraitTable Dataframe: A Trait-definition table. (EXAMPLE: WORK IN PROGRESS)
#' @param Outputdir Output directory
#' @param VctOutputIndividualColumns set VctOutputIndividualColumns=c("TS","SR","TS_RX","RX","LAB") to output all individual columns, default is 'VctOutputIndividualColumns=c()'
#' @keywords CreateUKBiobankPhentoypes
#' @references <Please contact mail@niekverweij.com, this package is not yet ready to be distributed.>
#' @importFrom parallel mclapply detectCores
#' @importFrom matrixStats rowMins
#' @importFrom matrixStats rowMaxs
#' @importFrom dplyr group_by mutate filter
#' @importFrom readstata13 save.dta13
#' @import data.table
#' @export
#' @examples
#' 
#' library(gdata) 
#' library(readstata13)
#' library(CreateUKBiobankPhentoypes)
#' Nvisits=3
#' visitreference=0
#' 
#' Outputdir<-getwd()
#' dfUKbioDataset<-as.data.frame(read.dta13("/path/to/UKbioDataset_file.dta",convert.dates = TRUE,convert.factors=FALSE))
#' 
#' hesin_file="/path/to/hesin-4july2016.tsv" 
#' hesin_diagicd10_file="/path/to/hesin_diagicd10-4july2016.tsv"
#' hesin_diagicd9_file="/path/to/hesin_diagicd9_22_03_2016.tsv"
#' hesin_oper_file="/path/to/hesin_oper-4july2016.tsv"
#' dfhesintables<-LoadHesinTable(dfUKbioDataset,hesin_file,hesin_diagicd10_file,hesin_diagicd9_file,hesin_oper_file)
#' 
#' ## load phenotype definitions. 
#' data(dfDefinitions)
#'
#' 
#' 
#'
#' #run with minimum variables:
#'CreateUKBiobankPhentoypes(Nvisits,visitreference,UKbioDataset,dfhesintables,dfDefinitions,Outputdir)
#' #which is the same as: 
#' CreateUKBiobankPhentoypes(Nvisits,visitreference,UKbioDataset,dfhesintables,dfDefinitions,Outputdir,VctOutputIndividualColumns=c(),VctAllUKBVDefinitionColumns=c("TS","SR","TS_RX","RX","LAB") )
#' 
#' #to output the variables based on RX (medication) in an *_RX file and exclude RX from the aggregate *_UKBV files run the following:
#' CreateUKBiobankPhentoypes(Nvisits,visitreference,UKbioDataset,dfhesintables,dfDefinitions,Outputdir,VctOutputIndividualColumns=c("RX",VctAllUKBVDefinitionColumns=c("TS","SR","TS_RX","LAB") )
################################################################################################################################################
##library(readstata13)
##library(dplyr)
##library(data.table) #
##library(parallel) ## native to R.
##library(gdata)   ## for excel
##library(matrixStats) ## for rowMins/max function
###########################################################
CreateUKBiobankPhentoypes<-function(Nvisits,visitreference,UKbioDataset,dfmaster_SQL_merge,dfDefinitions,Outputdir,VctOutputIndividualColumns=c() ){

  #### DEFINE COLUMNS
  VctAllUKBVDefinitionColumns=c("TS","SR","TS_RX","RX","LAB") #set this variable to a selection of columns (dfDefinition columns) to be outputted by the _UKBV variable, default is 'VctAllUKBVDefinitionColumns=c("TS","SR","TS_RX","RX","LAB")'
  VctAllHESINDefinitionColumns=c("ICD10CODES","ICD9CODES","OPERCODES")
  StrTSAgeColumn="TS_AGE_DIAG_COLNAME"
  StrSRcolumns<-c("READCODES","n_20001_","n_20002_","n_20003_","n_20004_")
  VctAllColumns=c("DEPENDENCY",VctAllUKBVDefinitionColumns,VctAllHESINDefinitionColumns,StrTSAgeColumn,StrSRcolumns)
  #dfmaster_SQL_merge<-dfhesintables
  ### OPTIONAL SETTINGS
  OutputIndividualItems=FALSE
  ###parallel=2 <-need to fix.
  ################################################################################################################################################
  ###### RUN following loop:
  ################################################################################################################################################
  
  
  
  visitdt=paste("ts_53_",visitreference,"_0",sep="")
  
  print("processing dfDefinitions")
  dfDefinitions<-ProcessDfDefinitions(dfDefinitions,VctAllColumns)
  ### TODO: check if n_52_0_0 and n_34_0_0 are available == month/year of birth 
  ### TODO: add function to check if columns in dfdefinitions are availabe. 
  ### TODO: check how many hits per ICD10code --> logfile. 
  ### TODO: evaluate what happens if uncommented and refvisit=2 where some ts_53 are missing.
  #dfmaster_TSDEATHMEDICD10_visitdtonly<-UKbioDataset[!is.na(UKbioDataset[visitdt]),]
  dfmaster_TSDEATHMEDICD10_visitdtonly<-UKbioDataset
  
  print("checking definitions; for(i in 1:nrow(dfDefinitions))")
  for(i in 1:nrow(dfDefinitions)) {
    row <- dfDefinitions[i,]
    row <-ConvertFactorsToStringReplaceNAInDf(row) ## convert factor->string
    row$DESCRIPTION<-gsub(pattern = '[^a-zA-Z0-9_ ()=|;,+-\\&]+',"",row$DESCRIPTION) ## remove non-alphanumeric except ()=|;,+-\\&
    ## check if columns are empty. 
        sd
    StrTrait<-row$TRAIT
    if (grepl("_",StrTrait)) {
      stop(paste("Trait: ",StrTrait," contains '_' !!,",sep=""))
      
    }
    
    ########################################################
    ########################################################
    ### HESIN DATA : exact dates.
    Strcatagory="HESIN"
    ########################################################
    
    ########################################################
    ####### OPERATION:
    #######################################################
    if( !is.na(row$OPERCODES)  ){
      #### SETTINGS:
      print(paste("   ..finding operationcodes",StrTrait))
      VctCodes<-unlist(strsplit(row$OPERCODES,","))
      StrColumnForHescodes<-c("oper4_1","oper4_2" )
      StrColumnnameForEventDate<-"opdate_1"
      
      ####################################
      ### HES: FOLLOW UP + HISTORY VARIABLES:
      StrName="Ohes"
      StrColumnnameForVisitDate<-visitdt
      StrDescription<-paste(row$DESCRIPTION,"- OPER(HES)",paste(VctCodes,collapse=","))
      
      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      Outcome_HES(dfmaster_SQL_merge,paste(StrTrait,StrName,sep="_"),StrDescription,VctCodes,StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile)
      
    }
    ########################################################
    ####### ICD10:
    ########################################################
    if( !is.na(row$ICD10CODES) ) {
      #### SETTINGS:
      print("    ..finding icd10diagnosiscodes")
      VctCodes<-unlist(strsplit(row$ICD10CODE,","))
      ####################################
      ### HESIN ICD10: FOLLOW UP + HISTORY VARIABLES:
      StrName="Dhes"
      StrColumnForHescodes<-c("diag_icd10_1","diag_icd10_2" )
      StrColumnnameForEventDate<-"epistart_1"
      StrDescription<-paste(row$DESCRIPTION,"- ICD10(HES) -",paste(VctCodes,collapse=","))
      
      StrColumnnameForVisitDate<-visitdt
      
      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      Outcome_HES(dfmaster_SQL_merge,paste(StrTrait,StrName,sep="_"),StrDescription,VctCodes,StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile)
      
      ####################################
      #### DEATH ICD10:
      StrName<-"DO"
      StrFieldcode<-"40001|40002"
      StrFieldDateColumn<-"ts_40000_0_0"
      StrDescription<-paste(row$DESCRIPTION," - ","Prim+sec Death Of Disorder ICD10(DB):",paste(VctCodes,collapse=","),sep="")
      
      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      Death_Masterset(dfmaster_TSDEATHMEDICD10_visitdtonly,StrTrait,StrDescription,StrName,VctCodes,StrFieldcode,StrFieldDateColumn,visitdt,StataOutputFile)
      
      ####################################
      #### DEATH ICD10:
      StrName<-"DOp"
      StrFieldcode<-"40001"
      StrFieldDateColumn<-"ts_40000_0_0"
      StrDescription<-paste(row$DESCRIPTION," - ","Prim Death Of Disorder ICD10(DB):",paste(VctCodes,collapse=","),sep="")
      
      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      Death_Masterset(dfmaster_TSDEATHMEDICD10_visitdtonly,StrTrait,StrDescription,StrName,VctCodes,StrFieldcode,StrFieldDateColumn,visitdt,StataOutputFile)
    }
    
    ########################################################
    ####### ICD9:
    ########################################################
    
    if( !is.na(row$ICD9CODES) ) {
      #### SETTINGS:
      print("    ..finding icd9diagnosiscodes")
      VctCodes<-unlist(strsplit(row$ICD9CODE,","))
      
      ####################################
      ### HESIN ICD9: FOLLOW UP + HISTORY VARIABLES:
      StrName="D9hes"
      StrColumnForHescodes<-c("diag_icd9_1","diag_icd9_2" )
      StrColumnnameForEventDate<-"epistart_1"
      StrDescription<-paste(row$DESCRIPTION,"- ICD9(HES) -",paste(VctCodes,collapse=","))
      StrColumnnameForVisitDate<-visitdt
      
      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      Outcome_HES(dfmaster_SQL_merge,paste(StrTrait,StrName,sep="_"),StrDescription,VctCodes,StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile)
      
    }
    
    ########################################################
    ####### TS(touchscreen//SR(selfreported) // RX // LAB: only based on visit dates. ..
    Strcatagory="UKBVisit"
    
    if(  !is.na(row$TS_AGE_DIAG_COLNAME)) {
      #### SETTINGS:
      print("   ..finding TS; creating age of diagnosis composit for TS")
      
      StrName<-"TSAd"
      VctTSAcolumns<-unlist(strsplit(row$TS_AGE_DIAG_COLNAME,","))
      StrDescription<-paste(row$DESCRIPTION,"- TS mean age of diagnosis in days -",paste(VctTSAcolumns,collapse =",")  )
      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      
      Query_TSAmin_masterset(dfmaster_TSDEATHMEDICD10_visitdtonly,StrTrait,StrDescription,StrName,VctTSAcolumns,visitdt,StataOutputFile)
    }
    
    
    ############### OUTPUT INDIVIDUAL ITEMS IF NEEDED: #############################
    for(StrDefinitionColumn in VctOutputIndividualColumns){
      print(StrDefinitionColumn)
      if(  !is.na(row[,StrDefinitionColumn])) {
        #### SETTINGS:
        print(paste("   ..finding",StrDefinitionColumn))
              VctTSconditions<-unlist(strsplit(row[,StrDefinitionColumn],","))
              StrDescription<-paste(row$DESCRIPTION,"-",StrDefinitionColumn,"-",paste(VctTSconditions,collapse =",")  )
              dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
              StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrDefinitionColumn,sep="_"),".dta",sep="")
              Query_Masterset2(dfmaster_TSDEATHMEDICD10_visitdtonly,StrTrait,StrDescription,StrDefinitionColumn,VctTSconditions,StataOutputFile,visitreference,Nvisits)
              
      }
    }
 
    ############### OUTPUT ALL IN UKBV VARIABLE:: #############################
    UKBV<-paste(as.vector(row[,VctAllUKBVDefinitionColumns]),collapse=",")
    
    if( nchar(UKBV)>= length(VctAllUKBVDefinitionColumns)*3 ) { ### CHECK IF UKBV columns  have something in it (*3 ="NA," for every column)
      print("   ..finding all ukbio visit terms.")
      StrName<-"UKBV"
      VctTSconditions<-unlist(strsplit(UKBV,","))
      VctTSconditions<-VctTSconditions[VctTSconditions!="NA"]
      
      StrDescription<-paste(row$DESCRIPTION,"-",StrName,"-",paste(VctTSconditions,collapse =",")  ) ### TOO LONG Description..
      StrDescription<-paste(row$DESCRIPTION,"-",StrName)
      
      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,"/ALL",sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/ALL","/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      Query_Masterset2(dfmaster_TSDEATHMEDICD10_visitdtonly,StrTrait,StrDescription,StrName,VctTSconditions,StataOutputFile,visitreference,Nvisits)
      
    }
    
    
    ##########################################
    ### MERGING Files: 
    ##########################################
    
    ###### MERGING ICD10/9 and death follow up.
    print(visitdt)
    filesEpisode<-list.files(path=paste(Outputdir,"/",visitdt,"/HESIN",sep=""),pattern=paste("^",StrTrait,"(_)(.*)dta$",sep=""), full.names=TRUE)
   
    
    if( length(filesEpisode)>0) {
      HEScodes=paste(row$ICD10CODES,row$ICD9CODES,row$OPERCODES,sep="|")
      OutputdirMerged=paste(Outputdir,"/",visitdt,"/HESIN/merged/",sep="")
      StataOutputFile=paste(OutputdirMerged,"/",StrTrait,"_merged.dta",sep="")
      
      if( file.exists(StataOutputFile )) {### CHECK IF FILE EXISTS: 
        print(paste(StataOutputFile," ... already exists!, skipping"))
      } else {
      dfMerged<-MultiMergeEpisodeData(filesEpisode,StrTrait,row$DESCRIPTION,HEScodes)
      dir.create(OutputdirMerged, showWarnings = FALSE, recursive = TRUE)
      if(is.data.frame(dfMerged)){      save.dta13(dfMerged,StataOutputFile,compress = TRUE) }
      }
    }
    
    #### MERGING ICD merged , UKBV and TSAmin
    fileICDmerged<-list.files(path=paste(Outputdir,"/",visitdt,"/HESIN/merged",sep=""),pattern=paste("^",StrTrait,"(_)(.*)dta$",sep=""), full.names=TRUE)
    fileUKBV<-list.files(path=paste(Outputdir,"/",visitdt,"/UKBVisit/ALL",sep=""),pattern=paste("^",StrTrait,"(_)(UKBV.*)dta$",sep=""), full.names=TRUE)
    fileTSA<-list.files(path=paste(Outputdir,"/",visitdt,"/UKBVisit/",sep=""),pattern=paste("^",StrTrait,"(_)(TSA.*)dta$",sep=""), full.names=TRUE)
    files=c(fileICDmerged,fileUKBV,fileTSA)
    print(files)
    if (length(files)>0){
      dir.create( paste(Outputdir,"/",visitdt,"/ALL","/",sep=""), showWarnings = FALSE, recursive = TRUE)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/ALL","/",paste(StrTrait,"ALL",sep="_"),".dta",sep="")
      
      if( file.exists(StataOutputFile )) { ### CHECK IF FILE EXISTS: 
        print(paste(StataOutputFile," ... already exists!, skipping"))
      } else {
      MultiMergeHESIN_UKBV(files,StrTrait,row$DESCRIPTION,HEScodes,StataOutputFile)
      }
      #  MultiMergeUKBV(files,StrTrait,row$DESCRIPTION,HEScodes,StataOutputFile)
    }
    ##########################################
    
  }
  
}





Outcome_HES<-function(dfmaster_SQL_merge,StrTrait,StrDescription,VctCodes,StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile){
  ### CHECK IF FILE EXISTS: 
  if( file.exists(StataOutputFile )) {
    print(paste(StataOutputFile," ... already exists!, skipping"))
    return(0)     
  }
  print(paste("Outcome_HES , checking HES:",StrTrait,":",StrDescription,"-",StrColumnnameForVisitDate) )
  
  DiagnosisVars=c("n_eid","HXn","HXd","FUn","FUd") ### variables you wil create and you want to keep later... 
  #### RENAME COLUMNS  
  names(dfmaster_SQL_merge)[names(dfmaster_SQL_merge)==StrColumnnameForEventDate]<-"event_date"
  names(dfmaster_SQL_merge)[names(dfmaster_SQL_merge)==StrColumnnameForVisitDate]<-"visit_date"
  ## remove records that have no visit date (relevant for  ts_53_1_0 ,  ts_53_2_0). 
  dfmaster_SQL_merge_tmp<-dfmaster_SQL_merge[!is.na(dfmaster_SQL_merge$visit_date),]

  ########### GREP CODE AT THE COLUMNS 'StrColumnForHescodes' ;; TODO grep should be ^: paste(sep="","^",VctCodes, collapse='|')
  #grepoper=c()
  #for (i in StrColumnForHescodes ){
  #  grepoper<-c(grepoper,grep(paste(VctCodes, collapse='|'), dfmaster_SQL_merge_tmp[,i], ignore.case=TRUE))
  #}
  grepoper <-unlist( mclapply(  dfmaster_SQL_merge_tmp[, StrColumnForHescodes], function(col) grep(paste(sep="","^",VctCodes, collapse='|'), col, ignore.case=TRUE),mc.cores =detectCores()/2 ) ) ## PARALLEL of the above. 

  dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_tmp[ unique(c(grepoper)),]
  dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_oper %>% group_by(n_eid) %>% filter (! duplicated(event_date)) ### remove duplicate dates.

  ## History vars: 
  dfmaster_SQL_merge_oper$HX  <-ifelse(dfmaster_SQL_merge_oper$visit_date >=  dfmaster_SQL_merge_oper$event_date, 1, 0)      ### SET HISTORY VARIABLE (0/1)
  dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_oper %>% group_by(n_eid) %>% mutate(HXn =sum(HX))       ### Count number of history events

  dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_oper %>%  group_by(n_eid) %>%  mutate(HX = max(HX) )  ### replace HX for every duplicating n_eid
  dfmaster_SQL_merge_oper$HXd <-ifelse(dfmaster_SQL_merge_oper$visit_date >= dfmaster_SQL_merge_oper$event_date, difftime(dfmaster_SQL_merge_oper$event_date, dfmaster_SQL_merge_oper$visit_date ,units="days") , NA) ### SET FUP time (days) - first event in history 
  dfmaster_SQL_merge_oper     <-suppressWarnings(  dfmaster_SQL_merge_oper %>%  group_by(n_eid) %>%  mutate(HXd = -max(abs(HXd), na.rm=T) ) ) ### replace HXd 
  dfmaster_SQL_merge_oper$HXd[is.infinite(dfmaster_SQL_merge_oper$HXd)]<-NA
  #  dfmaster_SQL_merge_oper$HXd1 <-ifelse(dfmaster_SQL_merge_oper$visit_date >= dfmaster_SQL_merge_oper$event_date, difftime(dfmaster_SQL_merge_oper$event_date, dfmaster_SQL_merge_oper$visit_date ,units="days") , NA) ### SET FUP time (days) - last event in history
  #  dfmaster_SQL_merge_oper     <-suppressWarnings( dfmaster_SQL_merge_oper %>%  group_by(n_eid) %>%  mutate(HXd1 = -min(abs(HXd1), na.rm=T) ) ) ### replace HXd 
  #  dfmaster_SQL_merge_oper$HXd1[is.infinite(dfmaster_SQL_merge_oper$HXd1)]<-NA
  
  
  ## follow up vars; 3170028 1834566 
  dfmaster_SQL_merge_oper$FU  <-ifelse(dfmaster_SQL_merge_oper$visit_date < dfmaster_SQL_merge_oper$event_date, 1, 0)   ### SET FOLLOW UP VARIABLE (0/1)
  dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_oper %>% group_by(n_eid) %>% mutate(FUn =sum(FU)  )  ### Count number of followup events
  dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_oper %>%  group_by(n_eid) %>%  mutate(FU = max(FU) )  ### replace FU for every duplicating n_eid
  dfmaster_SQL_merge_oper$FUd <-ifelse(dfmaster_SQL_merge_oper$visit_date < dfmaster_SQL_merge_oper$event_date, difftime(dfmaster_SQL_merge_oper$event_date, dfmaster_SQL_merge_oper$visit_date ,units="days") , NA) ### SET FUP time (days) - first event in future. 
  dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_oper %>%  group_by(n_eid) %>%  mutate(FUd = min(FUd, na.rm=T) )  ### replace HXd 
  #  dfmaster_SQL_merge_oper$FUd1 <-ifelse(dfmaster_SQL_merge_oper$visit_date < dfmaster_SQL_merge_oper$event_date, difftime(dfmaster_SQL_merge_oper$event_date, dfmaster_SQL_merge_oper$visit_date ,units="days") , NA) ### SET FUP time (days)  - last event in future. 
  #  dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_oper %>%  group_by(n_eid) %>%  mutate(FUd1 = max(FUd1, na.rm=T) )  ### replace HXd 
  
  #dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_oper %>% group_by(n_eid) %>% slice(which.min(FUd)) #### FILTER on unique n_eid by minimum FUD
  ## keep uniq
  dfmaster_SQL_merge_oper<-dfmaster_SQL_merge_oper[!duplicated(dfmaster_SQL_merge_oper[,DiagnosisVars]),][,DiagnosisVars] ## FILTER UNIQ + KEEP ONLY COLUMNS NEEDED.
  
  ### ADD DESCRIPTION // RENAME VARIABLES 
  attr(dfmaster_SQL_merge_oper, "var.labels") <- c("Identifier",paste(StrDescription, colnames(dfmaster_SQL_merge_oper[,-1]), sep = " - ") )
  attr(dfmaster_SQL_merge_oper, "var.labels")<-substr(attr(dfmaster_SQL_merge_oper, "var.labels"),0,80)
  colnames(dfmaster_SQL_merge_oper)[-1] <- paste(StrTrait, colnames(dfmaster_SQL_merge_oper[,-1]), sep = "_")
  save.dta13(dfmaster_SQL_merge_oper, StataOutputFile,compress = TRUE)
  return(dfmaster_SQL_merge_oper)
}  


Death_Masterset<-function(dfmaster_TSDEATHMEDICD10_visitdtonly,StrTrait,StrDescription,StrName,VctCodes,StrFieldcode,StrFieldDateColumn,visitdt,StataOutputFile){
  ### CHECK IF FILE EXISTS: 
  if( file.exists(StataOutputFile )) {
    print(paste(StataOutputFile," ... already exists!, skipping"))
    return(0)     
  }
  print(paste("Death_Masterset - ", StrTrait,":",StrDescription,":",StrFieldcode))
  
  StrColumnNames<-names(dfmaster_TSDEATHMEDICD10_visitdtonly[ grepl(StrFieldcode, names(dfmaster_TSDEATHMEDICD10_visitdtonly))])
  #grepoper=c()
  #for (i in StrColumnNames )
  #  grepoper<-c(grepoper,grep(paste(VctCodes, collapse='|'), dfmaster_TSDEATHMEDICD10_visitdtonly[,i], ignore.case=TRUE))
  grepoper <-unlist( mclapply(  dfmaster_TSDEATHMEDICD10_visitdtonly[, StrColumnNames], function(col) grep(paste(VctCodes, collapse='|'), col, ignore.case=TRUE),mc.cores =detectCores()/2 ) ) ### PARALLEL., the abovve is not parallelized 
  dfmaster_TSDEATHMEDICD10_visitdtonly_subset<-dfmaster_TSDEATHMEDICD10_visitdtonly[ unique(c(grepoper)),][c("n_eid",visitdt,StrFieldDateColumn)]
  
  if (is.data.frame(dfmaster_TSDEATHMEDICD10_visitdtonly_subset) && nrow(dfmaster_TSDEATHMEDICD10_visitdtonly_subset)>0){
    dfmaster_TSDEATHMEDICD10_visitdtonly_subset$grep<-1
    ### Calculate days to death_event
    dfmaster_TSDEATHMEDICD10_visitdtonly_subset$DthFud <-as.vector(difftime(dfmaster_TSDEATHMEDICD10_visitdtonly_subset[,StrFieldDateColumn], dfmaster_TSDEATHMEDICD10_visitdtonly_subset[,visitdt] ,units="days"))
    
    #### REMOVING IF DAYS ARE <0 !!!! TODO FIGURE OUT WHY !!!! 
    #n_eid with negative death:
    #1118475
    #1730573 <-- ts_40000_1_0 contains a later date thatn ts_40000_0_0 ?? ( which is the only case.)
    #1507226
    dfmaster_TSDEATHMEDICD10_visitdtonly_subset<-dfmaster_TSDEATHMEDICD10_visitdtonly_subset[dfmaster_TSDEATHMEDICD10_visitdtonly_subset$DthFud>=0,c(1,4,5)] ## keep only n_eid, death=1 and days until death
    
    ### ADD DESCRIPTION // RENAME VARIABLES 
    names(dfmaster_TSDEATHMEDICD10_visitdtonly_subset)<-c("n_eid",paste(StrTrait,StrName,"FUn",sep="_"),paste(StrTrait,StrName,"FUd",sep="_"))
    attr(dfmaster_TSDEATHMEDICD10_visitdtonly_subset, "var.labels") <- c("Identifier",StrDescription,paste(StrDescription,"- FU, Days"))
    attr(dfmaster_TSDEATHMEDICD10_visitdtonly_subset, "var.labels")<-substr(attr(dfmaster_TSDEATHMEDICD10_visitdtonly_subset, "var.labels"),0,80)

    ### filter on NA so when death!=1 (FUd can be missing as visit date for visit 2 is not always available )
    head(dfmaster_TSDEATHMEDICD10_visitdtonly_subset)
    dfmaster_TSDEATHMEDICD10_visitdtonly_subset<-dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ !is.na(dfmaster_TSDEATHMEDICD10_visitdtonly_subset[,2]),]
    head(dfmaster_TSDEATHMEDICD10_visitdtonly_subset)
    
    if(nrow(dfmaster_TSDEATHMEDICD10_visitdtonly_subset)==0){
      print("ERROR: could not find any variables")
      return(0)
    }
    
    save.dta13(dfmaster_TSDEATHMEDICD10_visitdtonly_subset, StataOutputFile,compress = TRUE)
    return(dfmaster_TSDEATHMEDICD10_visitdtonly_subset)
  }  else {
    print("ERROR: could not find any variables");
  }
}

Query_Masterset2<-function(dfmaster_TSDEATHMEDICD10_visitdtonly,StrTrait,StrDescription,StrName,VctTSconditions,StataOutputFile,visitreference,Nvisits){
  ### CHECK IF FILE EXISTS: 
  if( file.exists(StataOutputFile )) {
    print(paste(StataOutputFile," ... already exists!, skipping"))
    return(0)     
  }
  print(paste(StrTrait,":",StrDescription))
  
  VctVisits=0:(Nvisits-1)
  LstColumnnamesForVisitDates=paste("ts_53_",VctVisits,"_0",sep="")
  
  ### Loop over visits. 
  df=dfmaster_TSDEATHMEDICD10_visitdtonly[c("n_eid",LstColumnnamesForVisitDates)] 
  
  for (visit in VctVisits){
    visitcol=paste("_",visit,"_",sep="")
    print(visitcol)
    
    ### FIND MATCHES, store them in grepoper=c():
    grepoper=c()
    #### ======= argument ; can be _n3434_=1|2||3
    for (Condition in VctTSconditions[grep("=",VctTSconditions)]) {
      print( unlist(strsplit(Condition,split="=")) )
      StrTScolumn=unlist(strsplit(Condition,split="="))[1]
      StrTScodes=unlist(strsplit(Condition,split="="))[2]
      StrColumnNames<-names(dfmaster_TSDEATHMEDICD10_visitdtonly[ grepl(StrTScolumn, names(dfmaster_TSDEATHMEDICD10_visitdtonly)) ])
      StrColumnNames_pervisit<-StrColumnNames[grepl(visitcol,StrColumnNames)]
      if ( length(StrColumnNames) == 0)  { stop(paste("ERROR: at least one column missing, there may be more missing, please check!!!: ",paste(StrTScolumn,collapse=",") ) ) }
      
      for (i in StrColumnNames_pervisit ) {
        #for (code in unlist(strsplit(StrTScodes,split="\\|" )) ){
        #  grepoper<- c(grepoper,which(as.numeric(dfmaster_TSDEATHMEDICD10_visitdtonly[,i])==code))
        #}
        grepoper<- c(grepoper,which(as.numeric(dfmaster_TSDEATHMEDICD10_visitdtonly[,i]) %in% unlist(strsplit(StrTScodes,split="\\|" )) ))  ## faster.. but the above is similar to the code below.. should figure out how to merge below with this or something.. 
      }
    }
    #### ≥≥≥≥≥ = argument ; can be _n3434_≥1    ##### FIX ≤
    for (Condition in VctTSconditions[grep("≥",VctTSconditions)]) {
      print( unlist(strsplit(Condition,split="≥")) )
      StrTScolumn=unlist(strsplit(Condition,split="≥"))[1]
      StrTScodes=unlist(strsplit(Condition,split="≥"))[2]
      StrColumnNames<-names(dfmaster_TSDEATHMEDICD10_visitdtonly[ grepl(StrTScolumn, names(dfmaster_TSDEATHMEDICD10_visitdtonly))])
      StrColumnNames_pervisit<-StrColumnNames[grepl(visitcol,StrColumnNames)]
      if ( length(StrColumnNames) == 0)  { stop(paste("ERROR: at least one column missing, there may be more missing, please check!!!: ",paste(StrTScolumn,collapse=",") ) ) }
      
      for (i in StrColumnNames_pervisit ) {
        grepoper<- c(grepoper,which(as.numeric(dfmaster_TSDEATHMEDICD10_visitdtonly[,i])>=StrTScodes))
      }
    }
    
    #### >>>>> = argument ; can be _n3434_ > 1   ##### FIX <
    for (Condition in VctTSconditions[grep(">",VctTSconditions)]) {
      print( unlist(strsplit(Condition,split=">")) )
      StrTScolumn=unlist(strsplit(Condition,split=">"))[1]
      StrTScodes=unlist(strsplit(Condition,split=">"))[2]
      StrColumnNames<-names(dfmaster_TSDEATHMEDICD10_visitdtonly[ grepl(StrTScolumn, names(dfmaster_TSDEATHMEDICD10_visitdtonly))])
      StrColumnNames_pervisit<-StrColumnNames[grepl(visitcol,StrColumnNames)]
      if ( length(StrColumnNames) == 0)  { stop(paste("ERROR: at least one column missing, there may be more missing, please check!!!: ",paste(StrTScolumn,collapse=",") ) ) }
      
      for (i in StrColumnNames_pervisit ) {
        grepoper<- c(grepoper,which(as.numeric(dfmaster_TSDEATHMEDICD10_visitdtonly[,i])>StrTScodes))
      }
    }
    ##########################################
    ### / put it in dataframe:
    if (!length(grepoper)==0){
      df$grep<-NA
      df[unique(c(grepoper)),]$grep<-1
      names(df)[names(df)=="grep"]<-paste(StrName,"_",visit,sep="")
    } else {
      df$grep<-NA
      names(df)[names(df)=="grep"]<-paste(StrName,"_",visit,sep="")
    }
  }
  
  
  df<-df[ !apply(df[(ncol(df)-Nvisits+1):ncol(df)], 1, function(x) all(is.na(x))),]  ### filter on trait==1 ( df[(ncol(df)-3+1):ncol(df)] filters on TS_*_*. should substitute 3 by #visits)
  
  df$HXn<-0
  df$FUn<-0
  df$HXd<-NA
  df$FUd<-NA
  
  # loop over VctVisits 0 1 2
  for (visitnr in VctVisits){   ### replace NA with 0 for the TS_0/1/2 column if they attended the visit. ; leave NA if there is no visit
    Coltrait=paste(StrName,visitnr,sep="_")
    Colvisit=paste("ts_53_",visitnr,"_0",sep="")
    Colreference=paste("ts_53_",visitreference,"_0",sep="")  ## TODO: ADD IN OPTION TO USE CUSTOMIZED VECTOR WITH DATES
    
    df[ is.na(df[Coltrait]) & !is.na(df[Colvisit]) ,Coltrait]<-0 
    
    
    if (visitnr<=visitreference ){ ## TODO: ADD IN OPTION TO USE CUSTOMIZED VECTOR WITH DATES -- replace this loop with something more generic?
      HX_TF<- df[Coltrait]==1 & !is.na(df[Coltrait])
      df[ HX_TF , ]$HXn<-1+df[ HX_TF, ]$HXn
      HXd_TF<-!is.na(df[c(Colvisit)]) & !is.na(df[c(Colreference)]) & df[Colvisit] <= df[Colreference] & df[Coltrait]==1
      df[HXd_TF & is.na(df$HXd), ]$HXd <- df[HXd_TF & is.na(df$HXd), Colvisit]-df[HXd_TF & is.na(df$HXd), Colreference]  ### skip if HXd is filled in, to keep first observation.
      
    }
    
    if (visitnr>visitreference ){
      FU_TF<- df[Coltrait]==1 & !is.na(df[Coltrait])
      df[FU_TF, ]$FUn<- 1+ df[ FU_TF, ]$FUn
      FUd_TF<- !is.na(df[c(Colvisit)]) & !is.na(df[c(Colreference)]) & df[Colvisit] > df[Colreference] & df[Coltrait]==1
      df[FUd_TF & is.na(df$FUd), ]$FUd<- df[FUd_TF & is.na(df$FUd), Colvisit]-df[FUd_TF & is.na(df$FUd), Colreference]  ### skip if FUd is filled in, to keep first observation.
      
    }
    
  }
  
  df$FUdlv<-NA ## calculate days to last visit before FUd and discrepancies between Q's
  df$HXdlv<-NA ## calculate days to last visit before FUd and discrepancies between Q's
  
  df$ERR<-NA ## detect individuals with questionable answers (010, 100,101 etc.)
  df$ERRtmp<-df[,paste(StrName,"0",sep="_")] 
  
  for (visitnr in VctVisits){   ### replace NA with 0 for the TS_0/1/2 column if they attended the visit. ; leave NA if there is no visit
    Coltrait=paste(StrName,visitnr,sep="_")
    Colvisit=paste("ts_53_",visitnr,"_0",sep="")
    Colreference=paste("ts_53_",visitreference,"_0",sep="") 
    
    Colvisitreferencediff<-!is.na(df[c(Colvisit)]) & !is.na(df[c(Colreference)]) & !is.na(df$FUd) & (df[,Colvisit]-df[,Colreference]) < df$FUd
    df[Colvisitreferencediff,]$FUdlv<-df[Colvisitreferencediff,Colvisit ] - df[Colvisitreferencediff,Colreference ]
    
    
    Colvisitreferencediff<-!is.na(df[c(Colvisit)]) & !is.na(df[c(Colreference)]) & !is.na(df$HXd) & (df[,Colvisit]-df[,Colreference]) < df$HXd
    df[Colvisitreferencediff,]$HXdlv<-df[Colvisitreferencediff,Colvisit ] - df[Colvisitreferencediff,Colreference ]
    
    
    df[ !is.na(df[c(Colvisit)]) & df$ERRtmp > df[,Coltrait] ,"ERR"]<-1
    df[ !is.na(df[c(Colvisit)]) ,"ERRtmp"]<-df[ !is.na(df[c(Colvisit)]) ,Coltrait] 
    
  }
  
  
  
  if( is.data.frame(df) && nrow(df)>0){
    DiagnosisVars=c("n_eid",paste(StrName,visitreference,sep="_"),"HXn","HXd","HXdlv","FUn","FUd","FUdlv","ERR") ### variables you wil create and you want to keep later... 
    DiagnosisVarsDescription=c("Identifier","Yes/No on VisitDate","History cases answered yes","History days before first observation","History to RefDate in days", "FollowUp cases answered yes","FollowUp from RefDate in days","FollowUp; last followup before first occurence in days","Erroneous individuals, questionable answers")
    df<-df[,DiagnosisVars] ## KEEP ONLY COLUMNS NEEDED.
    names(df)[-1]<-c(paste(StrTrait,"_",StrName,sep="") ,paste(StrTrait,"_",StrName,"_",names(df[-1:-2]),sep=""))
    
    ### ADD DESCRIPTION // RENAME VARIABLES 
   # attr(df, "var.labels") <- c("Identifier",paste(StrDescription,DiagnosisVarsDescription[-1],colapse="",sep=" "))
    attr(df, "var.labels") <- c("Identifier",paste(StrDescription,DiagnosisVarsDescription[-1],colapse="",sep=" "))
    attr(df, "var.labels")<-substr(attr(df, "var.labels"),0,80)
    save.dta13(df, StataOutputFile,compress = TRUE)
    return(df)
  } else {
    stop(paste("ERROR: could not find any variables"),paste(VctTSconditions,collapse=","))
  }
  
}

CheckIfColumnsAreAvailable<-function(df,VctColumnames){
  #TODO check if columns exists before running functions!!
  return(0)
}

Query_TSAmin_masterset<-function(dfmaster_TSDEATHMEDICD10_visitdtonly,StrTrait,StrDescription,StrName,VctTSAcolumns,visitdt,StataOutputFile){
  #dfmaster_TSDEATHMEDICD10_visitdtonly<-UKbioDatasetm
  #StrTrait<-"a1def"
  #StrDescription<-"asdsad"
  #VctTSAcolumns<-"n_22152_"
  #visitdt<-paste("ts_53_",2,"_0",sep="")
  #StataOutputFile<-"/data_work/databases/ukbiobanks/dataset/"
  ### CHECK IF FILE EXISTS: 
  if( file.exists(StataOutputFile )) {
    print(paste(StataOutputFile," ... already exists!, skipping"))
    return(0)     
  }
  print(paste(StrTrait,":",StrDescription))
  
  
  dfmaster_TSDEATHMEDICD10_visitdtonly_subset<-dfmaster_TSDEATHMEDICD10_visitdtonly[ ,  grepl( paste(c("n_eid",visitdt,VctTSAcolumns), collapse='|') , names( dfmaster_TSDEATHMEDICD10_visitdtonly ) ) ]
  
  ### BirthDate:
  dfmaster_TSDEATHMEDICD10_visitdtonly_subset$Birthdate<- as.Date(paste('01', paste(dfmaster_TSDEATHMEDICD10_visitdtonly$n_52_0_0,dfmaster_TSDEATHMEDICD10_visitdtonly$n_34_0_0)), format='%d %m %Y')
  dfmaster_TSDEATHMEDICD10_visitdtonly_subset$Age<-(dfmaster_TSDEATHMEDICD10_visitdtonly_subset[,visitdt] - dfmaster_TSDEATHMEDICD10_visitdtonly_subset$Birthdate) /365.25
  dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAOveralMin<-NA ## set initial value
  
  if( ncol(dfmaster_TSDEATHMEDICD10_visitdtonly_subset)==5) {stop("Missing TSA column? cannot find anything; so at least 1 missing")}  ## 0 column ; error; missing column. 
  
  
  
  for (VctTSAcolumn in VctTSAcolumns){
    print(paste("mean-minimum age across columns:",VctTSAcolumn))
    
    TSAcolumns_TF=grepl( paste(c(VctTSAcolumn), collapse='|') , names( dfmaster_TSDEATHMEDICD10_visitdtonly_subset ) )
    
    if(length(which(TSAcolumns_TF))==1){ ## if only one visit is available no mean is possible; then duplicate it. 
      print("Only one visit is available.")
      dfmaster_TSDEATHMEDICD10_visitdtonly_subset[paste(VctTSAcolumn,"tmp_0",sep="")]<-dfmaster_TSDEATHMEDICD10_visitdtonly_subset[TSAcolumns_TF]
      TSAcolumns_TF=grepl( paste(c(VctTSAcolumn), collapse='|') , names( dfmaster_TSDEATHMEDICD10_visitdtonly_subset ) )
    }
    
    dfmaster_TSDEATHMEDICD10_visitdtonly_subset[TSAcolumns_TF][dfmaster_TSDEATHMEDICD10_visitdtonly_subset[TSAcolumns_TF] == -1]<-NA 
    dfmaster_TSDEATHMEDICD10_visitdtonly_subset[TSAcolumns_TF][dfmaster_TSDEATHMEDICD10_visitdtonly_subset[TSAcolumns_TF] == -3]<-NA
    dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAmin<-rowMins(as.matrix( dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ ,grepl(paste(VctTSAcolumn, collapse='|'),names(dfmaster_TSDEATHMEDICD10_visitdtonly_subset)) ]),na.rm=TRUE)
    
    dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAmax<-rowMaxs(as.matrix( dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ ,grepl(paste(VctTSAcolumn, collapse='|'),names(dfmaster_TSDEATHMEDICD10_visitdtonly_subset)) ]),na.rm=TRUE)
    dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSADiff<-dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAmax-dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAmin
    dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAmean<-rowMeans(as.matrix( dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ ,grepl(paste(VctTSAcolumn, collapse='|'),names(dfmaster_TSDEATHMEDICD10_visitdtonly_subset))]),na.rm=TRUE)
    ### filter outliers: difference between visits >5years.
    outliers<-dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSADiff>5 & !is.infinite(dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSADiff)
    
    if(length(which(outliers))>0){## check if outliers exist. 
      print("removing TSAd outliers based on age-difference>5")
      dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ outliers,]$TSAmean<-NA
    } 
    
    #dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ is.infinite(dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAmean),]$TSAmean<-NA
    ### add minimum to existing column;
    dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAOveralMin<-rowMins(as.matrix( cbind(dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAmean,dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAOveralMin)),na.rm=TRUE)
    dfmaster_TSDEATHMEDICD10_visitdtonly_subset[is.infinite(dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAOveralMin),]$TSAOveralMin<-NA

  }
  
  dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd<-as.numeric(round((dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAOveralMin - dfmaster_TSDEATHMEDICD10_visitdtonly_subset$Age)*365.25))
  
  
  
  # dfmaster_TSDEATHMEDICD10_visitdtonly_subset<-dfmaster_TSDEATHMEDICD10_visitdtonly_subset[dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAmean>=0,] # filter age>0?
  dfmaster_TSDEATHMEDICD10_visitdtonly_subset<-(dfmaster_TSDEATHMEDICD10_visitdtonly_subset[!is.na(dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd),c("n_eid","TSAd")])
  
  ### Recalculate to Days from Reference. 
  dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ abs(dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd)<365.25 & dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd>0, "TSAd"] <- dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ abs(dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd)<365.25 & dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd>0, "TSAd"] /2 # within 1 year of diagnosis
  dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ abs(dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd)<365.25 & dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd<0, "TSAd"] <- dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ abs(dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd)<365.25 & dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd<0, "TSAd"] /2 # within 1 year of diagnosis
  dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd <= -365.25 ,"TSAd"]<-dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd  <= -365.25 ,"TSAd"]+round(365.25/2) ## beyond one year in history
  dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd >= 365.25 ,"TSAd"]<-dfmaster_TSDEATHMEDICD10_visitdtonly_subset[ dfmaster_TSDEATHMEDICD10_visitdtonly_subset$TSAd  >= 365.25 ,"TSAd"]-round(365.25/2) ## beyond one year in future
  
  if( is.data.frame(dfmaster_TSDEATHMEDICD10_visitdtonly_subset) && nrow(dfmaster_TSDEATHMEDICD10_visitdtonly_subset)>0){
    
    names(dfmaster_TSDEATHMEDICD10_visitdtonly_subset)<-c("n_eid",paste(StrTrait,StrName,sep="_"))
    attr(dfmaster_TSDEATHMEDICD10_visitdtonly_subset, "var.labels") <- c("Identifier",StrDescription)
    attr(dfmaster_TSDEATHMEDICD10_visitdtonly_subset, "var.labels")<-substr(attr(dfmaster_TSDEATHMEDICD10_visitdtonly_subset, "var.labels"),0,80)
    save.dta13(dfmaster_TSDEATHMEDICD10_visitdtonly_subset, StataOutputFile,compress = TRUE)
    return(dfmaster_TSDEATHMEDICD10_visitdtonly_subset)
  } else {
    print(paste("WARNING: df empty; could not find any variables OR observations (could be due to visit 2 having no observations vs visit 1), maybe check: ",StrTrait,"///",paste(VctTSAcolumn,collapse=",")))
  }
}


#' Extracting Vars From larger UK Bio dataset
#'
#' Extracts vars from ukbiobank from larger set of variables based on selection. 
#'
#' @param UKbioDataset df
#' @param dfDefinitions_VarSelection df
#' @param StataOutputFile file  
#' @keywords ExtractVarsFromMasterSet CreateUKBiobankPhentoypes
#' @return None
#'
#' @examples
#' ExtractVarsFromMasterSet(UKbioDataset,dfDefinitions_VarSelection,StataOutputFile)
#'
#' @export
  
ExtractVarsFromMasterSet<-function(UKbioDataset,dfDefinitions_VarSelection,StataOutputFile){
  
  #OutputdirMerged=paste(Outputdir,"/merged/",sep="")
  #StataOutputFile=paste(OutputdirMerged,"/DataSelection.dta",sep="")
  
  ############################################################
  ############ EXTRACT SELECTION FROM MASTERSET:
  dfSQL_TraitTable_VarSelection_Selected<-dfDefinitions_VarSelection[dfDefinitions_VarSelection$Include==1,]
  #as.vector(dfSQL_TraitTable_VarSelection_Selected$Var)
  ### FOUND:
  KeyAvail=names(UKbioDataset) %in% as.vector(dfSQL_TraitTable_VarSelection_Selected$var)
  Avail=(names(UKbioDataset)[KeyAvail])
  dfmasterSelection<-UKbioDataset[KeyAvail]
  attr(dfmasterSelection,"var.labels")<-attr(UKbioDataset,"var.labels")[KeyAvail]
  attr(dfmasterSelection, "var.labels")<-substr(attr(dfmasterSelection, "var.labels"),0,80)
  
  ### NOT FOUND
  NotAvail=as.vector(dfSQL_TraitTable_VarSelection_Selected$Var[!as.vector(dfSQL_TraitTable_VarSelection_Selected$var) %in% names(UKbioDataset) ])
  
  print(paste("FOUND:",paste(Avail,collapse=", ")))
  print(paste("MISSING:",paste(NotAvail,collapse=", ")))
  
  save.dta13(dfmasterSelection,StataOutputFile,compress = TRUE)
  
}
