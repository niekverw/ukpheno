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
#' @param VctOutputIndividualColumns set VctOutputIndividualColumns=c("TS","SR","TS_RX","SR_RX","LAB") to output all individual columns, default is 'VctOutputIndividualColumns=c()'
#' @keywords CreateUKBiobankPhentoypes
#' @references <Please contact mail@niekverweij.com, this package is not yet ready to be distributed.>
#' @importFrom parallel mclapply detectCores
#' @importFrom matrixStats rowMins
#' @importFrom matrixStats rowMaxs
#' @importFrom matrixStats rowSums2
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
#' CreateUKBiobankPhentoypes(Nvisits,visitreference,UKbioDataset,dfhesintables,dfDefinitions,Outputdir,VctOutputIndividualColumns=c(),VctAllUKBVDefinitionColumns=c("TS","SR","TS_RX","SR_RX","LAB") )
#'
#' #to output the variables based on RX (medication) in an *_RX file and exclude RX from the aggregate *_UKBV files run the following:
#' CreateUKBiobankPhentoypes(Nvisits,visitreference,UKbioDataset,dfhesintables,dfDefinitions,Outputdir,VctOutputIndividualColumns=c("SR_RX",VctAllUKBVDefinitionColumns=c("TS","SR","TS_RX","LAB") )
################################################################################################################################################
# library(readstata13)
# library(dplyr)
# library(data.table) #
# library(parallel) ## native to R.
# library(gdata)   ## for excel
# library(matrixStats) ## for rowMins/max function
###########################################################
CreateUKBiobankPhentoypes<-function(Nvisits,
                                    visitreference,
                                    dfmaster_TSDEATHMEDICD10_visitdtonly=UKbioDataset,
                                    dfmaster_SQL_merge,
                                    dfgpclinical,dfgpscripts,
                                    dfDefinitions,Outputdir,VctOutputIndividualColumns=c(),
                                    merge_all_output=TRUE ){

  #### DEFINE COLUMNS
  VctAllUKBVDefinitionColumns=c("TS","SR","TS_RX","SR_RX","LAB") #set this variable to a selection of columns (dfDefinition columns) to be outputted by the _UKBV variable, default is 'VctAllUKBVDefinitionColumns=c("TS","SR","TS_RX","SR_RX","LAB")'
  VctAllHESINDefinitionColumns=c("ICD10CODES","ICD9CODES","OPCS4CODES","OPCS3CODES")
  VctGPColumns=c("READCODES","CTV3CODES","BNFCODES","DMDCODES")
  StrTSAgeColumn="TS_AGE_DIAG_COLNAME"
  StrSRcolumns<-c("n_20001_","n_20002_","n_20003_","n_20004_")
  VctAllColumns=c("DEPENDENCY",VctAllUKBVDefinitionColumns,VctAllHESINDefinitionColumns,VctGPColumns,StrTSAgeColumn,StrSRcolumns)

  ### OPTIONAL SETTINGS
  ###parallel=2 <-need to fix.
  ################################################################################################################################################
  ###### RUN following loop:
  ################################################################################################################################################
  # dfmaster_SQL_merge<-dfhesintables
  # i=1
  # visitreference=0
  # VctOutputIndividualColumns=c("TS","SR","TS_RX","SR_RX","LAB")

  visitdt=paste("ts_53_",visitreference,"_0",sep="")

  print("processing dfDefinitions")
  dfDefinitions<-ProcessDfDefinitions(dfDefinitions,VctAllColumns)
  ### TODO: check if n_52_0_0 and n_34_0_0 are available == month/year of birth
  ### TODO: add function to check if columns in dfdefinitions are availabe.
  ### TODO: check how many hits per ICD10code --> logfile.
  ### TODO: evaluate what happens if uncommented and refvisit=2 where some ts_53 are missing.
  #dfmaster_TSDEATHMEDICD10_visitdtonly<-UKbioDataset[!is.na(UKbioDataset[visitdt]),]
  #dfmaster_TSDEATHMEDICD10_visitdtonly<-dfmaster_TSDEATHMEDICD10_visitdtonly

  # running a test on cholesterol Touchscreen, shouldnt be characters..
  #if(nchar(as.character(unique(dfmaster_TSDEATHMEDICD10_visitdtonly[,"n_6153_0_0"]))[2])>2){print("ERROR: be aware that column entries are text, convert to numbers, defactor.."); return(0)}


  # write n_eids available from gp and hesin data
  dir.create( Outputdir, showWarnings =F,recursive=T)
  fwrite(data.frame(n_eid=unique(c(dfgpclinical$n_eid,dfgpscripts$n_eid))),paste0(Outputdir,"/n_eids.GP.txt"))
  fwrite(data.frame(n_eid=unique(c(dfmaster_SQL_merge$n_eid))),paste0(Outputdir,"/n_eids.HESIN.txt"))
  fwrite(data.frame(n_eid=unique(c(dfmaster_TSDEATHMEDICD10_visitdtonly$n_eid))),paste0(Outputdir,"/n_eids.ALL.txt"))

  print("checking definitions; for(i in 1:nrow(dfDefinitions))")
  for(i in 1:nrow(dfDefinitions)) {
    row <- dfDefinitions[i,]
    row <-ConvertFactorsToStringReplaceNAInDf(row) ## convert factor->string
    row$DESCRIPTION<-gsub(pattern = '[^a-zA-Z0-9_ ()=|;,+-\\&]+',"",row$DESCRIPTION) ## remove non-alphanumeric except ()=|;,+-\\&
    StrTrait <- row$TRAIT
    epidurfilter <- max(0,row$Minimum_Episode_duration,na.rm = T) # 0 by default in case its not filled out.
    include_secondary <- min(1,row$include_secondary,na.rm = T) # 1 by default in case its not filled out.

    if (grepl("_",StrTrait)) {
      stop(paste("Trait: ",StrTrait," contains '_' !!,",sep=""))

    }

    ########################################################
    ########################################################
    ### HESIN DATA : exact dates.
    Strcatagory="HESIN"
    ########################################################

    ########################################################
    ####### OPERATION:"OPCS4CODES","OPCS3CODES"
    #######################################################
    if( !is.na(row$OPCS4CODES)  ){
      #### SETTINGS:
      print(paste("   ..finding operationcodes",StrTrait))
      VctCodes<-unlist(strsplit(row$OPCS4CODES,","))


      StrColumnForHescodes<-c("oper4_1","oper4_2" )
      if(include_secondary==0) {
        StrColumnForHescodes<-c("oper4_1")
      }
      StrColumnnameForEventDate<-"opdate_1"
      ####################################
      ### HES: FOLLOW UP + HISTORY VARIABLES:
      StrName="O4hes"
      StrColumnnameForVisitDate<-visitdt
      StrDescription<-paste(row$DESCRIPTION,"- OPCS4(HES)",paste(VctCodes,collapse=","))

      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      Outcome_HES(dfmaster_SQL_merge,paste(StrTrait,StrName,sep="_"),StrDescription,VctCodes,epidurfilter,StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile)

    }


    if( !is.na(row$OPCS3CODES)  ){
      #### SETTINGS:
      print(paste("   ..finding operationcodes",StrTrait))
      VctCodes<-unlist(strsplit(row$OPCS3CODES,","))


      StrColumnForHescodes<-c("oper3_1","oper3_2" )
      if(include_secondary==0) {
        StrColumnForHescodes<-c("oper3_1")
      }
      StrColumnnameForEventDate<-"opdate_1"
      ####################################
      ### HES: FOLLOW UP + HISTORY VARIABLES:
      StrName="O3hes"
      StrColumnnameForVisitDate<-visitdt
      StrDescription<-paste(row$DESCRIPTION,"- OPCS3(HES)",paste(VctCodes,collapse=","))

      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      Outcome_HES(dfmaster_SQL_merge,paste(StrTrait,StrName,sep="_"),StrDescription,VctCodes,epidurfilter,StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile)

    }
    ########################################################
    ####### ICD10 LOOKUPS:
    ########################################################
    if( !is.na(row$ICD10CODES) ) {
      #### SETTINGS:
      print("    ..finding icd10diagnosiscodes")
      VctCodes<-unlist(strsplit(row$ICD10CODE,","))
      ####################################
      ### HESIN ICD10: FOLLOW UP + HISTORY VARIABLES:
      StrName="Dhes"

      StrColumnForHescodes<-c("diag_icd10_1","diag_icd10_2" )
      if(include_secondary==0) {
        StrColumnForHescodes<-c("diag_icd10_1")
      }
      StrColumnnameForEventDate<-"epistart_1"
      StrDescription<-paste(row$DESCRIPTION,"- ICD10(HES) -",paste(VctCodes,collapse=","))

      StrColumnnameForVisitDate<-visitdt

      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      Outcome_HES(dfmaster_SQL_merge,StrTrait=paste(StrTrait,StrName,sep="_"),StrDescription,VctCodes,epidurfilter=epidurfilter,StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile)

      ####################################
      #### DEATH ICD10:
      StrName<-"DO"
      StrFieldcode<-"40001|40002"
      StrFieldDateColumn<-"ts_40000_0_0"
      StrDescription<-paste(row$DESCRIPTION," - ","Prim+sec Death Of Disorder ICD10(DB):",paste(VctCodes,collapse=","),sep="")

      if(include_secondary==0) {
        StrFieldcode<-"40001"
        StrDescription<-paste(row$DESCRIPTION," - ","Prim!-!Sec Death Of Disorder ICD10(DB):",paste(VctCodes,collapse=","),sep="")
      }

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
      if(include_secondary==0) {
        StrColumnForHescodes<-c("diag_icd9_1" )
      }

      StrColumnnameForEventDate<-"epistart_1"
      StrDescription<-paste(row$DESCRIPTION,"- ICD9(HES) -",paste(VctCodes,collapse=","))
      StrColumnnameForVisitDate<-visitdt

      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      Outcome_HES(dfmaster_SQL_merge,paste(StrTrait,StrName,sep="_"),StrDescription,VctCodes,epidurfilter,StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile)

    }

    ########################################################
    ####### TS(touchscreen//SR(selfreported) // RX // LAB: only based on visit dates. ..
    ########################################################
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
              StrDescription<-paste(row$DESCRIPTION,"-",StrDefinitionColumn)
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


    ########################################################
    ####### READ CODES (V2 + V3): GP DATA; CTV3CODES, READCODES
    ########################################################
    Strcatagory = "GP"
    if( !is.na(row$READCODES)  ) {
      #### SETTINGS:
      print("    ..finding READCODES diagnosiscodes")
      VctCodes<-unlist(strsplit(row$READCODES,","))
      ####################################
      ### READ: FOLLOW UP + HISTORY VARIABLES:
      ###
      StrName="Read"
      StrColumnForHescodes<-c("read_2")
      StrColumnnameForEventDate<-"event_dt"
      StrDescription<-paste(row$DESCRIPTION,"- READ -",paste(VctCodes,collapse=","))
      StrColumnnameForVisitDate<-visitdt

      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      # fix epiend_1 is not available in GP, so added it to the dataframe.. but it is not useful here.
      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      Outcome_HES(dfmaster_SQL_merge = dfgpclinical,
                  StrTrait = paste(StrTrait,StrName,sep="_"),
                  StrDescription,VctCodes,epidurfilter,
                  StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile)

    }

    Strcatagory = "GP"
    if( !is.na(row$CTV3CODES)  ) {
      #### SETTINGS:
      print("    ..finding READCODES diagnosiscodes")
      VctCodes<-unlist(strsplit(row$CTV3CODES,","))
      ####################################
      ### READ: FOLLOW UP + HISTORY VARIABLES:
      ###
      StrName="CTV3"
      StrColumnForHescodes<-c("read_3")
      StrColumnnameForEventDate<-"event_dt"
      StrDescription<-paste(row$DESCRIPTION,"- CTV3 -",paste(VctCodes,collapse=","))
      StrColumnnameForVisitDate<-visitdt

      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      # fix epiend_1 is not available in GP, so added it to the dataframe.. but it is not useful here.
      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      Outcome_HES(dfmaster_SQL_merge = dfgpclinical,
                  StrTrait = paste(StrTrait,StrName,sep="_"),
                  StrDescription,VctCodes,epidurfilter,
                  StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile)

    }
    ########################################################
    ####### MEDICATION FROM READCODES +  DMDCODES + BNFCODES
    ########################################################
    Strcatagory = "GP"
    if( !is.na(row$BNFCODES)  ) {
      #### SETTINGS:
      print("    ..finding READCODES diagnosiscodes")
      VctCodes<-unlist(strsplit(row$BNFCODES,","))
      ####################################
      ### READ: FOLLOW UP + HISTORY VARIABLES:
      ###
      StrName="BNF"
      StrColumnForHescodes<-c("bnf_code")
      StrColumnnameForEventDate<-"event_dt"
      StrDescription<-paste(row$DESCRIPTION,"- BNF -",paste(VctCodes,collapse=","))
      StrColumnnameForVisitDate<-visitdt

      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      # fix epiend_1 is not available in GP, so added it to the dataframe.. but it is not useful here.
      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      Outcome_HES(dfmaster_SQL_merge = dfgpclinical,
                  StrTrait = paste(StrTrait,StrName,sep="_"),
                  StrDescription,VctCodes,epidurfilter,
                  StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile)

    }


    Strcatagory = "GP"
    if( !is.na(row$DMDCODES)  ) {
      #### SETTINGS:
      print("    ..finding READCODES diagnosiscodes")
      VctCodes<-unlist(strsplit(row$DMDCODES,","))
      ####################################
      ### READ: FOLLOW UP + HISTORY VARIABLES:
      ###
      StrName="DMD"
      StrColumnForHescodes<-c("dmd_code")
      StrColumnnameForEventDate<-"event_dt"
      StrDescription<-paste(row$DESCRIPTION,"- DMD -",paste(VctCodes,collapse=","))
      StrColumnnameForVisitDate<-visitdt

      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
      # fix epiend_1 is not available in GP, so added it to the dataframe.. but it is not useful here.
      dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
      Outcome_HES(dfmaster_SQL_merge = dfgpclinical,
                  StrTrait = paste(StrTrait,StrName,sep="_"),
                  StrDescription,VctCodes,epidurfilter,
                  StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile)

    }
    ##########################################
    ### MERGING Files:
    ##########################################

      ##########################################
      ###### MERGING ICD10/9 and death follow up.
      ##########################################
    print(visitdt)
    filesEpisode<-list.files(path=paste(Outputdir,"/",visitdt,"/HESIN",sep=""),pattern=paste("^",StrTrait,"(_)(.*)dta$",sep=""), full.names=TRUE,recursive = F)
    if( length(filesEpisode)>0) {
      HEScodes=paste(row$ICD10CODES,row$ICD9CODES,row$OPCS4CODES,row$OPCS3CODES,sep="|")
      OutputdirMerged=paste(Outputdir,"/",visitdt,"/HESIN/merged/",sep="")
      StataOutputFile=paste(OutputdirMerged,"/",StrTrait,"_HESIN.dta",sep="")

      if( file.exists(StataOutputFile )) {### CHECK IF FILE EXISTS:
        print(paste(StataOutputFile," ... already exists!, skipping"))
      } else {
      dfMerged<-MultiMergeEpisodeData(filesEpisode,StrTrait,row$DESCRIPTION,HEScodes)
      dir.create(OutputdirMerged, showWarnings = FALSE, recursive = TRUE)
      if(is.data.frame(dfMerged)){      save.dta13(dfMerged,StataOutputFile,compress = TRUE) }
      }
    }

      ##########################################
      ## > TODO MERGE ALL GP data
      ##########################################

    print(visitdt)
    filesEpisode<-list.files(path=paste(Outputdir,"/",visitdt,"/GP",sep=""),pattern=paste("^",StrTrait,"(_)(.*)dta$",sep=""), full.names=TRUE,recursive = F)
    if( length(filesEpisode)>0) {
      HEScodes=paste(row$READCODES,row$CTV3CODES,row$BNFCODES,row$BNFCODES,row$DMDCODES,sep="|")
      OutputdirMerged=paste(Outputdir,"/",visitdt,"/GP/merged/",sep="")
      StataOutputFile=paste(OutputdirMerged,"/",StrTrait,"_GP.dta",sep="")

      if( file.exists(StataOutputFile )) {### CHECK IF FILE EXISTS:
        print(paste(StataOutputFile," ... already exists!, skipping"))
      } else {
        dfMerged<-MultiMergeEpisodeData(filesEpisode,StrTrait,row$DESCRIPTION,HEScodes)
        dir.create(OutputdirMerged, showWarnings = FALSE, recursive = TRUE)
        if(is.data.frame(dfMerged)){      save.dta13(dfMerged,StataOutputFile,compress = TRUE) }
      }
    }
      ##########################################
      ###### MERGING ALL EPISODE DATA
      ##########################################
    print(visitdt)
    filesEpisode<-list.files(path=paste(Outputdir,"/",visitdt,"/HESIN",sep=""),pattern=paste("^",StrTrait,"(_)(.*)dta$",sep=""), full.names=TRUE,recursive = F)
    filesEpisode<-c(filesEpisode,list.files(path=paste(Outputdir,"/",visitdt,"/GP",sep=""),pattern=paste("^",StrTrait,"(_)(.*)dta$",sep=""), full.names=TRUE,recursive = F))

    if( length(filesEpisode)>0) {
      EPcodes=paste(row$ICD10CODES,row$ICD9CODES,row$OPCS4CODES,row$OPCS3CODES,row$READCODES,row$CTV3CODES,row$BNFCODES,row$DMDCODES,sep="|")
      OutputdirMerged=paste(Outputdir,"/",visitdt,"/ALL/Episodes/",sep="")
      dir.create(OutputdirMerged, showWarnings = FALSE, recursive = TRUE)
      StataOutputFile=paste(OutputdirMerged,"/",StrTrait,"_merged.dta",sep="")

      if( file.exists(StataOutputFile )) {### CHECK IF FILE EXISTS:
        print(paste(StataOutputFile," ... already exists!, skipping"))
      } else {
        dfMerged<-MultiMergeEpisodeData(filesEpisode,StrTrait,row$DESCRIPTION,EPcodes,suffix="_EP_")
        dir.create(OutputdirMerged, showWarnings = FALSE, recursive = TRUE)
        if(is.data.frame(dfMerged)){      save.dta13(dfMerged,StataOutputFile,compress = TRUE) }
      }
    }



    ##########################################
    #### MERGING EP merged , UKBV and TSAmin
    fileEPmerged<-list.files(path=paste(Outputdir,"/",visitdt,"/ALL/Episodes/",sep=""),pattern=paste("^",StrTrait,"(_)(.*)dta$",sep=""), full.names=TRUE)
    fileUKBV<-list.files(path=paste(Outputdir,"/",visitdt,"/UKBVisit/ALL",sep=""),pattern=paste("^",StrTrait,"(_)(UKBV.*)dta$",sep=""), full.names=TRUE)
    fileTSA<-list.files(path=paste(Outputdir,"/",visitdt,"/UKBVisit/",sep=""),pattern=paste("^",StrTrait,"(_)(TSA.*)dta$",sep=""), full.names=TRUE)
    files=c(fileEPmerged,fileUKBV,fileTSA)
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

  gc()
  ## summarize sources (death, hesin, gp (read), self reports from nurse (SR), touchscreen (TS) )
  df_phenotype_crosstabulation <- summarize_cross_phenotype_fields(paste0(Outputdir,"/",visitdt),UKbioDataset$n_eid) #dfmaster_TSDEATHMEDICD10_visitdtonly
  fwrite(df_phenotype_crosstabulation,f=paste0(Outputdir,"/",visitdt,"/crosstabulation.tsv"),row.names = F,sep = "\t")

  # merge traits
  if (merge_all_output){
    merge_stata_files(statafiles = list.files(paste0(Outputdir,"/",visitdt,"/ALL/"),full.names = T,pattern = "*.dta$" ),paste0(Outputdir,"/",visitdt,"/ALL.dta"))
  }


}




