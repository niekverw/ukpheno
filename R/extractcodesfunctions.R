

#' @export
Outcome_HES<-function(dfmaster_SQL_merge,StrTrait,StrDescription,VctCodes,epidurfilter=0,StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile){
  # epidurfilter=0
  # dfmaster_SQL_merge=dfgpclinical
  # dfmaster_SQL_merge=dfhesintables
  ### CHECK IF FILE EXISTS:
  if( file.exists(StataOutputFile )) {
    print(paste(StataOutputFile," ... already exists!, skipping"))
    return(0)
  }
  print(paste("extracting outcomes from episode data, checking codes:",StrTrait,":",StrDescription,"-",StrColumnnameForVisitDate, "--",paste(StrColumnForHescodes,collapse=",") ) )

  DiagnosisVars=c("n_eid","HXn","HXd","HXt","HXto","FUn","FUd","FUt","FUto") ### variables you wil create and you want to keep later...
  ### FUt/ FUto = NOT FOR GP DATA, only HESIN
  #### RENAME COLUMNS
  names(dfmaster_SQL_merge)[names(dfmaster_SQL_merge)==StrColumnnameForEventDate]<-"event_date"
  names(dfmaster_SQL_merge)[names(dfmaster_SQL_merge)==StrColumnnameForVisitDate]<-"visit_date"
  ## remove records that have no visit date (relevant for  ts_53_1_0 ,  ts_53_2_0).
  dfmaster_SQL_merge_tmp<-dfmaster_SQL_merge[!is.na(dfmaster_SQL_merge$visit_date),]

  ########### GREP CODE AT THE COLUMNS 'StrColumnForHescodes' ;; TODO grep should be ^: paste(sep="","^",VctCodes, collapse='|')
  #grepoper=c()
  #for (i in StrColumnForHescodes ){
  #  grepoper<-c(grepoper,grep(paste(VctCodes, collapse='|'), dfmaster_SQL_merge_tmp[,i], ignore.case=FALSE))
  #}
  if (length(StrColumnForHescodes)==1){StrColumnForHescodes = c(StrColumnForHescodes,StrColumnForHescodes)}
  grepoper <-unlist( mclapply(  dfmaster_SQL_merge_tmp[, c(StrColumnForHescodes,StrColumnForHescodes)], function(col) grep(paste(sep="","^",VctCodes, collapse='|'), col, ignore.case=FALSE),mc.cores =detectCores()/2 ) ) ## PARALLEL of the above.

  if(length(grepoper) ==0) {print("    nothing here");return(0) }

  dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_tmp[ unique(c(grepoper)),]
  ## filter Episode Durations on a minimum duration beforehand
  if (epidurfilter>0){
    dfmaster_SQL_merge_oper     <- dfmaster_SQL_merge_oper[abs(round(difftime( dfmaster_SQL_merge_oper$epiend_1, dfmaster_SQL_merge_oper$event_date ,units="days") )) >= epidurfilter,]
  }
    # remove duplicates
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
  ## Episode Duration.
  dfmaster_SQL_merge_oper$HXt <-ifelse(dfmaster_SQL_merge_oper$visit_date >= dfmaster_SQL_merge_oper$event_date, round(difftime( dfmaster_SQL_merge_oper$event_date , dfmaster_SQL_merge_oper$epiend_1 ,units="days")) , NA) ### SET FUP time (days) - first event in history
  dfmaster_SQL_merge_oper     <-suppressWarnings(  dfmaster_SQL_merge_oper %>%  group_by(n_eid) %>%  mutate(HXt = -max(abs(HXt), na.rm=T) ) ) ### replace HXd
  dfmaster_SQL_merge_oper$HXt[is.infinite(dfmaster_SQL_merge_oper$HXt)]<-NA

  # ## total episode dur (tt): ### NOT FOR GP DATA
  dfmaster_SQL_merge_oper$episodedur <- abs(round(difftime( dfmaster_SQL_merge_oper$event_date , dfmaster_SQL_merge_oper$epiend_1 ,units="days")))
  dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_oper %>% group_by(n_eid) %>% mutate(HXto = sum(ifelse (visit_date >= event_date, episodedur ,0 ) ))
  #
  ## follow up vars; 3170028 1834566
  dfmaster_SQL_merge_oper$FU  <-ifelse(dfmaster_SQL_merge_oper$visit_date < dfmaster_SQL_merge_oper$event_date, 1, 0)   ### SET FOLLOW UP VARIABLE (0/1)
  dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_oper %>% group_by(n_eid) %>% mutate(FUn =sum(FU)  )  ### Count number of followup events
  dfmaster_SQL_merge_oper     <- suppressWarnings(dfmaster_SQL_merge_oper %>%  group_by(n_eid) %>%  mutate(FU = max(FU) ))  ### replace FU for every duplicating n_eid
  dfmaster_SQL_merge_oper$FUd <-ifelse(dfmaster_SQL_merge_oper$visit_date < dfmaster_SQL_merge_oper$event_date, difftime(dfmaster_SQL_merge_oper$event_date, dfmaster_SQL_merge_oper$visit_date ,units="days") , NA) ### SET FUP time (days) - first event in future.
  dfmaster_SQL_merge_oper     <-suppressWarnings( dfmaster_SQL_merge_oper %>%  group_by(n_eid) %>%  mutate(FUd = min(FUd, na.rm=T) ) ) ### replace HXd
  dfmaster_SQL_merge_oper$FUd[is.infinite(dfmaster_SQL_merge_oper$FUd)]<-NA
  #  dfmaster_SQL_merge_oper$FUd1 <-ifelse(dfmaster_SQL_merge_oper$visit_date < dfmaster_SQL_merge_oper$event_date, difftime(dfmaster_SQL_merge_oper$event_date, dfmaster_SQL_merge_oper$visit_date ,units="days") , NA) ### SET FUP time (days)  - last event in future.
  #  dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_oper %>%  group_by(n_eid) %>%  mutate(FUd1 = max(FUd1, na.rm=T) )  ### replace HXd
  dfmaster_SQL_merge_oper$FUt <-ifelse(dfmaster_SQL_merge_oper$visit_date < dfmaster_SQL_merge_oper$event_date, round(difftime( dfmaster_SQL_merge_oper$epiend_1, dfmaster_SQL_merge_oper$event_date ,units="days") ) , NA) ### SET FUP time (days) - first event in future.
  dfmaster_SQL_merge_oper     <-suppressWarnings(dfmaster_SQL_merge_oper %>%  group_by(n_eid) %>%  mutate(FUt = min(FUt, na.rm=T) ))  ### replace HXd
  # ## total episode dur (tt):
  dfmaster_SQL_merge_oper     <-dfmaster_SQL_merge_oper %>% group_by(n_eid) %>% mutate(FUto = sum(ifelse (visit_date < event_date, episodedur ,0 ) ))


  ## keep uniq
  dfmaster_SQL_merge_oper <- dfmaster_SQL_merge_oper[!duplicated(dfmaster_SQL_merge_oper[,DiagnosisVars]),][,DiagnosisVars] ## FILTER UNIQ + KEEP ONLY COLUMNS NEEDED.

  print(paste("   #found ",nrow(dfmaster_SQL_merge_oper) ))

  ### ADD DESCRIPTION // RENAME VARIABLES
  attr(dfmaster_SQL_merge_oper, "var.labels") <- c("Identifier",paste(StrDescription, colnames(dfmaster_SQL_merge_oper[,-1]), sep = " - ") )
  attr(dfmaster_SQL_merge_oper, "var.labels")<-substr(attr(dfmaster_SQL_merge_oper, "var.labels"),0,80)
  colnames(dfmaster_SQL_merge_oper)[-1] <- paste(StrTrait, colnames(dfmaster_SQL_merge_oper[,-1]), sep = "_")
  save.dta13(dfmaster_SQL_merge_oper, StataOutputFile,compress = TRUE)
  return(dfmaster_SQL_merge_oper)
}

#' @export
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
  #  grepoper<-c(grepoper,grep(paste(VctCodes, collapse='|'), dfmaster_TSDEATHMEDICD10_visitdtonly[,i], ignore.case=FALSE))
  grepoper <-unlist( mclapply(  dfmaster_TSDEATHMEDICD10_visitdtonly[, StrColumnNames], function(col) grep(paste(VctCodes, collapse='|'), col, ignore.case=FALSE),mc.cores =detectCores()/2 ) ) ### PARALLEL., the abovve is not parallelized
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


#' @export
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
      # Condition=VctTSconditions[grep("=",VctTSconditions)][1]
      print( unlist(strsplit(Condition,split="=")) )
      StrTScolumn=unlist(strsplit(Condition,split="="))[1]
      StrTScodes=unlist(strsplit(Condition,split="="))[2]
      StrColumnNames<-names(dfmaster_TSDEATHMEDICD10_visitdtonly[ grepl(StrTScolumn, names(dfmaster_TSDEATHMEDICD10_visitdtonly)) ])
      StrColumnNames_pervisit<-StrColumnNames[grepl(visitcol,StrColumnNames)]
      if ( length(StrColumnNames) == 0)  { stop(paste("ERROR: at least one column missing, there may be more missing, please check!!!: ",paste(StrTScolumn,collapse=",") ) ) }

      for (i in StrColumnNames_pervisit ) {
        # i = StrColumnNames_pervisit[1]
        #for (code in unlist(strsplit(StrTScodes,split="\\|" )) ){
        #  grepoper<- c(grepoper,which(as.numeric(dfmaster_TSDEATHMEDICD10_visitdtonly[,i])==code))
        #}
        print(paste("   #visit=",visit, "; #Condition=",Condition, "; #i=",i, "; #found ",sum(as.numeric(dfmaster_TSDEATHMEDICD10_visitdtonly[,i]) %in% unlist(strsplit(StrTScodes,split="\\|" )))) )
        grepoper<- c(grepoper,which(as.numeric(dfmaster_TSDEATHMEDICD10_visitdtonly[,i]) %in% unlist(strsplit(StrTScodes,split="\\|" )) ))  ## faster.. but the above is similar to the code below.. should figure out how to merge below with this or something..
      }
    }
    # ≥ argument  can be _n3434_≥1   , FIX ≤

    for (Condition in VctTSconditions[ grep("≥",VctTSconditions) ]) {
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

#' @export
CheckIfColumnsAreAvailable<-function(df,VctColumnames){
  #TODO check if columns exists before running functions!!
  return(0)
}


#' @export
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


#' #' Extracting Vars From larger UK Bio dataset
#' #'
#' #' Extracts vars from ukbiobank from larger set of variables based on selection.
#' #'
#' #' @param UKbioDataset df
#' #' @param dfDefinitions_VarSelection df
#' #' @param StataOutputFile file
#' #' @keywords ExtractVarsFromMasterSet CreateUKBiobankPhentoypes
#' #' @return None
#' #'
#' #' @examples
#' #' ExtractVarsFromMasterSet(UKbioDataset,dfDefinitions_VarSelection,StataOutputFile)
#' #'
#' #' @export
#' ExtractVarsFromMasterSet<-function(UKbioDataset,dfDefinitions_VarSelection,StataOutputFile){
#'
#'   #OutputdirMerged=paste(Outputdir,"/merged/",sep="")
#'   #StataOutputFile=paste(OutputdirMerged,"/DataSelection.dta",sep="")
#'
#'   ############################################################
#'   ############ EXTRACT SELECTION FROM MASTERSET:
#'   dfSQL_TraitTable_VarSelection_Selected<-dfDefinitions_VarSelection[dfDefinitions_VarSelection$Include==1,]
#'   #as.vector(dfSQL_TraitTable_VarSelection_Selected$Var)
#'   ### FOUND:
#'   KeyAvail=names(UKbioDataset) %in% as.vector(dfSQL_TraitTable_VarSelection_Selected$var)
#'   Avail=(names(UKbioDataset)[KeyAvail])
#'   dfmasterSelection<-UKbioDataset[KeyAvail]
#'   attr(dfmasterSelection,"var.labels")<-attr(UKbioDataset,"var.labels")[KeyAvail]
#'   attr(dfmasterSelection, "var.labels")<-substr(attr(dfmasterSelection, "var.labels"),0,80)
#'
#'   ### NOT FOUND
#'   NotAvail=as.vector(dfSQL_TraitTable_VarSelection_Selected$Var[!as.vector(dfSQL_TraitTable_VarSelection_Selected$var) %in% names(UKbioDataset) ])
#'
#'   print(paste("FOUND:",paste(Avail,collapse=", ")))
#'   print(paste("MISSING:",paste(NotAvail,collapse=", ")))
#'
#'   save.dta13(dfmasterSelection,StataOutputFile,compress = TRUE)
#'
#' }
