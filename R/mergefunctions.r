#library(matrixStats)
#' @export
min_nv<-function(x){ ## select value in the last half of x, on the basis of the first half of x
  x.lookup=x[1:(length(x)/2) ]
  x.return=x[(1+(length(x)/2)):(length(x)) ]
  #m <- suppressWarnings(min(x,na.rm=TRUE))
  i.m <- which.min(x.lookup)
  m <- x.return[i.m]
  if (length(m)==0){
    m=NA
  }
  return(m)
}



### MERGING FUNCTION:
#' @export
MultiMergeEpisodeData <- function(filesEpisode,StrTrait,StrDescription,HEScodes,suffix="_HES_"){
  #  filesEpisode<-c("/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/HESIN/CAD_D9hes.dta", "/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/HESIN/CAD_Dhes.dta", "/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/HESIN/CAD_DO.dta", "/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/HESIN/CAD_Ohes.dta")
  # filesEpisode<-c("/data_work/databases/ukbiobanks/Phenotypes/CAD_definitions_testNewscript/output/ts_53_0_0/HESIN/AF_D9hes.dta","/data_work/databases/ukbiobanks/Phenotypes/CAD_definitions_testNewscript/output/ts_53_0_0/HESIN/AF_Dhes.dta","/data_work/databases/ukbiobanks/Phenotypes/CAD_definitions_testNewscript/output/ts_53_0_0/HESIN/AF_DO.dta","/data_work/databases/ukbiobanks/Phenotypes/CAD_definitions_testNewscript/output/ts_53_0_0/HESIN/AF_DOp.dta")
  #filesEpisode<-c("/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_0_0/HESIN/CVD_Dhes.dta","/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_0_0/HESIN/CVD_DO.dta","/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_0_0/HESIN/CVD_DOp.dta")
  #  filesEpisode<-c("/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/HESIN/HxCb_D9hes.dta","/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/HESIN/HxCb_Dhes.dta")
  #  StrTrait<-"Ca"
  #  StrDescription<-"asd"
  #  HEScodes<-"Z825"
  # filesEpisode<-c("/data_work/databases/ukbiobank/9628_pvdh/HFmortality/ts_53_0_0/HESIN/HfMortHfinclQ_D9hes.dta",
  #  "/data_work/databases/ukbiobank/9628_pvdh/HFmortality/ts_53_0_0/HESIN/HfMortHfinclQ_Dhes.dta",
  #  "/data_work/databases/ukbiobank/9628_pvdh/HFmortality/ts_53_0_0/HESIN/HfMortHfinclQ_DO.dta",
  #  "/data_work/databases/ukbiobank/9628_pvdh/HFmortality/ts_53_0_0/HESIN/HfMortHfinclQ_DOp.dta")
  # suffix="_HES_"
  print(paste("merging",filesEpisode))

  LstDf = lapply(filesEpisode, function(x){    as.data.frame(read.dta13(x,convert.dates = TRUE))   })
  VctDescriptions<-unlist(lapply(LstDf,function(x){ varlabel(x) } ))
  dfVctDescriptions <- unique(cbind(names(VctDescriptions),unname(VctDescriptions)))

  dfMerged<-Reduce(function(x,y) {merge(x,y,by="n_eid",all=TRUE)}, LstDf)

  ### temporary rename primary cause of death vars: CAD2_DOp_FUn CAD2_DOp_FUd so they are not taken into account in the following loop
  VctTmpColNamesDOp<-names(dfMerged)[ grepl( "_DOp_FUd" , names( dfMerged ) ) | grepl( "_DOp_FUn" , names( dfMerged ) ) ]
  names(dfMerged)[ grepl( "_DOp_FUd" , names( dfMerged ) ) | grepl( "_DOp_FUn" , names( dfMerged ) ) ]<-c("DOpn","DOpd")


  ###################################
  ##### Combine HES:
  ###HX
  if(nrow(dfMerged)==0){
    print("no data found, returning 0 HES, exiting merge function")
    return(0)
  }

  if ( length(which(grepl( "_HXn" , names( dfMerged ) )))>1){
    dfMerged$HXn<-rowSums(dfMerged[ , grepl( "_HXn$" , names( dfMerged ) ) ] , na.rm=T)
    df.1 <- as.matrix(dfMerged[ , grepl( "_HXd$" , names( dfMerged ) ) ])
    HXd <- rowMins( df.1 , na.rm=T)
    HXd[HXd==Inf]<-NA
    dfMerged$HXd<-HXd

    df.2 <- as.matrix(dfMerged[ , grepl( "_HXt$" , names( dfMerged ) ) ])
    HXt <- apply(cbind(df.1,df.2), 1, FUN=min_nv)
    dfMerged$HXt<-HXt

    dfMerged$HXto <- rowSums2(as.matrix(dfMerged[ , grepl( "_HXto$" , names( dfMerged ) ) ]),na.rm=T)

  } else if (length(which(grepl( "_HXn$" , names( dfMerged ) ))) == 1) {
    dfMerged$HXn<-dfMerged[ , grepl( "_HXn$" , names( dfMerged ) ) ]
    dfMerged$HXd<-dfMerged[ , grepl( "_HXd$" , names( dfMerged ) ) ]
    dfMerged$HXt<-dfMerged[ , grepl( "_HXt$" , names( dfMerged ) ) ]
    dfMerged$HXto<-dfMerged[ , grepl( "_HXto$" , names( dfMerged ) ) ]
  }
  ### FU
  if (length(which(grepl( "_FUn" , names( dfMerged ) )))>1){

    dfMerged$FUn<-rowSums(dfMerged[ , grepl( "_FUn$" , names( dfMerged ) ) ] , na.rm=T)
    df.1 <- dfMerged[ , grepl( "_FUd$" , names( dfMerged )) ]
    FUd<-rowMins( as.matrix(df.1) , na.rm=T)
    FUd[FUd==Inf]<-NA
    dfMerged$FUd<-FUd


    df.1 <- df.1[ ,!grepl( "DO_FUd$" , names( df.1 )) ]
    df.2 <- as.matrix(dfMerged[ , grepl( "_FUt$" , names( dfMerged ) ) ]) # select HXt using minimum HXd
    FUt <- apply(cbind(df.1,df.2), 1, FUN=min_nv)
    dfMerged$FUt<-FUt

    dfMerged$FUto <- rowSums2(as.matrix(dfMerged[ , grepl( "_FUto$" , names( dfMerged ) ) ]),na.rm=T)

  } else if (length(which(grepl( "_FUn$" , names( dfMerged ) ))) == 1) {
    dfMerged$FUn<-dfMerged[ , grepl( "_FUn$" , names( dfMerged ) ) ]
    dfMerged$FUd<-dfMerged[ , grepl( "_FUd$" , names( dfMerged ) ) ]
    dfMerged$FUt<-dfMerged[ , grepl( "_FUt$" , names( dfMerged ) ) ]
    dfMerged$FUt<-dfMerged[ , grepl( "_FUtot$" , names( dfMerged ) ) ]
  }
  ########################################

  #rename back primary DO.

  names(dfMerged)[names(dfMerged) %in% c("DOpn","DOpd")]<-VctTmpColNamesDOp

  ### Add _HES to the variable name of HX/FU
  newcols = c("HXn","HXd","HXt","HXto","FUn","FUd","FUt","FUto")
  names(dfMerged)[names(dfMerged) %in% newcols] <- paste0(StrTrait,suffix,names(dfMerged)[names(dfMerged) %in% newcols ])

  ## remove indiviual HES data
  VctRemoveCols<-!grepl( "_DO_FUd" , names( dfMerged ) ) & !grepl( "_Read_" , names( dfMerged ) ) & !grepl( "_DOp_FUd" , names( dfMerged ) )  & !grepl( "_Dhes_" , names( dfMerged ) ) & !grepl( "_D9hes_" , names( dfMerged ) ) & !grepl( "_Ohes_" , names( dfMerged ) )
  dfMerged<-dfMerged[ , VctRemoveCols  ] ## remove individual followup data.

  names(dfMerged) %in% dfVctDescriptions[,1]

  dfMerged[!names(dfMerged) %in% "n_eid", ]

  attr(dfMerged , "var.labels") <-  c("identifier",rep(paste(StrDescription," - HES Union:",HEScodes,sep=""),ncol(dfMerged)-1))
  attr(dfMerged, "var.labels") <- substr(attr(dfMerged, "var.labels"),0,80)


  #
  # #### IF =1
  # if (length(which(grepl( "_FUn" , names( dfMerged ) )))>0){
  #   attr(dfMerged, "var.labels") <- c(VctDescriptions[VctRemoveCols[1:length(VctDescriptions)]],paste(StrDescription," - ",cnames[1:4]," - HES Union:",HEScodes,sep=""))
  #   attr(dfMerged, "var.labels") <- substr(attr(dfMerged, "var.labels"),0,80)
  # } else {
  #   attr(dfMerged, "var.labels") <- c(VctDescriptions[VctRemoveCols[1:length(VctDescriptions)]])
  #   attr(dfMerged, "var.labels") <- substr(attr(dfMerged, "var.labels"),0,80)
  # }

  return(dfMerged)
}

#' @export
MultiMergeHESIN_UKBV <- function(files,StrTrait,StrDescription,HEScodes,StataOutputFile){
  #  StataOutputFile="/data_work/databases/ukbiobanks/Phenotypes/CAD_definitions_testNewscript/output/CAD2_testNewFU.dta"
  #  files<-c("/data_work/databases/ukbiobanks/Phenotypes/CAD_definitions_testNewscript/output/ts_53_0_0/HESIN/merged/CAD2_merged.dta","/data_work/databases/ukbiobanks/Phenotypes/CAD_definitions_testNewscript/output/ts_53_0_0/UKBVisit/ALL/CAD2_UKBV.dta"  ,"/data_work/databases/ukbiobanks/Phenotypes/CAD_definitions_testNewscript/output/ts_53_0_0/UKBVisit/CAD2_TSAd.dta"   )
  #  StrTrait<-"CAD2"
  #  HEScode<-"III"

  #  files<-c("/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_0_0/UKBVisit/ALL/DM_UKBV.dta","/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_0_0/UKBVisit//DM_TSAd.dta","/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_0_0/HESIN/merged/DM_merged.dta")
  #  files<-c("/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/UKBVisit/ALL/DM_UKBV.dta","/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/HESIN/merged/DM_merged.dta")
  # files<-c( "/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/HESIN/merged/CvaBh_merged.dta","/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/UKBVisit/ALL/CvaBh_UKBV.dta" )
  #StrDescription<-"asd"
  #StrTrait<-"CvaBh"
  #files<-c("/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/HESIN/merged/Ht_merged.dta", "/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/UKBVisit/ALL/Ht_UKBV.dta"  ,"/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/UKBVisit//Ht_TSAd.dta"     )
  #  StrTrait<-"DM"
  #  files<- "/data_work/databases/ukbiobanks/Phenotypes/test/ts_53_2_0/UKBVisit/ALL/htmeds_UKBV.dta"
  #  StrTrait<-"HF"


  LstDf = lapply(files, function(x){    as.data.frame(read.dta13(x,convert.dates = TRUE))   })
  VctDescriptions<-unlist(lapply(LstDf,function(x){ varlabel(x) } ))
  dfVctDescriptions <- unique(cbind(names(VctDescriptions),unname(VctDescriptions)))

  dfMerged<-Reduce(function(x,y) {merge(x,y,by="n_eid",all=TRUE)}, LstDf)

  StrTrait_HXd<-paste(StrTrait,"_EP_HXd",sep="")
  StrTrait_HXn<-paste(StrTrait,"_EP_HXn",sep="")
  StrTrait_HXt<-paste(StrTrait,"_EP_HXt",sep="")
  StrTrait_HXto<-paste(StrTrait,"_EP_HXto",sep="")
  StrTrait_FUd<- paste(StrTrait,"_EP_FUd",sep="")
  StrTrait_FUn<- paste(StrTrait,"_EP_FUn",sep="")
  StrTrait_FUt<- paste(StrTrait,"_EP_FUt",sep="")
  StrTrait_FUto<- paste(StrTrait,"_EP_FUto",sep="")
  StrTrait_UKBV_HXn<-paste(StrTrait,"_UKBV_HXn",sep="")
  StrTrait_UKBV_HXdlv<-paste(StrTrait,"_UKBV_HXdlv",sep="")
  StrTrait_UKBV_HXd<-paste(StrTrait,"_UKBV_HXd",sep="")
  StrTrait_UKBV_FUd<-paste(StrTrait,"_UKBV_FUd",sep="")
  StrTrait_UKBV_FUdlv<-paste(sep="",StrTrait,"_UKBV_FUdlv")
  StrTrait_UKBV_FUn<-paste(sep="",StrTrait,"_UKBV_FUn")
  StrTrait_TSAd<-paste(StrTrait,"_TSAd",sep="")

  VctAllColumns<-c(StrTrait_HXd,StrTrait_HXn,StrTrait_HXt,StrTrait_HXto,StrTrait_FUd,StrTrait_FUn,StrTrait_FUt,StrTrait_FUto,StrTrait_UKBV_HXn,StrTrait_UKBV_HXdlv,StrTrait_UKBV_HXd,StrTrait_UKBV_FUd,StrTrait_UKBV_FUdlv,StrTrait_UKBV_FUn,StrTrait_TSAd)
  VctTmpColumns<-c("TmpTSAdFuture","TmpTSAdHist")

  VctDummys<-c()
  for (StrColumn in c(VctAllColumns,VctTmpColumns)) {
    if (!StrColumn %in% names(dfMerged)) {
      print(paste(StrColumn," is missing; creating dummy"))
      dfMerged[,StrColumn]<-NA
      VctDummys<-c(VctDummys,StrColumn)
    }
  }


  #### Filter TSAd!=. out if there are no other other observations/evidence (e.g. T2D vs T1D have a shared TSAd column).
  ExcludeTSAd_TF<- -which( names(dfMerged) %in% c(StrTrait_TSAd,"n_eid")  )
  dfMerged<-dfMerged[rowSums( is.na( dfMerged[ , ExcludeTSAd_TF] ))!= ncol(dfMerged[ , ExcludeTSAd_TF]), ]

  ##################
  #### HISTORY:
  #### HX
  dfMerged$HX<-NA
  dfMerged[ dfMerged[,StrTrait_HXn]>0 & !is.na(dfMerged[,StrTrait_HXn]),"HX" ] <-1 ### HES history
  dfMerged[ dfMerged[,StrTrait_UKBV_HXn]>0 & !is.na(dfMerged[,StrTrait_UKBV_HXn]),"HX" ] <-1 ### UKBV history
  dfMerged[ dfMerged[,StrTrait_TSAd]<0 & !is.na(dfMerged[,StrTrait_TSAd]),"HX" ] <-1  #WHAT TO DO WHEN TSAd IS FILLED IN AFTER BASELINE BUT HAPPENS BEFORE THAT VISIT/BASELINE so that UKBV_HXn==., but actually TSAd indicates it should be 1?

  dfMerged[which(dfMerged[,StrTrait_HXn] >1),"HX"]<-dfMerged[which(dfMerged[,StrTrait_HXn] >1),StrTrait_HXn] ## replace HX with number of HES events. if >1

  ### HISTORY-days  Algorithm
  #   HXd= RowMins(HES_HXd,TSAd_if_History)
  #   replace with  UKBV_HXdmean if HXd !=.
  #   replace with UKBV_HXdmean if UKBV_HXd<HXd
  #   replace with NA if HXd!=NA & StrTrait_UKBV_HXdlv ==NA (because that would mean someone was scored positive on a questionaire, but the first age-of-diagnosis could not be determined; However there may be a HES, or weirdly TSAd in the future (altough this is very rare.) after the questionaire, but this would be no first age of diagnosis)
  #
  ###
  ### select TSAd if TSAd<-0 @TmpTSAdHist
  TF<- !is.na(dfMerged[,StrTrait_TSAd])& dfMerged[ , StrTrait_TSAd]<=0
  dfMerged[TF,"TmpTSAdHist"] <-dfMerged[  TF,StrTrait_TSAd]
  ### rowMins (HES_history , TSAd_if_history)
  dfForRowmins<-as.matrix(dfMerged[ , c(StrTrait_HXd,"TmpTSAdHist")] )
  #### NOTE TO SELF, WHY THIS LOOP????
  if ( is.null(nrow( dfForRowmins[ rowSums(is.na(dfForRowmins))!=2 , ] )) ) {
    HXd<-rowMins( dfForRowmins , na.rm=T)
  } else if (   nrow( dfForRowmins[ rowSums(is.na(dfForRowmins))!=2 , ] )==0 ) {  #### CHECK IF ROWSUMS is not empty (nrow ==0) {
    HXd<-as.vector(dfForRowmins[,1])
  } else{
    HXd<-rowMins( dfForRowmins , na.rm=T)
  }
  HXd[HXd==Inf]<-NA
  dfMerged$HXd<-HXd

  ### Replace HXd with days until after the first occurence in history and before the next last visit for which there is no other observation.  0 --<here>--- 1 ------- 1
  NoICD10_TSAdButAvailableVisitdates_TF<-is.na(dfMerged$HXd) & !is.na(dfMerged[,StrTrait_UKBV_HXd]) & !is.na(dfMerged[,StrTrait_UKBV_HXdlv])
  dfMerged[NoICD10_TSAdButAvailableVisitdates_TF,]$HXd<-dfMerged[NoICD10_TSAdButAvailableVisitdates_TF,StrTrait_UKBV_HXd] + (dfMerged[NoICD10_TSAdButAvailableVisitdates_TF,StrTrait_UKBV_HXdlv]  - dfMerged[NoICD10_TSAdButAvailableVisitdates_TF,StrTrait_UKBV_HXd])/2
  ### if currently defined history is after visit that someone said yes, replace with the between date of the visits: ## 0 --<here>--- 1 ----event--- 1
  CheckVisitDates_TF<-!is.na(dfMerged[,StrTrait_UKBV_HXdlv]) & !is.na(dfMerged[,StrTrait_UKBV_HXd]) &!is.na(dfMerged$HXd) & dfMerged[,StrTrait_UKBV_HXd] < dfMerged$HXd
  dfMerged[CheckVisitDates_TF,"HXd"]<- dfMerged[CheckVisitDates_TF,StrTrait_UKBV_HXd] + (dfMerged[CheckVisitDates_TF,StrTrait_UKBV_HXdlv]  - dfMerged[CheckVisitDates_TF,StrTrait_UKBV_HXd])/2
  ### if currently defined history is after visit that someone said yes but there is no visit before, replace with NA since you don't have any evidence of age of diagnosis, altough you have an HES event in the history after someone answerd 'yes' on selfreported.  1 ----- 1 ------- 1
  CheckVisitDates_TF<-is.na(dfMerged[,StrTrait_UKBV_HXdlv]) & !is.na(dfMerged[,StrTrait_UKBV_HXd]) &!is.na(dfMerged$HXd) & dfMerged[,StrTrait_UKBV_HXd] < dfMerged$HXd
  dfMerged[CheckVisitDates_TF,"HXd"]<-NA

  ##################
  ### FOLLOW UP-days Algorithm
  #   FUd=Min( HES_FUd TSAd_future(if there is no history & TSAd is in the future) )
  #   replace with  UKBV_FUdmean if FUd !=. & Hx !=.
  #   replace with UKBV_FUdmean if UKBV_FUd<FUd
  #
  ### Select TSAd in future when no history is present.
  TF<-!is.na(dfMerged[,StrTrait_TSAd]) & dfMerged[ , StrTrait_TSAd]>0 & is.na(dfMerged[,"HX" ])
  dfMerged[ TF,"TmpTSAdFuture"] <-dfMerged[  TF,StrTrait_TSAd]
  ### RowMins( HES_FUd TSAd_future(if there is no history & TSAd is in the future) )
  dfForRowmins<-as.matrix(dfMerged[ , c(StrTrait_FUd,"TmpTSAdFuture")])

  if ( is.null(nrow( dfForRowmins[ rowSums(is.na(dfForRowmins))!=2 , ] )) ) { ### als er maar 1 row is geeft nrow NULL terug, maar doet rowMins() het nog wel...
    FUd<-rowMins( dfForRowmins  , na.rm=T)
  } else if ( nrow( dfForRowmins[ rowSums(is.na(dfForRowmins))!=2 , ] )==0  ) {  #### CHECK IF ROWSUMS is not empty (nrow ==0), dan doet rowMins() het niet..
    FUd<-as.vector(dfForRowmins[,1])
  } else{
    FUd<-rowMins( dfForRowmins  , na.rm=T)
  }
  FUd[FUd==Inf]<-NA
  dfMerged$FUd<-FUd
  ### Include the middle between first occurence at visit during self-reported and last visit without
  ### reported diagnosis for cases that have no FUd by already (by HESIN or by TSAd) (these individuals could have been
  ### ignored above or contain TSAd outliers that were  filtered out) but also no history.
  NoHistory_TF<- (dfMerged[,StrTrait_HXn]==0 | is.na(dfMerged[,StrTrait_HXn])) &
    (dfMerged[,StrTrait_UKBV_HXn]==0 | is.na(dfMerged[,StrTrait_UKBV_HXn]))
  NoFuturedate_TF<-is.na(dfMerged$FUd)
  NoTSAd_TF<-is.na(dfMerged[,StrTrait_TSAd])
  NNA_UKBV_FUd<-!is.na(dfMerged[,StrTrait_UKBV_FUd]) & !is.na(dfMerged[,StrTrait_UKBV_FUd])

  # 1440079 ; 1556857
  #asd<-dfMerged[NoHistory_TF & NoFuturedate_TF & NoTSAd_TF & NNA_UKBV_FUd,]
  dfMerged[NoHistory_TF & NoFuturedate_TF & NoTSAd_TF & NNA_UKBV_FUd,]$FUd<-
    dfMerged[NoHistory_TF & NoFuturedate_TF & NoTSAd_TF & NNA_UKBV_FUd, StrTrait_UKBV_FUdlv] +
    ( dfMerged[NoHistory_TF & NoFuturedate_TF & NoTSAd_TF & NNA_UKBV_FUd, StrTrait_UKBV_FUd]  -
        dfMerged[NoHistory_TF & NoFuturedate_TF & NoTSAd_TF & NNA_UKBV_FUd, StrTrait_UKBV_FUdlv] )/2


  ###OLD if currently defined followup is after visit that someone said yes, replace with the between date of the visits:
  #CheckVisitDates_TF<-!is.na(dfMerged[,StrTrait_UKBV_FUdlv]) & !is.na(dfMerged[,StrTrait_UKBV_FUd]) &!is.na(dfMerged$FUd) & dfMerged[,StrTrait_UKBV_FUd] < dfMerged$FUd
  #dfMerged[CheckVisitDates_TF,"FUd"]<- dfMerged[CheckVisitDates_TF,StrTrait_UKBV_FUd] + (dfMerged[CheckVisitDates_TF,StrTrait_UKBV_FUdlv]  - dfMerged[CheckVisitDates_TF,StrTrait_UKBV_FUd])/2

  ### if currently defined followup is after visit that someone said yes, replace with the between date of the visits; BUT ONLY WHEN THERE IS NO HISTORY:
  CheckVisitDates_TF<-!is.na(dfMerged[,StrTrait_UKBV_FUdlv]) & !is.na(dfMerged[,StrTrait_UKBV_FUd]) &!is.na(dfMerged$FUd) & dfMerged[,StrTrait_UKBV_FUd] < dfMerged$FUd
  dfMerged[CheckVisitDates_TF & NoHistory_TF,"FUd"]<- dfMerged[CheckVisitDates_TF& NoHistory_TF,StrTrait_UKBV_FUd] + (dfMerged[CheckVisitDates_TF& NoHistory_TF,StrTrait_UKBV_FUdlv]  - dfMerged[CheckVisitDates_TF& NoHistory_TF,StrTrait_UKBV_FUd])/2


  #### FU:
  dfMerged$FU<-NA
  dfMerged[which(dfMerged$FUd>0) ,"FU"]<-1 #which(dfMerged$FUd>0) # !is.na(dfMerged$FUd) & dfMerged$FUd>0

  UKBV_FU_UKBVisitOnly_TF<-which ( dfMerged[,StrTrait_UKBV_FUn]>0 & is.na(dfMerged$HX) & is.na(dfMerged$FU) ) ### some individuals have no date-column, so they dont have a days;
  dfMerged[ UKBV_FU_UKBVisitOnly_TF ,"FU"]<-1
  dfMerged[which(dfMerged[,StrTrait_FUn] >1),"FU"]<-dfMerged[which(dfMerged[,StrTrait_FUn] >1),StrTrait_FUn] ## replace FU with number of HES events. if >1

  #### REMOVE THE DUMMY VARIABLES
  dfMerged<-dfMerged[,!names(dfMerged) %in% c(VctDummys), drop=F]

  #### ADD ANY VARIABLE
  dfMerged$ANY<-1
  ### check which sources were used:

  #### FIX NAMES:
  # names(dfMerged)[length(names(dfMerged)):(length(names(dfMerged))-4)]<-paste(StrTrait,"_",names(dfMerged)[length(names(dfMerged)):(length(names(dfMerged))-4)] , sep="")
  newcols <- c("HXd","FUd","FUn","FU","ANY","HX","HXn")
  names(dfMerged)[names(dfMerged) %in% newcols ] <- paste(StrTrait,"_", names(dfMerged)[names(dfMerged) %in% newcols] , sep="")

  dfVctDescriptions
  dfVctDescriptions[ dfVctDescriptions[,1] %in% names(dfMerged) ,]

  attr(dfMerged, "var.labels") <- c("identifier",rep(paste(StrDescription," - Combined"),ncol(dfMerged)-1))

  #
  # attr(dfMerged, "var.labels") <-
  #   c(VctDescriptions,
  #     paste(StrDescription," - ",names(dfMerged)[  (length(names(dfMerged))-4):(length(names(dfMerged)) )]," - Combined",sep="")
  #   )
  # attr(dfMerged, "var.labels")<-substr(attr(dfMerged, "var.labels"),0,80)

  print(names(dfMerged))
  save.dta13(dfMerged, StataOutputFile,compress = TRUE)
  return(dfMerged)
}



#' @export
MultiMergeHESIN_medication_UKBV <-function(files,StrTrait,StrDescription,HEScodes,StataOutputFile){


}





