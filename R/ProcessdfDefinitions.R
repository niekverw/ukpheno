ConvertFactorsToStringReplaceNAInDf<-function (df) {
  df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE) ## CHANGE Factors to strings; everything is now a string.
  df[df==""]  <- NA ### REPLACE EMPTY WITH NA.
  return(df)
}

pasteRemoveNA <- function(..., sep = " ", collapse = NULL, na.rm = F) {
  if (na.rm == F)
    paste(..., sep = sep, collapse = collapse)
  else
    if (na.rm == T) {
      paste.na <- function(x, sep) {
        x <- gsub("^\\s+|\\s+$", "", x)
        ret <- paste(na.omit(x), collapse = sep)
        is.na(ret) <- ret == ""
        return(ret)
      }
      df <- data.frame(..., stringsAsFactors = F)
      ret <- apply(df, 1, FUN = function(x) paste.na(x, sep))

      if (is.null(collapse))
        ret
      else {
        paste.na(ret, sep = collapse)
      }
    }
}

CheckDuplicateTRAITS<-function(df){
  if(length(unique(duplicated(df["TRAIT"])))>1){stop("TRAIT column contains duplicate ID's")}
}




PreProcessDfDefinitions<-function(df,VctAllColumns){
## df<-dfDefinitions
  ## for the names: remove everything between dots (R converts symbols to dots "(,.-)/" etc )
  names(df)<-gsub( " *\\..*?\\. *", "", names(df) )
  ## remove everything between brackets
  df[,VctAllColumns]<- data.frame(apply(df[,VctAllColumns], 2, function(y) gsub( " *\\(.*?\\) *", "", y)) )

  ### remove dots(.): names(df)
  df[,VctAllColumns]<- data.frame(apply(df[,VctAllColumns],2,function(x) gsub(".", "", x, fixed = TRUE)))
  ### remove spaces:
  df[,VctAllColumns]<- data.frame(apply(df[,VctAllColumns],2,function(x) gsub(" ", "", x, fixed = TRUE)))
  ### remove trailing commas:
  trim.commas <- function (x) gsub("(?<=[\\,])\\,*|^\\,+|\\,+$", "", x, perl=TRUE)
  df[,VctAllColumns]<- data.frame(apply(df[,VctAllColumns],2,function(x) trim.commas(x)))

  df<-ConvertFactorsToStringReplaceNAInDf(df) #### CONVERT FACTOR TO STRING

  return(df)
}

FillInSRdefinitions<-function(df,Var="SR",cols=c("n_20001_","n_20002_","n_20004_") ) {
  ## fill in SR
  df[,Var]<-as.character(df[,Var])
  df[is.na( as.character(df[,Var])) ,][,Var] <- ""

 # df [ is.na(as.character( df[,Var] )) %in% "NA"  ,]

  cols<- cols[cols %in% colnames(df)] ### check if available.
  for(col in cols){
    #print(col)
    Columnmatches<-gsub ( ",","|",paste(col,unlist(df[col]),sep="=" ))
    Columnmatches [grepl("NA",Columnmatches)]<-"" ### remove NAs..
    df[,Var]<-  paste(df[,Var],Columnmatches ,sep="," )
  }

  ### trim commas:
  trim.commas <- function (x) gsub("(?<=[\\,])\\,*|^\\,+|\\,+$", "", x, perl=TRUE)
  df[,Var]<-trim.commas(df[,Var])


  df<-ConvertFactorsToStringReplaceNAInDf(df)

  return(df)
}

CovertMednamesToUkbcoding<- function(StrRx){
  #StrRx<-"phenformin,metformin,buformin,glibenclamide,chlorpropamide,tolbutamide,glibornuride,tolazamide,carbutamide,glipizide,gliquidone,gliclazide,metahexamide,glisoxepide,glimepiride,acetohexamide,glymidine,acarbose,miglitol,voglibose,troglitazone,rosiglitazone,pioglitazone,sitagliptin,vildagliptin,saxagliptin,alogliptin,linagliptin,gemigliptin,repaglinide,nateglinide,exenatide,pramlintide,benfluorex,liraglutide,mitiglinide,dapagliflozin,lixisenatide,canagliflozin,empagliflozin,albiglutide,dulaglutide"
  StrRx<-as.character(StrRx)
  if(is.na(StrRx)) { return(NA)}
  VctRXstrings<-unlist(strsplit(StrRx,","))
  #VctRXstrings<-strsplit(df[!is.na(df$n_20003_),]$n_20003_,",")[[1]]
  StrRxCodes<-paste(unique(unlist(lapply(VctRXstrings,  function(x) dfCodesheetTreatment[,"UKB.Coding"] [ grep(x,dfCodesheetTreatment[,"Meaning"] ,ignore.case=TRUE )]  ))),collapse=",")
  return(StrRxCodes)
}

CovertReadcodesToUkbcoding<- function(StrRx){
 # StrRx<-"f3,f4,ft"
  StrRx<-as.character(StrRx)
  if(is.na(StrRx)) { return(NA)}
  VctRXstrings<-unlist(strsplit(StrRx,","))
  #VctRXstrings<-strsplit(df[!is.na(df$n_20003_),]$n_20003_,",")[[1]]
  StrRxCodes<-paste(unique(unlist(lapply(VctRXstrings,  function(x) dfCodesheetTreatment[,"UKB.Coding"] [ grep(paste("^", x,sep=""),dfCodesheetTreatment[,"READ.CODE"] ,ignore.case=TRUE )]  ))),collapse=",")
  return(StrRxCodes)
}


#### TODO:
ReduceRedundancyDf<- function(df){ ### NOT really nessesary
  return(df)
}



# DfDefinitions<-read.table("/Users/niekverw/Downloads/ex",sep="\t",header=T)
#columns<-c("ICD10CODES","ICD9CODES","OPERCODES","TOUCHSCREEN","TS_AGE_DIAG_COLNAME","SELFREPORTED","MEDICATION","LAB")
# print(dfDefinitions)
#dfDefinitionstmp2<-ProcessDfDefinitions(dfDefinitions,columns)

#' ProcessDfDefinitions
#'
#' Process definitions, for input
#'
#' @param dfDefinitions df
#' @param VctAllColumns Vct
#' @keywords ExtractVarsFromMasterSet CreateUKBiobankPhentoypes ProcessDfDefinitions
#' @return None
#'
#' @examples
#' #
#' #This function processes an excel file with definitions and is automtically performed in CreateUKBiobankPhentoypes().
#' #It can be usefull to run this function as a check prior to running CreateUKBiobankPhentoypes.
#' #
#' #VctAllColumns contains all column names of interest, so that it can ignore everything else.
#' #20001, 20002 and 20004 go into SR
#' #READCODES and 20003 is parsed into RX
#' #
#' #
#' #
#' VctAllColumns<-  c("TS", "SR", "TS_RX", "RX", "LAB", "ICD10CODES", "ICD9CODES", "OPERCODES", "TS_AGE_DIAG_COLNAME", "READCODES", "n_20001_",    "n_20002_", "n_20003_", "n_20004_", "DEPENDENCY")
#' ProcessDfDefinitions(dfDefinitions,VctAllColumns)
#'
#' @export
ProcessDfDefinitions<-function(df,VctAllColumns=c("TS", "SR", "TS_RX", "RX", "LAB", "ICD10CODES", "ICD9CODES", "OPERCODES", "TS_AGE_DIAG_COLNAME", "READCODES", "n_20001_",    "n_20002_", "n_20003_", "n_20004_", "DEPENDENCY")){
 # wb <- loadWorkbook(dfDefinitions_file)
 # dfDefinitions =readWorksheet(wb, sheet = 1)
 # df<- dfDefinitions  #  df<- dfDefinitions2
  # VctAllColumns<-  c("TS", "SR", "TS_RX", "RX", "LAB", "ICD10CODES", "ICD9CODES", "OPERCODES", "TS_AGE_DIAG_COLNAME", "READCODES", "n_20001_",    "n_20002_", "n_20003_", "n_20004_", "DEPENDENCY")

  if(nrow(df)==1 ) {stop("please have more than 1 phenotype definition.")} ## check if excel file has more than 1 row.

  df<-PreProcessDfDefinitions(df,VctAllColumns)
  df<-FillInSRdefinitions(df,"SR",c("n_20001_","n_20002_","n_20004_"))

  ### LOOKUP NAMES OF MEDICATION and put UKBIO.CODES in RX
  df$n_20003_<-unlist(lapply( df$n_20003_, CovertMednamesToUkbcoding))
  df<-FillInSRdefinitions(df,"RX",c("n_20003_"))

  ### LOOKUP READ.CODES and put UKBIO.CODES in RX
  df$n_20003_<-unlist(lapply( df$READCODES, CovertReadcodesToUkbcoding))
  df<-FillInSRdefinitions(df,"RX",c("n_20003_"))

  CheckDuplicateTRAITS(df) # check duplicateids.
  ### df = excel matrix. 1 hij loopt een voor een over elke rij heen,
  # 2 zoekt per dependency in die rij de bijpassende rijen voor elke dependency  erbij en plakt die naast elkaar (inclusief de dependencies van de dependencies).
  # 4) dan delete hij de depenencies die hij ingevuld heeft # dit op repeat tot dat er geen dependencies meer zijn en alles is ingevuld.
  repeat {
    for(i in 1:nrow(df)) {
      row <- df[i,]

      if(!is.na(row$DEPENDENCY)){
        VctDEPENDENCYs<-unlist(strsplit(row$DEPENDENCY,","))

        for (StrDEPENDENCY in VctDEPENDENCYs) {
          if(row$TRAIT == StrDEPENDENCY) {stop("Dependency is same as trait.")}
          targetrow<-df[df$TRAIT==StrDEPENDENCY,]
          if(nrow(targetrow)==0){stop(paste('Dependency: "',StrDEPENDENCY,'" not found in traits ',row$TRAIT,sep="")) }

          if( is.na(targetrow["DEPENDENCY"])){
            for(col in VctAllColumns){
              Vctcol<-unique( unlist(strsplit( c(df[i,col],df[df$TRAIT==StrDEPENDENCY,col]) ,",")) )

              df[i,col]<-pasteRemoveNA(Vctcol ,collapse=",",na.rm=T)
            }
            # remove DEPENDENCY that was just filled in:
            #df[i,"DEPENDENCY"]<-gsub(paste(StrDEPENDENCY,sep=""),"",df[i,"DEPENDENCY"],fixed=TRUE,ignore.case=FALSE)

            # remove DEPENDENCY that was just filled in:
            LstTmpDependencies<- unlist(strsplit(VctDEPENDENCYs,","))
            df[i,"DEPENDENCY"]<-paste( LstTmpDependencies [! LstTmpDependencies  %in%  StrDEPENDENCY] ,sep="",collapse = ",")

            df[i,"DEPENDENCY"]<-gsub("^,*|(?<=,),|,*$", "", df[i,"DEPENDENCY"], perl=T)
            if(df[i,"DEPENDENCY"]==""){df[i,"DEPENDENCY"]<-NA} ## if empty replace with NA
          }
        }
      }
    }
    if( length(unique(is.na(df$DEPENDENCY)))==1 ) break
  }
  df<-ConvertFactorsToStringReplaceNAInDf(df)
  return(df)

  #write.table(df,paste(dfDefinitions_file,".processed.tsv",sep=""),quote = FALSE,row.names = FALSE,sep="\t")
}
