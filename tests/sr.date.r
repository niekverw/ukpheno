library(data.table)






convert_selfreport_to_episodedata <- function(UKbioDataset,field_sr_diagnosis = "20002",field_sr_date = "20009"){
  df = UKbioDataset
  # field_sr_diagnosis = "20002"
  # field_sr_date = "20009"
  # 
  visits = sum(grepl("53_", names(df)))
  
  df <- df[,c("n_eid","n_34_0_0","n_52_0_0", names(df)[grepl("53_", names(df))],
              names(df)[grepl(field_sr_diagnosis, names(df))], 
              names(df)[grepl(field_sr_date, names(df))]
  )]
  daysinyear=365.25
  # add age and convert visit-date field
  for (v in 0:(visits-1)){
    print(v)
    agefield = paste0("age_",v)
    visitdatefield = paste0("ts_53_",v,"_0")
    
    df[,visitdatefield] <- as.Date (df[,visitdatefield],format="%d%b%Y")
    df[,agefield] <- as.vector((df[,visitdatefield]- as.Date(paste0(df[,"n_34_0_0"],"/",df[,"n_52_0_0"],"/14"),format = "%Y/%m/%d")) / daysinyear)
  }
  
  dfout <-  matrix(ncol=3, nrow=0)
  for (v in 0:(visits-1)){
    
    diagfields = names(df)[grepl(paste0("n_",field_sr_diagnosis,"_",v),names(df))]
    diagagefields = names(df)[grepl(paste0("n_",field_sr_date,"_",v),names(df))]
    
    for (i in 1:length(diagfields)){
      agefield = paste0("age_",v)
      visitdatefield = paste0("ts_53_",v,"_0")
      diagfield = diagfields[i]
      diagagefield = diagagefields[i]
      
      print(paste0((diagfield), "- ", diagagefield))
      
      dfsub <- df[!is.na(df[,agefield] ) & !is.na(df[,diagfield] ) ,c("n_eid", visitdatefield,agefield,diagfield,diagagefield)]
      dfsub[,diagagefield][dfsub[,diagagefield] <0] <- NA
      
      dfsub$eventdate = dfsub[,visitdatefield] - (dfsub[,diagagefield]*daysinyear)
      
      dfout <- rbind(dfout,as.matrix(dfsub[,c("n_eid","eventdate",diagfield)]))
    }
    
  }
  
  dfout <- as.data.frame(dfout,stringsAsFactors=F)
  names(dfout) <- c("n_eid","eventdate",paste0(field_sr_diagnosis))
  dfout$eventdate <- as.Date(dfout$eventdate,"%Y-%m-%d")
  dfout[,field_sr_diagnosis] <- as.integer(dfout[,field_sr_diagnosis])
  
  
  # deduplicate, min/max/mean/sd
  library(dplyr)
  dfout_extrastats<- dfout %>% group_by(n_eid,!!as.name("20002")) %>%
    mutate(mindt = min(eventdate, na.rm = TRUE),maxdt = max(eventdate, na.rm = TRUE),meandt = mean(eventdate, na.rm = TRUE))
  
  dfout_extrastats$diffdt <- (dfout_extrastats$maxdt - dfout_extrastats$mindt)/daysinyear
  dfout_extrastats[dfout_extrastats$diffdt>10 ,"meandt"] <- NA
  
  dfout$eventdate <- dfout_extrastats$meandt
  dfout <- dfout[!dfout$`20002` %in% 99999,]
  dfout[,"epiend_1"] <- dfout$eventdate
  
  dfout <- dfout[!duplicated(dfout),]
  
  dfout <- merge(dfout, df[,c("n_eid",paste0("ts_53_",0:2,"_0"))], by.x = "n_eid",  by.y = "n_eid")
  head(dfout)
  return(dfout)
}

f="test.200002_20009.tsv"
df <- data.frame(fread(f))
df_sr20001 <- convert_selfreport_to_episodedata(df,field_sr_diagnosis = "20004",field_sr_date = "20007") # cancer
df_sr20002 <- convert_selfreport_to_episodedata(df,field_sr_diagnosis = "20002",field_sr_date = "20009") # noncancer
# age of diagnosis for medication, 20003, not available. 
df_sr20004 <- convert_selfreport_to_episodedata(df,field_sr_diagnosis = "20004",field_sr_date = "20011") # operation
  
  

# dfmaster_SQL_merge<-dfhesintables
# i=1
# visitreference=0
# VctOutputIndividualColumns=c("TS","SR","TS_RX","SR_RX","LAB")
visitdt=paste("ts_53_",visitreference,"_0",sep="")

row <- dfDefinitions[41,]

StrTrait <- row$TRAIT
epidurfilter <- max(0,row$Minimum_Episode_duration,na.rm = T) # 0 by default in case its not filled out.
include_secondary <- min(1,row$include_secondary,na.rm = T) # 1 by default in case its not filled out.

Strcatagory = "SR"
if( !is.na(row$n_20002_.noncancer.)  ) {
  #### SETTINGS:
  print("    ..finding READCODES diagnosiscodes")
  VctCodes<-unlist(strsplit(row$n_20002_.noncancer.,","))
  ####################################
  ### READ: FOLLOW UP + HISTORY VARIABLES:
  ###
  StrName="20002"
  StrColumnForHescodes<-c("20002")
  StrColumnnameForEventDate<-"eventdate"
  StrDescription<-paste(row$DESCRIPTION,"- BNF -",paste(VctCodes,collapse=","))
  StrColumnnameForVisitDate<-visitdt
  
  #dir.create( paste( Outputdir,"/",visitdt,"/",Strcatagory,sep=""), showWarnings =F,recursive=T)
  #StataOutputFile= paste(Outputdir,"/",visitdt,"/",Strcatagory,"/",paste(StrTrait,StrName,sep="_"),".dta",sep="")
  StataOutputFile="~/test20002.dta"
  # fix epiend_1 is not available in GP, so added it to the dataframe.. but it is not useful here.
  diags20002 <- Outcome_HES(dfmaster_SQL_merge = df_sr20002,
              StrTrait = paste(StrTrait,StrName,sep="_"),
              StrDescription,VctCodes,epidurfilter,
              StrColumnForHescodes,StrColumnnameForEventDate,StrColumnnameForVisitDate,StataOutputFile)
  
}