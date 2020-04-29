#' @export
merge_stata_files <- function(statafiles=list.files(dir.ukbpheno_all,full.names = T ),f.out,n_eids=NULL){
  df.all <- data.frame(n_eid=0)
  attr(df.all,"var.labels") <- "Identifier"
  df.all_attr = data.frame(colname="n_eid",attr="Identifier")



  for (f in statafiles ){ 
    print(f) 
    trait=strsplit(basename(f),split="_")[[1]][1]
    columns_to_keep = c("n_eid", paste0(trait,"_DO_FUn"),paste0(trait,"_DOp_FUn"), paste0(trait,"_HX"), paste0(trait,"_HXd"), paste0(trait,"_FU"), paste0(trait,"_FUd"), paste0(trait,"_ANY"))
    
    d <- read.dta13(f) 
    d_attr <- as.data.frame(cbind(names(d),attr(d,"var.labels")))
    names(d_attr) <- c("colname","attr")

    d_attr <- d_attr[names(d) %in% columns_to_keep,]
    d <- d[,names(d) %in% columns_to_keep]
    
    df.all_attr <- unique(rbind(d_attr,df.all_attr))
    df.all = merge(df.all,d,by="n_eid",all=TRUE)

    gc()
  }
  Vctattrs <- df.all_attr[match(names(df.all), df.all_attr$colname),"attr"]
  attr(df.all,"var.labels") <- as.character(Vctattrs)
  print(paste("writing:",f.out))
  #write_dta(df.all, f.out, version = 14)
  save.dta13(df.all, f.out,compress = FALSE)
  return(df.all)
}

#' @export
summarize_cross_phenotype_fields <- function(dir.ukbpheno,lst.n_eids){
  searchsources=c("DO","DOp","Dhes","D9hes","O4hes","O3hes","TS","SR","TS_RX","SR_RX","Read","CTV3","BNF","DMD",
                  "UKBV", "HESIN","GP","ALL")

  files.all <- list.files(dir.ukbpheno,full.names = T,recursive = T )
  dfpheno <- as.data.frame(cbind(files.all,unlist(lapply(files.all,basename))))
  dfpheno <- cbind(dfpheno,do.call(rbind,stringi::stri_split_fixed(str = dfpheno$V2, pattern = "_", n = 2)))
  dfpheno$`2` <-  sub(".dta$","",dfpheno$`2`)
  names(dfpheno)<- c("path","file","pheno","source")
  dfpheno <- dfpheno[dfpheno$source %in% searchsources,]
  dfpheno <- dfpheno[!dfpheno$file %in% "ALL.dta",]
  dfcrosscounts <- data.frame()
  for (p in unique(dfpheno$pheno)) {
    print(p)
    dfphenocomparison <- data.frame(n_eid=-9)
    sources <- dfpheno[(dfpheno$pheno %in% p),"source"]
    for (s in sources) {
      path = as.character(dfpheno[(dfpheno$pheno %in% p) & dfpheno$source == s,"path"])
      d <- as.data.frame(cbind(read.dta13(path)$n_eid,1))
      names(d) <- c("n_eid",s)
      dfphenocomparison = merge(dfphenocomparison,d,by = "n_eid",all = T)
    }
    dfphenocomparison = dfphenocomparison[dfphenocomparison$n_eid != -9,]
    dfphenocomparison[is.na(dfphenocomparison)] <- 0
    dfphenocomparison = dfphenocomparison[dfphenocomparison$n_eid %in% lst.n_eids,] # restrict to n_eids after exclusions.

    summarize_matrix <- function(dfphenocomparison,sources){
      dfcounttable <- data.frame()
      for (s in sources){
        for (ss in sources){
          d <- melt(table(dfphenocomparison[,s] ,dfphenocomparison[,ss]) )
          d$s <-s
          d$ss <- ss
          d$pheno <- p
          names(d) <- c("source1_factor","source2_factor","count","source1_name","source2_name","pheno")

          dfcounttable=rbind(dfcounttable,d)
        }
      }
      return(dfcounttable)
    }
    if(nrow(dfphenocomparison)>0){
      dfcrosscounts <- rbind(dfcrosscounts,summarize_matrix(dfphenocomparison,sources))
    }
  }
  # reorder
  dfcrosscounts= dfcrosscounts[,c("pheno","source1_name","source2_name","source1_factor","source2_factor","count")]
  return(dfcrosscounts)
}



# summarize_cases_for_codes(codes="H3,X101n,XaEIV,XaEIW,XaEIY,XaIND,XaK8Q,XaN4a,Xac33,H3y0,Xa35l,H3122,X101i,XaZd1,H0611,H31,H3120,H3121,H312z,Hyu31,X101j,X101l,X102z,XSDOK,XaDtP,H3,H3y,H3z,XaK8Q,Xaa7B,H3y,14B3,X101i,Xa35l,XaZd1,XE2Pp",
#                           df=dfgpclinical,
#                           diagcols=c("read_3") ){
# VctCodes = strsplit(codes,split=",")[[1]]
# if (length(diagcols)==1){diagcols = c(diagcols,diagcols)}
# grepoper <-unlist( mclapply(  df[, diagcols], function(c) grep(paste(sep="","^",VctCodes, collapse='|'), c, ignore.case=FALSE),mc.cores =detectCores()/2 ) ) ## PARALLEL of the above.
#
# lstsummary <- list(all= length(unique(df[grepoper,]$n_eid ) ))
#
# for (c in VctCodes){
#   print(c)
#   grepoper <-unlist( mclapply(  df[, diagcols], function(c) grep(paste(sep="","^",c, collapse='|'), c, ignore.case=FALSE),mc.cores =detectCores()/2 ) ) ## PARALLEL of the above.
#   lstsummary[[c]]<- length(unique(df[grepoper,]$n_eid ))
# }
#   print(lstsummary)
# }
#


