#' @export
merge_stata_files <- function(statafiles=list.files(dir.ukbpheno_all,full.names = T ),f.out){

  dfs <- list()
  for (f in statafiles ){  dfs[[f]] = read.dta13(f) }

  df.all <- data.frame(n_eid=0)
  attr(df.all,"var.labels") <- "Identifier"
  df.all_attr = data.frame(colname="n_eid",attr="Identifier")
  for (d in dfs){
    d_attr <- as.data.frame(cbind(names(d),attr(d,"var.labels")))
    names(d_attr) <- c("colname","attr")
    df.all_attr <- unique(rbind(d_attr,df.all_attr))
    df.all = merge(df.all,d,by="n_eid",all=TRUE)
  }
  Vctattrs <- df.all_attr[match(names(df.all), df.all_attr$colname),"attr"]
  attr(df.all,"var.labels") <- as.character(Vctattrs)
  print(paste("writing:",f.out))
  #write_dta(df.all, f.out, version = 14)
  save.dta13(df.all, f.out,compress = TRUE)
  return(df.all)
}


#' @export
summarize_cross_phenotype_fields <- function(dir.ukbpheno,lst.n_eids){
  searchsources=c("DO","DOp","Dhes","D9hes","Ohes","TS","SR","TS_RX","SR_RX","ALL","Read","CTV3","BNF","DMD")

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
    dfcrosscounts <- rbind(dfcrosscounts,summarize_matrix(dfphenocomparison,sources))

  }
  # reorder
  dfcrosscounts= dfcrosscounts[,c("pheno","source1_name","source2_name","source1_factor","source2_factor","count")]
  return(dfcrosscounts)
}
