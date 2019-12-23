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

