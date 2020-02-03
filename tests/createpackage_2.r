library(data.table)
library("devtools")
library(roxygen2)
#devtools::install_github("klutometis/roxygen")

d="/data_work/bitbucket/ukbio/CreateUKBiobankPhentoypes"
d="~/tmpukb/ukpheno"
setwd(d)
roxygen2::roxygenise() # just
devtools::document() #srouce + compile

detach(name= "CreateUKBiobankPhentoypes", unload=TRUE)
rm(dfDefinitions)
rm(CreateUKBiobankPhentoypes)

devtools::install(d)
library(CreateUKBiobankPhentoypes)
####
#remove.packages("CreateUKBiobankPhentoypes")
?CreateUKBiobankPhentoypes



### save example dataset:
dfDefinitions_file_tsv="data/dfDefinitions.tsv"
dfDefinitions<-data.frame(fread(dfDefinitions_file_tsv))
dfDefinitions$Comment<-""
save(dfDefinitions, file="data/dfDefinitions.RData")
### save codesheet for treatment:
#....dfCodesheetREAD_SR.Coding<-data.frame(fread())

# coding sheets:
## READCODES V2 TO UKB SELFREPORTED
#save(dfCodesheetREAD_SR.Coding, file="data/dfCodesheetREAD_SR.Coding.RData")




build()

install.packages( "/data_work/bitbucket/ukbio/CreateUKBiobankPhentoypes_0.21.tar.gz",type="source",repos=NULL)
install.packages( "CreateUKBiobankPhentoypes_0.24.tar.gz",type="source",repos=NULL)

library(CreateUKBiobankPhentoypes)

####

remove.packages("CreateUKBiobankPhentoypes")
unload(CreateUKBiobankPhentoypes)
CreateUKBiobankPhentoypes(3,2 ,UKbioDataset,dfmaster_SQL_merge=dfmaster_SQL_merge,dfSQL_TraitTable,Outputdir )

data(dfDefinitions)
data(dfCodesheetREAD_UKB.Coding)
