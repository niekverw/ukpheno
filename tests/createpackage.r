library(data.table)
library("devtools")
library(roxygen2)

setwd("/data_work/bitbucket/ukbio-github/ukpheno")
document()
rm(dfDefinitions)
install("/data_work/bitbucket/ukbio-github/ukpheno")
library(CreateUKBiobankPhentoypes)
####
#remove.packages("CreateUKBiobankPhentoypes")

### save example dataset:
dfDefinitions_file_tsv="/data_work/bitbucket/ukbio-github/ukpheno/data/dfDefinitions.tsv"
dfDefinitions<-data.frame(fread(dfDefinitions_file_tsv))
dfDefinitions$Comment<-""
save(dfDefinitions, file="data/dfDefinitions.RData")
### save codesheet for treatment:
#....dfCodesheetTreatment<-data.frame(fread())
save(dfCodesheetTreatment, file="data/dfCodesheetTreatment.RData")
##### 

build()

install.packages( "/data_work/bitbucket/ukbio/CreateUKBiobankPhentoypes_0.21.tar.gz",type="source",repos=NULL)
install.packages( "CreateUKBiobankPhentoypes_0.21.tar.gz",type="source",repos=NULL)

library(CreateUKBiobankPhentoypes)

####

remove.packages("CreateUKBiobankPhentoypes")
unload(CreateUKBiobankPhentoypes)
CreateUKBiobankPhentoypes(3,2 ,UKbioDataset,dfmaster_SQL_merge=dfmaster_SQL_merge,dfSQL_TraitTable,Outputdir )

data(dfDefinitions)
data(dfCodesheetTreatment)
