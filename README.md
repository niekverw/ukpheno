# ukpheno
Creates dichotomous phenotypes for UKbio using all available data (ICD/oper/medication/self reports/age of diagnosis/visit-dates etc). The current output includes variables on history, study visit, future, time-to-first-event. happy to help. Mail@niekverweij.com.

Used in: 
- JAMA Cardiol. 2018 Aug 1;3(8):693-702. doi: 10.1001/jamacardio.2018.1717.
- Circ Res. 2018 Feb 2;122(3):433-443. (https://www.ncbi.nlm.nih.gov/pubmed/29212778)
- J Am Heart Assoc. 2018 Apr 5;7(8). pii: e008341. doi: 10.1161/JAHA.117.008341.
- Sci Rep. 2018 Apr 11;8(1):5817. doi: 10.1038/s41598-018-24002-0.
- Nat Commun. 2018 Mar 1;9(1):898. doi: 10.1038/s41467-018-03395-6.
- Nat Commun. 2018 Mar 7;9(1):987. doi: 10.1038/s41467-018-03252-6.
- J Am Heart Assoc. 2018 Jan 22;7(2). pii: e007621. doi: 10.1161/JAHA.117.007621.
- J Am Coll Cardiol. 2017 Jul 25;70(4):506-507. doi: 10.1016/j.jacc.2017.05.044
- Sci Rep. 2017 Jun 5;7(1):2761. doi: 10.1038/s41598-017-03062-8. <- please cite if software is used. 



## Input: 

- 1 datatable of UK Bioank in standard STATA format 
- 4 tables containing hospital records of 1) Primary ICD10 diagnoses, 2) Secondary ICD10, 3) ICD9 and 4) OPERATION CODES. 
- 1 table containing the diagnoses: example https://github.com/niekverw/ukpheno/blob/master/data/dfDefinitions.tsv

## example:
```
library(CreateUKBiobankPhentoypes)
library(readstata13)
library(data.table)

UKbioDataset_file = "/path/to/file.dta"
hesin_file="/path/to/hesin_2018-04-17.tsv"
hesin_diagicd10_file="/path/to/hesin_diagicd10_2018-04-17.tsv"
hesin_diagicd9_file="/path/to/hesin_diagicd9_2018-04-17.tsv"
hesin_oper_file="/path/to/hesin_oper4_2018-04-17.tsv"
dfDefinitions_file = "/path/to/https://github.com/niekverw/ukpheno/blob/master/data/dfDefinitions.tsv"
Outputdir="/path/to/output"

print("load definition table")
dfDefinitions = data.frame(fread(dfDefinitions_file))
write.table(ProcessDfDefinitions(dfDefinitions),paste(dfDefinitions_file,".check.tsv",sep=""),sep="\t",quote=FALSE,row.names = FALSE) # used to debug your definitions.
print("load dataframe ukbiobank")
UKbioDataset <- as.data.frame(read.dta13(UKbioDataset_file,convert.dates = TRUE))
print("load hesin")
dfhesintables<-LoadHesinTable(UKbioDataset,hesin_file,hesin_diagicd10_file,hesin_diagicd9_file,hesin_oper_file)

print("constructing diagnoses for baseline visit. 
CreateUKBiobankPhentoypes(Nvisits=3,
                          visitreference=0,
                          UKbioDataset,
                          dfhesintables,
                          dfDefinitions,
                          Outputdir,
                          VctOutputIndividualColumns=c("TS","SR","TS_RX","RX","LAB")
                          )
```
