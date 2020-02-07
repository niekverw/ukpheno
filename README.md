# ukpheno
Creates dichotomous phenotypes for UKbio and a composite time-to-event variable  using ICD/oper/medication/self reports/age of diagnosis/visit-dates etc. 
The current output includes variables on history, study visit, future, time-to-first-event, episode duration. If ICD10/9 is used for follow-up, it's possible to change the baseline date to for example an age of diagosis for studying disease trajectory, for ICD9/10 codes its possible to set a treshold for episode duration. mail@niekverweij.com.

Used in: 
- JAMA Cardiol. 2018 Aug 1;3(8):693-702. doi: 10.1001/jamacardio.2018.1717. Associations of Combined Genetic and Lifestyle Risks With Incident Cardiovascular Disease and Diabetes in the UK Biobank Study.
- J Am Heart Assoc. 2018 Apr 5;7(8). pii: e008341. doi: 10.1161/JAHA.117.008341. Heart Rate Recovery 10 Seconds After Cessation of Exercise Predicts Death.
- Sci Rep. 2018 Apr 11;8(1):5817. doi: 10.1038/s41598-018-24002-0. Causal Pathways from Blood Pressure to Larger Qrs Amplitudes a Mendelian Randomization Study.
- Nat Commun. 2018 Mar 1;9(1):898. doi: 10.1038/s41467-018-03395-6. Genetic study links components of the autonomous nervous system to heart-rate profile during exercise.
- Nat Commun. 2018 Mar 7;9(1):987. doi: 10.1038/s41467-018-03252-6.Genome-wide analysis yields new loci associating with aortic valve stenosis.
- Circ Res. 2018 Feb 2;122(3):433-443. https://www.ncbi.nlm.nih.gov/pubmed/29212778. Identification of 64 Novel Genetic Loci Provides an Expanded View on the Genetic Architecture of Coronary Artery Disease.

- J Am Heart Assoc. 2018 Jan 22;7(2). pii: e007621. doi: 10.1161/JAHA.117.007621. Relationship of Arterial Stiffness Index and Pulse Pressure With Cardiovascular Disease and Mortality.
- J Am Coll Cardiol. 2017 Jul 25;70(4):506-507. doi: 10.1016/j.jacc.2017.05.044. Telomere Length and Risk of Cardiovascular Disease and Cancer.  
- Sci Rep. 2017 Jun 5;7(1):2761. doi: 10.1038/s41598-017-03062-8. Identification of 15 novel risk loci for coronary artery disease and genetic risk of recurrent events, atrial fibrillation and heart failure. **<- please cite if software is used.**



## Input: 

- 1 datatable of UK Bioank in standard STATA format 
- 4 tables containing hospital records of 1) Primary ICD10 diagnoses, 2) Secondary ICD10, 3) ICD9 and 4) OPERATION CODES. 
- 2 tables containing hospital primary care data, 1) clinical, 2)scripts. 
- 1 table containing the diagnoses: example https://github.com/niekverw/ukpheno/blob/master/data/dfDefinitions.tsv

## example:
```
library(CreateUKBiobankPhentoypes)
library(readstata13) # make sure to install latest github version
library(data.table)

# set paths to data sources
# ukbiobank main dataset, converted to stata using standard scripts provided by UKB
UKbioDataset_file = "/path/to/file.dta"

# v1 hesin tables (can be loaded using LoadHesinTable (depricated, but still working)
# hesin_file="/path/to/hesin_2018-04-17.tsv"
# hesin_diagicd10_file="/path/to/hesin_diagicd10_2018-04-17.tsv"
# hesin_diagicd9_file="/path/to/hesin_diagicd9_2018-04-17.tsv"
# hesin_oper_file="/path/to/hesin_oper4_2018-04-17.tsv"

# v2 hesin tables (current)
fhesin <- "/data/scratch/project/medicinal_graphology/hes_tables/hesin.txt"
fhesin_oper <- "/data/scratch/project/medicinal_graphology/hes_tables/hesin_oper.txt"
fhesin_diag <- "/data/scratch/project/medicinal_graphology/hes_tables/hesin_diag.txt"


gp_clinical_file = "pathto/gpclinical.txt"
gp_scripts_file = "pathto/gpscripts.txt"

dfDefinitions_file = "/path/to/https://github.com/niekverw/ukpheno/blob/master/data/dfDefinitions.tsv"
Outputdir="/path/to/output"


# load data 
print("load definition table")
dfDefinitions = data.frame(fread(dfDefinitions_file))
# paste(c("n_eid",sapply(get_allvarnames(dfDefinitions_processed),function(x) paste0("*_",x,"_*"))),collapse=" ")
# write.table(ProcessDfDefinitions(dfDefinitions),paste(dfDefinitions_file,".check.tsv",sep=""),sep="\t",quote=FALSE,row.names = FALSE) # used to debug your definitions.
print("load dataframe ukbiobank")
UKbioDataset <-  as.data.frame(read.dta13(UKbioDataset_file,convert.dates = TRUE,convert.factors=F))

# names(UKbioDataset)[grepl('^s_',names(UKbioDataset) )] <- paste0("t",names(UKbioDataset)[grepl('^s_',names(UKbioDataset) )])
# UKbioDataset$ts_53_0_0 =  as.Date(UKbioDataset$ts_53_0_0, "%d%b%Y")
# UKbioDataset$ts_53_1_0 =  as.Date(UKbioDataset$ts_53_1_0, "%d%b%Y")
# UKbioDataset$ts_53_2_0 =  as.Date(UKbioDataset$ts_53_2_0, "%d%b%Y")
# UKbioDataset$ts_40000_0_0 =  as.Date(UKbioDataset$ts_40000_0_0, "%d%b%Y")

print("load hesin")
dfhesintables <- LoadHesinTable_v2(UKbioDataset,fhesin,fhesin_diag,fhesin_oper)

dfgpclinical <- loadGPTable(UKbioDataset,gp_clinical_file, 
                            cols_tokeep=c("eid","event_dt","read_2","read_3"),
                            cols_rename=c("n_eid","event_dt","read_2","read_3"))

dfgpscripts <- loadGPTable(UKbioDataset,gp_scripts_file, 
                            cols_tokeep=c("eid","issue_date","bnf_code","dmd_code"),
                            cols_rename=c("n_eid","event_dt","bnf_code","dmd_code"))


# pull out the data. 
print("constructing diagnoses for baseline visit. 
CreateUKBiobankPhentoypes(Nvisits=3,
                          visitreference=0,
                          UKbioDataset,
                          dfhesintables,
                          dfgpclinical,dfgpscripts,
                          dfDefinitions,
                          Outputdir,
                          VctOutputIndividualColumns=c("TS","SR","TS_RX","SR_RX")
                          )
```



