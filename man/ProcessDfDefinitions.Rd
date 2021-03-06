% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ProcessdfDefinitions.R
\name{ProcessDfDefinitions}
\alias{ProcessDfDefinitions}
\title{ProcessDfDefinitions}
\usage{
ProcessDfDefinitions(df, VctAllColumns = c("TS", "SR", "TS_RX", "SR_RX",
  "LAB", "ICD10CODES", "ICD9CODES", "OPCS4CODES", "OPCS3CODES",
  "TS_AGE_DIAG_COLNAME", "READCODES", "CTV3CODES", "BNFCODES", "DMDCODES",
  "n_20001_", "n_20002_", "n_20003_", "n_20004_", "DEPENDENCY"),
  VctColstoupper = c("ICD10CODES", "ICD9CODES", "OPCS4CODES",
  "OPCS3CODES"), fill_dependencies = T)
}
\arguments{
\item{VctAllColumns}{Vct}

\item{dfDefinitions}{df}
}
\value{
None
}
\description{
Process definitions, for input
}
\examples{
#
#This function processes an excel file with definitions and is automtically performed in CreateUKBiobankPhentoypes().
#It can be usefull to run this function as a check prior to running CreateUKBiobankPhentoypes.
#
#VctAllColumns contains all column names of interest, so that it can ignore everything else.
#20001, 20002 and 20004 go into SR
#READCODES and 20003 is parsed into RX
#
#
#
VctAllColumns<-  c("TS", "SR", "TS_RX", "SR_RX", "LAB", "ICD10CODES", "ICD9CODES", "OPCS4CODES","OPCS3CODES", "TS_AGE_DIAG_COLNAME", "READCODES","CTV3CODES","BNFCODES","DMDCODES", "n_20001_",    "n_20002_", "n_20003_", "n_20004_", "DEPENDENCY")
ProcessDfDefinitions(dfDefinitions,VctAllColumns)

}
\keyword{CreateUKBiobankPhentoypes}
\keyword{ExtractVarsFromMasterSet}
\keyword{ProcessDfDefinitions}
