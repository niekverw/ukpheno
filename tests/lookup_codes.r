
library(readxl)
fcoding.xls="data/all_lkps_maps.xlsx"
dfCodesheet.icd9_icd10 <- as.data.frame(read_xlsx(fcoding.xls,sheet="icd9_icd10"))
dfCodesheet.read_v2_icd9 <- as.data.frame(read_xlsx(fcoding.xls,sheet="read_v2_icd9"))
dfCodesheet.read_v2_icd10 <- as.data.frame(read_xlsx(fcoding.xls,sheet="read_v2_icd10"))

dfCodesheet.read_v2_lkp <- as.data.frame(read_xlsx(fcoding.xls,sheet="read_v2_lkp"))
dfCodesheet.read_v2_lkp <- as.data.frame(dfCodesheet.read_v2_lkp %>% arrange(read_code,term_code))
#dfCodesheet.read_v2_lkp <- dfCodesheet.read_v2_lkp[dfCodesheet.read_v2_lkp$term_code==0,]
dfCodesheet.read_v2_drugs_lkp<- as.data.frame(read_xlsx(fcoding.xls,sheet="read_v2_drugs_lkp"))
dfCodesheet.read_v2_lkp <- rbind(dfCodesheet.read_v2_lkp[,c(1,3)],dfCodesheet.read_v2_drugs_lkp[,1:2])
dfCodesheet.read_v2_lkp <- as.data.frame(dfCodesheet.read_v2_lkp %>% group_by(read_code) %>%  summarize(text = str_c(term_description, collapse = "/")))
#dfCodesheet.read_v2_lkp$read_code <- gsub("\\.","", dfCodesheet.read_v2_lkp$read_code)
dfCodesheet.read_v2_lkp <- dfCodesheet.read_v2_lkp %>% unique()  %>% arrange(read_code)

dfCodesheet.icd9_icd10 <- as.data.frame(read_xlsx(fcoding.xls,sheet="icd9_icd10"))


# ICD9 depscription, dfCodesheet.icd9_lkp
dfCodesheet.icd9_lkp <- as.data.frame(read_xlsx(fcoding.xls,sheet="icd9_lkp")) # certainly not complete!
dfCodesheet.ICD9.coding87 <- data.frame(fread("data/ICD9.coding87.tsv"))[,1:2]
names(dfCodesheet.ICD9.coding87)<- c("ICD9","DESCRIPTION_ICD9")
dfCodesheet.icd9_lkp <- rbind(dfCodesheet.ICD9.coding87[,1:2],dfCodesheet.icd9_lkp[!dfCodesheet.icd9_lkp$ICD9 %in% dfCodesheet.ICD9.coding87$ICD9 ,])  %>% unique()  %>% arrange(ICD9)


dfCodesheet.icd10_lkp <- as.data.frame(read_xlsx(fcoding.xls,sheet="icd10_lkp"))[,c(2,5)] # not complete (e.g. X*)
names(dfCodesheet.icd10_lkp) <- c("ICD10","DESCRIPTION")
dfCodesheet.ICD10.coding19 <- data.frame(fread("data/ICD10.coding19.tsv"))[,1:2]
names(dfCodesheet.ICD10.coding19) <- c("ICD10","DESCRIPTION")
dfCodesheet.icd10_lkp <- rbind(dfCodesheet.ICD10.coding19[,1:2],dfCodesheet.icd10_lkp[!dfCodesheet.icd10_lkp$ICD10 %in% dfCodesheet.ICD10.coding19$coding ,]) %>% unique() %>% arrange(ICD10)


dfCodesheet.OPCS4.coding240 <- data.frame(fread("data/OPCS4.coding240.tsv"))
dfCodesheet.opcs4_lkp<- dfCodesheet.OPCS4.coding240

#E8801
#dfCodesheet.read_v2_lkp[duplicated(dfCodesheet.read_v2_lkp$read_code)==TRUE,]
unique(dfCodesheet.read_v2_lkp$read_code)
#####


############################################################

library(stringr)

print("load definition table")
df = data.frame(fread(dfDefinitions_file))

#' @export
CovertMednamesToUkbcoding<- function(StrRx){
  #StrRx<-"phenformin,metformin,buformin,glibenclamide,chlorpropamide,tolbutamide,glibornuride,tolazamide,carbutamide,glipizide,gliquidone,gliclazide,metahexamide,glisoxepide,glimepiride,acetohexamide,glymidine,acarbose,miglitol,voglibose,troglitazone,rosiglitazone,pioglitazone,sitagliptin,vildagliptin,saxagliptin,alogliptin,linagliptin,gemigliptin,repaglinide,nateglinide,exenatide,pramlintide,benfluorex,liraglutide,mitiglinide,dapagliflozin,lixisenatide,canagliflozin,empagliflozin,albiglutide,dulaglutide"
  StrRx<-as.character(StrRx)
  if(is.na(StrRx)) { return(NA)}
  VctRXstrings<-unlist(strsplit(StrRx,","))
  #VctRXstrings<-strsplit(df[!is.na(df$n_20003_),]$n_20003_,",")[[1]]
  StrRxCodes<-paste(unique(unlist(lapply(VctRXstrings,  function(x) dfCodesheetREAD_SR.Coding[,"UKB.Coding"] [ grep(x,dfCodesheetREAD_SR.Coding[,"Meaning"] ,ignore.case=TRUE )]  ))),collapse=",")
  return(StrRxCodes)
}

#' @export
convert.coding<- function(Str,
                          from.code="READ.CODE",
                          to.code="UKB.Coding",
                          lookuptable=dfCodesheetREAD_SR.Coding,ignore.case=FALSE){
  # Str<-"f3,f4,ft"
  Str<-as.character(Str)
  if(is.na(Str)) { return(NA)}
  VctStr<-unlist(strsplit(Str,","))
  #VctRXstrings<-strsplit(df[!is.na(df$n_20003_),]$n_20003_,",")[[1]]
  c <- paste(unique(unlist(lapply(VctStr,  function(x)
    lookuptable[,to.code] [ grep(paste("^", x,sep=""),lookuptable[,from.code] ,ignore.case=ignore.case )]
    ))),collapse=",")

  return(c)
}

add.description.to.codes <- function(Str,code.id="UKB.Coding",description.id="Meaning",description.lookuptable=dfCodesheetREAD_SR.Coding,ignore.case=FALSE,firstcodeonly=TRUE) {
  if(is.na(Str) | Str =="NA"){return(Str)}
  Str<-as.character(Str)
  if(is.na(Str)) { return(NA)}
  VctStr<-unlist(strsplit(Str,","))


  c<- sapply(VctStr,  function(x){
    x.d <- description.lookuptable[,description.id] [ grep(paste("^", x,sep=""),description.lookuptable[,code.id] ,ignore.case=ignore.case )]
    if(length(x.d)==0){return(paste0(x," (NA)"))}
    x.d <- str_replace_all(x.d,  "[^/[:^punct:]]", "") # replace all symbols to not mess up downstream things.
    if(firstcodeonly==TRUE){
      x.d <- x.d[1]
      }
    x.d <- paste0(x.d,collapse=" /")
    x.d <- paste0(x," (",x.d,")")
    x.d
  },USE.NAMES = F)

  c <- paste(unique(unname(c)),collapse=",")
  return(c)
}



convert_definition_column <- function(source_col=df$READCODES,target_col=df$n_20003_,
                                      from.code="READ.CODE",to.code="UKB.Coding",lookuptable=dfCodesheetREAD_SR.Coding,
                                      description.code.id=NULL,description.id="Meaning",description.lookuptable=NULL) {

  source_col <- PreProcessDfDefinitions(  data.frame(source_col=source_col,tmp=rep("NA",length(source_col))),VctAllColumns = c("source_col","tmp"))[,1]
  if(!is.null(target_col)){
    target_col <- PreProcessDfDefinitions(  data.frame(target_col=target_col,tmp=rep("NA",length(target_col))),VctAllColumns = c("target_col","tmp"))[,1]
  } else{
    target_col <- rep("NA",length(source_col))
  }

  c1 <- paste(target_col, unlist(lapply(source_col, convert.coding,from.code=from.code,to.code=to.code,lookuptable=lookuptable)),sep=",")
  c2 <- unlist(lapply(c1,function(x) {  x = unique(strsplit(x,"," )[[1]]); if(length(x)==1 & x[1] =="NA"){ return("NA")} else{ return( paste(x[x != "NA"],collapse=",") )} }))
  # add description?
  if(is.null(description.lookuptable)){description.lookuptable=lookuptable}
  if(is.null(code.id)) { code.id =to.code }
  if(!is.null(description.id)){
    #code.id=to.code
    c2 <- sapply( c2, add.description.to.codes,
                  code.id=code.id,
                  description.id=description.id,
                  description.lookuptable=description.lookuptable,USE.NAMES = F)

  }
  return(c2)
}

expand_clean_codes <- function(col=df$ICD10CODES, from.code="ALT_CODE",description.id='DESCRIPTION',lookuptable = dfCodesheet.icd10_lkp,add_description=T){
  # col=df$ICD10CODES
  # lookuptable=dfCodesheet.icd10_lkp
  # from.code="ALT_CODE"
  # description.id='DESCRIPTION'

  col <- PreProcessDfDefinitions(  data.frame(col=col,tmp=rep("NA",length(col))),VctAllColumns = c("col","tmp"))[,1]

  to.code="self"
  lookuptable$self <- lookuptable[,from.code]
  #c <- unlist(lapply(col, convert.coding,from.code=from.code,to.code=to.code,lookuptable=lookuptable))

  c1 <- paste(col, unlist(lapply(col, convert.coding,from.code=from.code,to.code=to.code,lookuptable=lookuptable)),sep=",")
  c2 <- unlist(lapply(c1,function(x) {  x = unique(strsplit(x,"," )[[1]]); if(length(x)==1 & x[1] =="NA"){ return("NA")} else{ return( paste(x[x != "NA"],collapse=",") )} }))

  if(add_description==T){
    c2 <- sapply( c2, add.description.to.codes,
                  code.id=from.code,
                  description.id=description.id,
                  description.lookuptable=lookuptable,USE.NAMES = F)
  }
  return(c2)
}


df = data.frame(fread(dfDefinitions_file))
df$n_20003_ <- convert_definition_column(source_col = df$READCODES,
                                         target_col = df$n_20003_,
                                         fromcode="READ.CODE",to.code="UKB.Coding",description.id="Meaning")#,lookuptable = dfCodesheetREAD_SR.Coding)

convert_definition_column(source_col = df$n_20003_ ,
                          target_col = df$n_20003_,
                          from.code="UKB.Coding",to.code="UKB.Coding",lookuptable = dfCodesheetREAD_SR.Coding,
                          description.code=NULL,description.id="Meaning",description.lookuptable=NULL)


convert_definition_column(source_col = df$ICD10CODES ,
                          target_col = df$READCODES,
                          from.code="icd10_code",to.code="read_code",lookuptable = dfCodesheet.read_v2_icd10,
                          description.code=NULL,description.id="text",description.lookuptable=dfCodesheet.read_v2_lkp) #,lookuptable = dfCodesheetREAD_SR.Coding)


convert_definition_column(source_col = df$ICD9CODES ,
                          target_col = df$READCODES,
                          from.code="icd9_code",to.code="read_code",lookuptable = dfCodesheet.read_v2_icd9,
                          description.code=NULL,description.id="text",description.lookuptable=dfCodesheet.read_v2_lkp) #,lookuptable = dfCodesheetREAD_SR.Coding)

convert_definition_column(source_col = df$ICD10CODES ,
                          target_col = df$ICD9CODES,
                          from.code="icd10_code",to.code="icd9_code",lookuptable = dfCodesheet.read_v2_icd9,
                          description.code=NULL,description.id="text",description.lookuptable=dfCodesheet.read_v2_lkp) #,lookuptable = dfCodesheetREAD_SR.Coding)



# expand & clean icd9/10 codes.
expand_clean_codes(col =df$ICD10CODES, from.code="ICD10",description.id='DESCRIPTION',lookuptable = dfCodesheet.icd10_lkp,add_description=T)
expand_clean_codes(col =df$ICD9CODES, from.code="ICD9",description.id='DESCRIPTION_ICD9',lookuptable = dfCodesheet.icd9_lkp,add_description=T)
expand_clean_codes(col =df$OPERCODES, from.code="coding",description.id='meaning',lookuptable = dfCodesheet.opcs4_lkp,add_description=T)
expand_clean_codes(col =df$READCODES, from.code="read_code",description.id='text',lookuptable = dfCodesheet.read_v2_lkp,add_description=T)


# suggest READ codes based on all ICD10,9,oper,read,CTVT ....

col=expand_clean_codes(col =df$ICD10CODES, from.code="ICD10",description.id='DESCRIPTION',lookuptable = dfCodesheet.icd10_lkp,add_description=T)
icd10read<-convert_definition_column(source_col = col,
                          target_col = NULL,
                          from.code="icd10_code",to.code="read_code",lookuptable = dfCodesheet.read_v2_icd10,
                          description.code=NULL,description.id="term_description",description.lookuptable=NULL) #,lookuptable = dfCodesheetREAD_SR.Coding)

col=expand_clean_codes(col =df$ICD9CODES, from.code="ICD9",description.id='DESCRIPTION_ICD9',lookuptable = dfCodesheet.icd9_lkp,add_description=T)
icd9read <- convert_definition_column(source_col = col,
                                target_col = NULL,
                                from.code="icd9_code",to.code="read_code",lookuptable = dfCodesheet.read_v2_icd9,
                                description.code=NULL,description.id="text",description.lookuptable=dfCodesheet.read_v2_lkp) #,lookuptable = dfCodesheetREAD_SR.Coding)

dfCodesheet.read_v2_icd10[grep("-",dfCodesheet.read_v2_icd10$icd10_code),]


## collapse codes.


#######################################################################
#
dfDefs <- dfDefinitions
colnames(dfDefs) <- colnames(ProcessDfDefinitions(dfDefs,fill_dependencies = F))

library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Hello Shiny!"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
      #tabsetPanel(
        tabPanel("select", DT::dataTableOutput("oxids"))
       # )
      ),
    mainPanel(
      tabsetPanel(
        tabPanel("original ", DT::dataTableOutput("ox")),
        tabPanel("cleaned",  DT::dataTableOutput("ox_clean"))
      )
    )
  )
)


server <- function(input, output, session) {
  rdef <- reactiveValues(dforigin=dfDefs,
                         dforiginsliced=dfDefs,
                         dfclean=ProcessDfDefinitions(dfDefs,fill_dependencies = F))

  #y <- function() x
  #x$Date = Sys.time() + seq_len(nrow(x))
  output$oxids = DT::renderDataTable({
    rdef$dforigin[,c("TRAIT","DESCRIPTION")]
  },options = list(pageLength = 100), editable = TRUE, selection='single') #,class = 'nowrap stripe compact',selection = 'none',

  output$ox = DT::renderDataTable({
    rdef$dforiginsliced
    },options = list(pageLength = 100), editable = TRUE, selection='single') #,class = 'nowrap stripe compact',selection = 'none',

  # output$ox = DT::renderDataTable({
  #   x
  #   },options = list(pageLength = 100), editable = TRUE, selection='single') #,class = 'nowrap stripe compact',selection = 'none',
  #
  # proxy = DT::dataTableProxy('ox')
  #
  # observeEvent(input$OdfDefinitions_original_cell_edit, {
  #   info = input$OdfDefinitions_original_cell_edit
  #   i = info$row
  #   j = info$col
  #   v = info$value
  #   x[i, j] <<- DT::coerceValue(v, x[i, j])
  #   DT::replaceData(proxy, x, resetPaging = FALSE)  # important
  # })


  #
  observeEvent(input$oxids_rows_selected,  ignoreInit=TRUE,{
    i <- input$oxids_rows_selected
    print(i)
    rdef$dfclean <<- t(ProcessDfDefinitions(rdef$dforigin[i,],fill_dependencies = F))
    rdef$dforiginsliced <<- t(rdef$dforigin[i,])
  })

  output$ox_clean = DT::renderDataTable({
    rdef$dfclean
  },options = list(pageLength = 100), selection = 'none', editable = FALSE)


}

shinyApp(ui, server)


