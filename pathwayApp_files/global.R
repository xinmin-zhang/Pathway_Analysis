#source("renv.restore.R")
#Modified for Pathway App, in pathwayApp_files/ subfolder...

options(stringsAsFactors=F)
options(ggrepel.max.overlaps = Inf)

suppressPackageStartupMessages({
	library(shiny)
	library(shinythemes)
	library(shinyjs)
	library(plotly)
	#library(heatmaply)
	library(reshape2)
	#library(dplyr)
	#library(tidyr)
	#library(tibble)
	#library(stringr)
	#library(ggplot2)
	library(tidyverse)
	library(gplots)
	library(ggpubr)
	library(gridExtra)
	library(ggrepel)
	#library(dendextend)
	library(DT)
	#library(data.table)
	library(RColorBrewer)
	library(pheatmap)
  options(rgl.useNULL = TRUE)
	library(rgl)
	library(car)
	library(colourpicker)
	library(VennDiagram)
	#library(pathview)
	library(factoextra)
	#library(cluster)
	#library(NbClust)
	library(Mfuzz)
	library(openxlsx)
	#library(rmarkdown)
	library(visNetwork)
	library(cowplot)
	library(circlize)
	library(ComplexHeatmap)
	library(svglite)
	library(pathview)
	library(shinyjqui)
	library(ggrastr)
	library(biomaRt)
	library(networkD3)
	library(Hmisc)
	library(stringi)
	library(plyr)
	library(ggExtra)
  library(png)
  library(psych)
  library(broom)
  library(fgsea)
  #library(DEGreport)
  library(rclipboard)
})



#load("db/hgnc.RData")
#load("db/kegg.pathways.RData")
#load("db/gmtlist.RData")
#load("db/mouse_rat_genes_map2_human.RData")
homologs=readRDS("db/Homologs.rds") #cross species gene symbol mapping file from Ensembl Gene 110
#source("PC_Covariates.R")


homolog_mapping<-function(genelist, species1, species2, homologs) {
  if (species2=="human") {
    genelist2=toupper(genelist)
  } else {
    genelist2=str_to_title(genelist)
  }
  if (tolower(species1) %in% c("human", "mouse", "rat")) {
    species1=tolower(species1)
    lookup<-homologs%>%filter(Species1==species1, Species2==species2)
    df<-data.frame(genelist, genelist2)%>%left_join(lookup%>%transmute(genelist=symbol1, mapped_symbol=symbol2))%>%
      mutate(mapped_symbol=ifelse(is.na(mapped_symbol), genelist2, mapped_symbol))
    genelist2<-df$mapped_symbol
  }
  return(genelist2)
}

options(shiny.maxRequestSize = 40*1024^2)  #upload files up to 30 Mb


mycss <- "select ~ .selectize-control .selectize-input {
        max-height: 200px;
        overflow-y: auto;
        }
    .shiny-notification{
    position: fixed;
    top: 33%;
    left: 33%;
    right: 33%;
    }"

saved_projects = read.csv("data/saved_projects.csv")
projects = saved_projects$ProjectID
names(projects) = saved_projects$ShortNames

pub_projects=NULL
if (file.exists("unlisted/projects_public.csv")) { #load selected projects from unlisted folder for Venn across projects
  Pinfo=read.csv("unlisted/projects_public.csv")
  pub_projects=Pinfo$ProjectID
  names(pub_projects)=Pinfo$ShortName
}




footer_text = '
<hr>
<div align="center" style="font-size:11px">Pathway Tools Developed by:
Benbo Gao, Xinmin Zhang and Baohong Zhang
</div>
'

config=NULL
server_dir=NULL
test_dir=NULL
gmt_file_info=NULL
system_info=NULL
if (file.exists("config_PV.csv")) { #load optional configuration file
  config=read_csv("config_PV.csv")
  N=match("server_dir", config$category)
  if (!is.na(N)) {server_dir=config$value[N]}
  N=match("test_dir", config$category)
  if (!is.na(N)) {test_dir=config$value[N]}
  N=match("geneset_api", config$category) #if using internal api for geneset
  if (!is.na(N)) {footer_text=str_replace_all(footer_text, "//bxngs.com/bxomics/api", config$value[N])}
  N=match("gmt_file_info", config$category)
  if (!is.na(N)) {gmt_file_info=config$value[N]}
  N=match("system_info", config$category)
  if (!is.na(N)) {system_info=config$value[N]}
  #browser() #debug
}
