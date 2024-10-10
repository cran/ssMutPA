## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ssMutPA)

## ----echo = TRUE, results = 'hide',eval=FALSE---------------------------------
#  install.packages("ssMutPA")
#  library(ssMutPA)

## ----warning=TRUE, paged.print=TRUE-------------------------------------------
#load the mutation annotation file
mut_path <- system.file("extdata","mutation_data.Rdata",package = "ssMutPA")
load(mut_path)
#perform the function 'get_mut_status'
mut_status<-get_mut_status(mutation_data,nonsynonymous=TRUE,TCGA=TRUE,mut_rate=0)
#view the first five lines of mut_status matrix
mut_status[1:5,1:5]

## ----warning=TRUE, paged.print=TRUE-------------------------------------------
#Method of obtaining data
data(mut_status)
net_path <- system.file("extdata","ppi_network.Rdata",package = "ssMutPA")
load(net_path)
pathway_path<-system.file("extdata","kegg_323_gmt.Rdata",package = "ssMutPA")
load(pathway_path)
samp_name<-c("TCGA-32-1979-01A","TCGA-32-2494-01A")
examp_data<-mut_status[,samp_name]
#perform the function 'get_RWR_ES'
Path_ES<-get_RWR_ES(examp_data,net_data=ppi_network,pathway_data=kegg_323_gmt,BC_Num=12436)
#view the first six lines of pathway enrichment profile
head(Path_ES)

## ----warning=TRUE, paged.print=TRUE-------------------------------------------
#Load sample mutation data
surv_path <- system.file("extdata","sur.Rdata",package = "ssMutPA")
load(surv_path)
data(Path_ES)
#Perform the function `get_samp_class`
res<-get_samp_class(Path_ES,sur,seed_num=5,cox_pval=0.05,min.nc = 2,max.nc =5)
#view the label of samples
res$sample_class[1:10]

## ----fig.height=6, fig.width=8,warning=FALSE,results='hold'-------------------
#Load the data
data(Path_ES,sample_class,Path_Name)
#perform the function `get_heatmap`.
get_heatmap(Path_ES,Path_name=Path_Name,samp_class=sample_class)

## ----fig.height=6, fig.width=8, warning=FALSE, results='hold'-----------------
#Get the data of ROC curve
data(Path_ES,sample_class)
#perform the function `mountain_plot`
mountain_plot(data=Path_ES,sample_class=sample_class,Path_name=rownames(Path_ES)[c(12,20,74,103,113,123,138,151,188)])

## ----fig.height=6, fig.width=8,warning=FALSE,results='hold'-------------------
#load the data
data(dot_data)
#perform the function `dotplot`.
dotplot(dot_data)

## ----fig.height=4, fig.width=8,warning=FALSE,results='hold'-------------------
#load the data
mut_path <- system.file("extdata","maffile.txt",package = "ssMutPA")
maf<-maftools::read.maf(mut_path ,isTCGA = FALSE)
pathway_path <- system.file("extdata","kegg_323_gmt.Rdata",package = "ssMutPA")
load(pathway_path)
data(samp_class_onco,mut_onco,sur_onco)
#draw a waterfall plot
# \donttest{Oncoplot(maf,samp_class_onco,sur_onco,mut_onco,kegg_323_gmt,"IL-17 signaling pathway")}

