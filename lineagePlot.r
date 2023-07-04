#import necessary packages
library(Seurat)
library(ggplot2)
library(scales)
library(stringr)
library(RColorBrewer)
library(repr)
library(doParallel)
library("viridis")  
library(patchwork)
library(ggsci)
source("~/dualproject/lineageData/common_lineage_functions.r")
library(tidyverse)
library(ggsignif)
library(FSelectorRcpp)
library(Seurat)

#read merged datasets
EBmerge = readRDS("/data1/home/gdpeng/chengchen/dualproject/scripts/EBcombined.rds")
EBmerge@meta.data$cellBC = rownames(EBmerge@meta.data)
EBmerge@meta.data = EBmerge@meta.data %>% select(orig.ident, nCount_RNA,
                                                nFeature_RNA, S.Score, G2M.Score,
                                                Phase, cell_type, cellBC)

#color: similar cell type with similar colors, different cell type with different colors
cell_type = c('Heart','Neuron', 'Gut','Blood', 'CPM','PGC-like','Endothelium')
#colors = pal_d3()(7) #use the d3 series color

colors = c('#1F77B4FF','#FF7F0EFF','#2CA02CFF','#D62728FF','#9467BDFF','#8C564BFF','#E377C2FF')
cl_cols = setNames(colors, cell_type)

#read amplicon dataset and preproces
amplicon = "D10b"
singlecell=EBmerge@meta.data
cpf1_mutated=TRUE 
cas9_mutated=TRUE
cpf1_first=FALSE
remove_single_type=FALSE
filename1 = paste0('~/dualproject/lineageData/mannualSharedCellV5/', amplicon, 'allCassiopeiaAllele.csv')
allele_table = read.csv(filename1, row.names=1)
##preprocessing
#both lineage and cell type information
LinTrans = allele_to_lineage(allele_table, singlecell)
#preprocessing, default, keep only muatated cas9 and cpf1 scars only
LinTrans = lintrans_preprocess(LinTrans, cpf1_mutated, cas9_mutated)

#obtatin the final dataframe for information gain calculation
linIG = lintrans_to_linig(LinTrans, cl_cols,cpf1_first, remove_single_type)
head(linIG)
