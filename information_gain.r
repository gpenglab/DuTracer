#This script provide basic function to calculate information gain and normalized information gain from lineage plots
amplicon = "D8b"
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
linIGx = linIG[,-1]

#remove cell information and reorder the columns to the actual gene editing orders
linIGx = linIGx[, c("grp","cas9","cpf1","cell_type")]

##check the informataiton gain within one sample
#Note: information_gain calculate from FSelectorRcpp package use nat as basic unit, should transform into bit as unit
#information gain calculate the explained cell type shannon, which is equivalent the explained cell differentiation
#check the information gain change for each clone, each level
#define a function to obtain both information gain and cell type information gain ratio

ig_grp_within_sample = function(linIGx) {
    linGrps_split = split(linIGx, list(linIGx$grp))
    linGrps_igs = do.call(rbind, lapply(linGrps_split, function(linGrp) {
        ig = information_gain(formula = cell_type ~ ., linGrp[,-1]) #remove grp column
        ig$importance = (log2(exp(1)))*(ig$importance)
        ig$linGrp = linGrp$grp[1]
        celltype_entropy <- sum(-prop.table(table(linGrp[,"cell_type"])) * log2(prop.table(table(linGrp[,"cell_type"])))) 
        if (celltype_entropy>0) {
            ig$celltype_ratio = ig$importance/celltype_entropy
        } else {
            ig$celltype_ratio = 0
        }
        return(ig)
} ))
    return(linGrps_igs)
}

##check the information gains of different samples
amplicons = c("D8a","D10a","D8b","D10b","D12","D14")

#NOTE: the barcode entropy here are level-wise, instead of calculate independently
#NOTE: currently, we only calculate the barcodeset instead of calculate individual barcodesets
#be careful with the column name order

#data filtering requirement
singlecell=EBmerge@meta.data
cpf1_mutated=TRUE 
cas9_mutated=TRUE
cpf1_first=FALSE #change the order of mutation events
#cpf1_first=TRUE
remove_single_type=FALSE

igs = data.frame(attributes=character(), importance=double(),celltype_ratio=double(), BC_ratio=double(),
                BCentropy=double(),sample=character())

for (amplicon in amplicons) {
    filename1 = paste0('~/dualproject/lineageData/mannualSharedCellV5/', amplicon, 'allCassiopeiaAllele.csv')
    allele_table = read.csv(filename1, row.names=1)
    ##preprocessing
    #both lineage and cell type information
    LinTrans = allele_to_lineage(allele_table, singlecell)
    #preprocessing, default, keep only muatated cas9 and cpf1 scars only
    LinTrans = lintrans_preprocess(LinTrans, cpf1_mutated, cas9_mutated)

    #obtatin the final dataframe for information gain calculation
    linIG = lintrans_to_linig(LinTrans, cl_cols,cpf1_first, remove_single_type)
    linIG = linIG[,-1]
    ig = information_gain(formula = cell_type ~ ., linIG)
    ig1 = information_gain(formula = cell_type ~ ., linIG, type='gainratio')
    ig$importance = (log2(exp(1)))*(ig$importance)
    celltype_entropy <- sum(-prop.table(table(linIG[,"cell_type"])) * log2(prop.table(table(linIG[,"cell_type"]))))

    if (celltype_entropy>0) {
        ig$celltype_ratio = ig$importance/celltype_entropy #normalized by celltype entropy
    } else {
        ig$celltype_ratio = 0
    }
    ig$BC_ratio = ig1$importance #gain ratio of barcodes,normalized by barcode entropy

    grpBCentropy =  sum(-prop.table(table(linIG[,"grp"])) * log2(prop.table(table(linIG[,"grp"])))) #actually all barcode combinations
    cas9BCentropy =  sum(-prop.table(table(linIG[,"cas9"])) * log2(prop.table(table(linIG[,"cas9"]))))
    cpf1BCentropy =  sum(-prop.table(table(linIG[,"cpf1"])) * log2(prop.table(table(linIG[,"cpf1"]))))
    ig$BCentropy = c(grpBCentropy, cpf1BCentropy, cas9BCentropy)

    ig$sample = amplicon
    igs = rbind(igs, ig)
}
colnames(igs) = c("bc_type", "information_gain", "normalized_ig", "normalizedBCig", "BCentropy","samples")
igs$samples = factor(igs$samples, levels=c("D8a","D10a","D8b","D10b","D12","D14"))
igs$bc_type = factor(igs$bc_type, levels=c("grp","cas9","cpf1"))
