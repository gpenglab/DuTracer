#preprocess
#read the allele_table and related parameters
amplicon = "D8a"
singlecell=EBmerge@meta.data
cpf1_mutated=FALSE #also calcute cpf_mutated=TRUE, not shown here
cas9_mutated=FALSE
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
linIG = linIG[,-1]

#check the observed normalized_ig
celltype_entropy <- sum(-prop.table(table(linIG$cell_type)) * log2(prop.table(table(linIG$cell_type)))) 

ig_observed = information_gain(formula = cell_type ~ cpf1, linIG) #only check the different ablitiy in the last level
        ig_observed$importance = (log2(exp(1)))*(ig_observed$importance)        
        if (celltype_entropy>0) {
        ig_observed$normalized_ig = ig_observed$importance/celltype_entropy
        } else {
        ig_observed$normalized_ig = 0
        }


#perform simulation to check the distribution
simulated_igs = vector()
linIG_simulated = linIG
n_sim = 1:1000
for (i in n_sim) {
        numbers = nrow(linIG)
        linIG_simulated$cell_type = sample(linIG$cell_type, numbers, replace=F) 
        ig = information_gain(formula = cell_type ~ cpf1, linIG_simulated) #remove grp column
        ig$importance = (log2(exp(1)))*(ig$importance)        
        if (celltype_entropy>0) {
        ig$normalized_ig = ig$importance/celltype_entropy
        } else {
        ig$normalized_ig = 0
        }
        simulated_igs = append(simulated_igs, ig$normalized_ig)
} 
simulated_igs_table = data.frame(n_sim, simulated_igs)  
