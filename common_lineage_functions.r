#waiting for adding more description and printing examples and neccessary hints
#some functions are required further improvement to ensure idenpent usages
#some functions are tedious and ugly, need improvement

#PART1: miscellaneous small scripts
#--------------------------------------------------------------
#calculate the shannon diversity for the distribution of different categories
#input: a column with different categories
#output: shannon diversity
shannon = function (freq) {
    sum_freq = sum(freq)
    shan = 0
    for (i in freq) {
        index = -(i/sum_freq)*log2(i/sum_freq)
        shan = shan+index
    }
    return(shan)
}

#a summarge function pass to summarise function to aggregate all scar information
scar_merge = function(x) {
    indel = ""
    for (i in x) {
        indel  = paste(indel, i, sep="|")
    }
    return(indel)
}
scar_ig_merge = function(x) {
    indel = ""
    for (i in x) {
        indel  = paste(indel, i, sep="-")
    }
    return(indel)
}

#PART2: some functions to handle calculation of information gains
#--------------------------------------------------------------
#obtain both the lineage information and barcode information
#input: allele_table of DuTracer and single cell meta data
#output: a dataframe contain both single cell meta information and lineage information
allele_to_lineage = function(allele_table, sc_meta) {
    LinTrans = inner_join(allele_table, sc_meta, by="cellBC")
    #LinTrans = LinTrans[grep("[ID]", LinTrans$allele),] will introduce significant biase to the dataset
    LinTrans = LinTrans %>% select(cellBC, intBC, r1, r2, r3, r4, lineageGrp, cell_type)

    #transform the integration barcode for better representation
    intBCnames = sapply(1:length(unique(LinTrans$intBC)), function(x) paste0("bc",x))
    setIntBCnames = setNames(intBCnames, unique(LinTrans$intBC))             
    LinTrans$intBC = setIntBCnames[LinTrans$intBC]

    #transform the target site for better representation
    LinTrans$r1 = gsub("[ATCG]", "", LinTrans$r1)
    LinTrans$r2 = gsub("[ATCG]", "", LinTrans$r2)
    LinTrans$r3 = gsub("[ATCG]", "", LinTrans$r3)
    LinTrans$r4 = gsub("[ATCG]", "", LinTrans$r4)
    return(LinTrans)
}

#--------------------------------------------------------------
#select_barcode = "CTACGCTTGTCCGA[91:12D]"
#define a function to calculate the information gain of a specific barcode
#target_barcode is a data.frame for specific target set
info_gain_barcode = function(target_barcode, select_barcode) {
     #calculate those shannons with selected barcodes
    selected_barcode_cells = target_barcode %>% filter( grepl(select_barcode, barcodeset, fixed=T) ) %>% 
        select(type) %>% table() %>% as.data.frame()
    positive = sum(selected_barcode_cells$Freq)
    positive_shannon = shannon(selected_barcode_cells$Freq[selected_barcode_cells$Freq != 0])
     
     #calculate those shannons without selected barcodes
    selected_barcode_cells = target_barcode %>% filter( !grepl(select_barcode, barcodeset, fixed=T) ) %>% 
        select(type) %>% table() %>% as.data.frame()
    negative = sum(selected_barcode_cells$Freq)
    negative_shannon = shannon(selected_barcode_cells$Freq[selected_barcode_cells$Freq != 0])
    son_shannon = positive_shannon*positive/nrow(target_barcode)+negative_shannon*negative/nrow(target_barcode)
    overall_barcode_cells = as.data.frame(table(target_barcode$type))
    father_shannon = shannon(overall_barcode_cells$Freq[overall_barcode_cells$Freq != 0])
    info_gain = round(father_shannon-son_shannon, 5)
    return(list(select_barcode, info_gain))
}


#--------------------------------------------------------------
#select_barcodes should contains both mutated and non-mutated
#In order to compare different separation effect of different targets, calculate the total inforamtion gain using all cells irrespective of their lineage group origin
#Note: LinTrans should be a dataframe for a specific sample, while linGrp should be a specific number representing a lineage group
IG_all__barcodes_all_cells = function(LinTrans, cpf1_mutated=FALSE, cas9_mutated=FALSE) {
    #filter those allele and cell type information for a specific lineage group
    lineage_selected = LinTrans %>% mutate(alltarget=paste("LG", lineageGrp, intBC, r1, r2, r3, r4, sep="_")) %>% 
          mutate(cas9 = paste("LG", lineageGrp, intBC, r3, r4, sep="_"), cpf1=paste("LG", lineageGrp, intBC, r1, r2, sep="_")) %>% 
          mutate(r1 = paste("LG", lineageGrp, intBC, r1, sep="_"), r2 = paste("LG", lineageGrp,intBC, r2, sep="_"), 
          r3=paste("LG", lineageGrp,intBC, r3, sep="_"), r4=paste("LG", lineageGrp,intBC, r4, sep="_")) 
    if (cas9_mutated) {
        lineage_selected = lineage_selected[grep("[ID]",lineage_selected$cas9),]
    }
    #if (cpf1_mutated=="mutated") {
    if (cpf1_mutated) {
        lineage_selected = lineage_selected[grep("[ID]", lineage_selected$cpf1),]
    }

    #waiting to achieve remove those clones with only one cell types

    #calculate father shannon
    target_barcode = lineage_selected %>% select(cellBC, r1, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r1), type = unique(cell_type))
    overall_barcode_cells = as.data.frame(table(target_barcode$type))
    father_shannon = shannon(overall_barcode_cells$Freq[overall_barcode_cells$Freq != 0])

    #calculate target r1, only considering the information_gain in all cells without requiring mutated cpf1 or mutated cas9
    target_barcode = lineage_selected %>% select(cellBC, r1, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r1), type = unique(cell_type))
    #checking the barcode according to single integration barcode
    #selected those barcodes with both mutation and no-mutation to calculate information gain
    #barcodes = unique(grep("[ID]", lineage_selected$r1, value=TRUE))
    barcodes = unique(lineage_selected$r1)
    info_gain_r1 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_r1 = rbind(info_gain_r1, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_r1) = c('barcode', 'info_gain')

    #calculate target r2
    target_barcode = lineage_selected %>% select(cellBC, r2, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r2), type = unique(cell_type))

    barcodes = unique(lineage_selected$r2)
    info_gain_r2 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_r2 = rbind(info_gain_r2, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_r2) = c('barcode', 'info_gain')

    #calculate target r3
    target_barcode = lineage_selected %>% select(cellBC, r3, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r3), type = unique(cell_type))

    barcodes = unique(lineage_selected$r3)
    info_gain_r3 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_r3 = rbind(info_gain_r3, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_r3) = c('barcode', 'info_gain')

    #calculate target r4
    target_barcode = lineage_selected %>% select(cellBC, r4, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r4), type = unique(cell_type))

    barcodes = unique(lineage_selected$r4)
    info_gain_r4 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_r4 = rbind(info_gain_r4, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_r4) = c('barcode', 'info_gain')

    #calculate target cas9
    target_barcode = lineage_selected %>% select(cellBC, cas9, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(cas9), type = unique(cell_type))

    barcodes = unique(lineage_selected$cas9)
    info_gain_cas9 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_cas9 = rbind(info_gain_cas9, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_cas9) = c('barcode', 'info_gain')

    #calculate target cpf1
    target_barcode = lineage_selected %>% select(cellBC, cpf1, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(cpf1), type = unique(cell_type))

    barcodes = unique(lineage_selected$cpf1)
    info_gain_cpf1 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_cpf1 = rbind(info_gain_cpf1, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_cpf1) = c('barcode', 'info_gain')

    #calculate all target
    target_barcode = lineage_selected %>% select(cellBC, alltarget, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(alltarget), type = unique(cell_type))

    barcodes = unique(lineage_selected$alltarget)
    info_gain_alltarget = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_alltarget = rbind(info_gain_alltarget, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_alltarget) = c('barcode', 'info_gain')

    #add target label for separte plot
    if (nrow(info_gain_r1)>0) {
        info_gain_r1$target = "r1"
    }
    if (nrow(info_gain_r2)>0) {
        info_gain_r2$target = "r2"
    }
    if (nrow(info_gain_r3)>0) {
        info_gain_r3$target = "r3"
    }
    if (nrow(info_gain_r4)>0) {
        info_gain_r4$target = "r4"
    }
    if (nrow(info_gain_cpf1)>0) {
        info_gain_cpf1$target = "cpf1"
    }
    if (nrow(info_gain_cas9)>0) {
        info_gain_cas9$target = "cas9"
    }
    if (nrow(info_gain_alltarget)>0) {
        info_gain_alltarget$target = "alltarget"
    }

    #build a datafram to store all infomation gain for four targets
    info_gain_r1 = info_gain_r1 %>% arrange(info_gain)
    info_gain_r2 = info_gain_r2 %>% arrange(info_gain)
    info_gain_r3 = info_gain_r3 %>% arrange(info_gain)
    info_gain_r4 = info_gain_r4 %>% arrange(info_gain)
    
    info_gain_all = rbind(info_gain_r1, info_gain_r2, info_gain_r3, info_gain_r4)
    info_gain_all = info_gain_all[!duplicated(info_gain_all$barcode),]
    info_gain_all$target = factor(info_gain_all$target, levels=c("r1", "r2", "r3", "r4"))
    info_gain_all$barcode = factor(info_gain_all$barcode, levels=info_gain_all$barcode)
    #build a data.frame to store all information gain for either cas9 or cpf1
    info_gain_cas9 = info_gain_cas9 %>% arrange(info_gain)
    info_gain_cpf1 = info_gain_cpf1 %>% arrange(info_gain)
    info_gain_cas9cpf1 = rbind(info_gain_cas9, info_gain_cpf1)
    info_gain_cas9cpf1 = info_gain_cas9cpf1[!duplicated(info_gain_cas9cpf1$barcode),]
    info_gain_cas9cpf1$target = factor(info_gain_cas9cpf1$target, levels=c("cpf1", "cas9"))
    info_gain_cas9cpf1$barcode = factor(info_gain_cas9cpf1$barcode, levels=info_gain_cas9cpf1$barcode)
    
    #build a data.frame to store information gain for all targets
    info_gain_alltarget = info_gain_alltarget %>% arrange(info_gain)
    info_gain_alltarget$barcode = factor(info_gain_alltarget$barcode, levels=info_gain_alltarget$barcode)

    return(list(info_gain_all, info_gain_cas9cpf1, father_shannon, info_gain_alltarget))
}
#output: a list, the first data.frame contains the information gain for each target (r1, r2, r3, r4);
#the second dataframe contains IG for each cas9 or cpf1 combination
#the third one is father shannon index before splitting
#the fourth one contains information gain for each target combinations



#--------------------------------------------------------------
#linGrp = 17, should select a specific lineage group
#Note: LinTrans should be a dataframe for a specific sample, while linGrp should be a specific number representing a lineage group
info_gain_all_targets = function(LinTrans, linGrp) {
    #filter those allele and cell type information for a specific lineage group
    lineage_selected = LinTrans %>% filter(lineageGrp==linGrp) %>% 
          mutate(cas9 = paste(intBC, r3, r4, sep="_"), cpf1=paste(intBC, r1, r2, sep="_")) %>% 
          mutate(r1 = paste(intBC, r1, sep="_"), r2 = paste(intBC, r2, sep="_"), r3=paste(intBC, r3, sep="_"), r4=paste(intBC, r4, sep="_")) 

    #calculate father shannon
    target_barcode = lineage_selected %>% select(cellBC, r1, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r1), type = unique(cell_type))
    overall_barcode_cells = as.data.frame(table(target_barcode$type))
    father_shannon = shannon(overall_barcode_cells$Freq[overall_barcode_cells$Freq != 0])

    #calculate target r1
    target_barcode = lineage_selected %>% select(cellBC, r1, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r1), type = unique(cell_type))

    barcodes = unique(grep("[ID]", lineage_selected$r1, value=TRUE))
    info_gain_r1 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_r1 = rbind(info_gain_r1, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_r1) = c('barcode', 'info_gain')

    #calculate target r2
    target_barcode = lineage_selected %>% select(cellBC, r2, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r2), type = unique(cell_type))

    barcodes = unique(grep("[ID]", lineage_selected$r2, value=TRUE))
    info_gain_r2 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_r2 = rbind(info_gain_r2, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_r2) = c('barcode', 'info_gain')

    #calculate target r3
    target_barcode = lineage_selected %>% select(cellBC, r3, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r3), type = unique(cell_type))

    barcodes = unique(grep("[ID]", lineage_selected$r3, value=TRUE))
    info_gain_r3 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_r3 = rbind(info_gain_r3, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_r3) = c('barcode', 'info_gain')

    #calculate target r4
    target_barcode = lineage_selected %>% select(cellBC, r4, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r4), type = unique(cell_type))

    barcodes = unique(grep("[ID]", lineage_selected$r4, value=TRUE))
    info_gain_r4 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_r4 = rbind(info_gain_r4, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_r4) = c('barcode', 'info_gain')

    #calculate target cas9
    target_barcode = lineage_selected %>% select(cellBC, cas9, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(cas9), type = unique(cell_type))

    barcodes = unique(grep("[ID]", lineage_selected$cas9, value=TRUE))
    info_gain_cas9 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_cas9 = rbind(info_gain_cas9, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_cas9) = c('barcode', 'info_gain')

    #calculate target cpf1
    target_barcode = lineage_selected %>% select(cellBC, cpf1, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(cpf1), type = unique(cell_type))

    barcodes = unique(grep("[ID]", lineage_selected$cpf1, value=TRUE))
    info_gain_cpf1 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_cpf1 = rbind(info_gain_cpf1, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_cpf1) = c('barcode', 'info_gain')

    #add target label for separte plot
    if (nrow(info_gain_r1)>0) {
        info_gain_r1$target = "r1"
    }
    if (nrow(info_gain_r2)>0) {
        info_gain_r2$target = "r2"
    }
    if (nrow(info_gain_r3)>0) {
        info_gain_r3$target = "r3"
    }
    if (nrow(info_gain_r4)>0) {
        info_gain_r4$target = "r4"
    }
    if (nrow(info_gain_cpf1)>0) {
        info_gain_cpf1$target = "cpf1"
    }
    if (nrow(info_gain_cas9)>0) {
        info_gain_cas9$target = "cas9"
    }

    #build a datafram to store all infomation gain for four targets
    info_gain_r1 = info_gain_r1 %>% arrange(info_gain)
    info_gain_r2 = info_gain_r2 %>% arrange(info_gain)
    info_gain_r3 = info_gain_r3 %>% arrange(info_gain)
    info_gain_r4 = info_gain_r4 %>% arrange(info_gain)
    
    info_gain_all = rbind(info_gain_r1, info_gain_r2, info_gain_r3, info_gain_r4)
    info_gain_all = info_gain_all[!duplicated(info_gain_all$barcode),]
    info_gain_all$target = factor(info_gain_all$target, levels=c("r1", "r2", "r3", "r4"))
    info_gain_all$barcode = factor(info_gain_all$barcode, levels=info_gain_all$barcode)
    #build a data.frame to store all information gain for either cas9 or cpf1
    info_gain_cas9 = info_gain_cas9 %>% arrange(info_gain)
    info_gain_cpf1 = info_gain_cpf1 %>% arrange(info_gain)
    info_gain_cas9cpf1 = rbind(info_gain_cas9, info_gain_cpf1)
    info_gain_cas9cpf1 = info_gain_cas9cpf1[!duplicated(info_gain_cas9cpf1$barcode),]
    info_gain_cas9cpf1$target = factor(info_gain_cas9cpf1$target, levels=c("cpf1", "cas9"))
    info_gain_cas9cpf1$barcode = factor(info_gain_cas9cpf1$barcode, levels=info_gain_cas9cpf1$barcode)
    
    return(list(info_gain_all, info_gain_cas9cpf1, father_shannon))
}
#output: a list, the first data.frame contains the information gain for each target (r1, r2, r3, r4);
#the second dataframe contains IG for each cas9 or cpf1 combination
#the third one is father shannon index before splitting



#--------------------------------------------------------------
#In order to compare different separation effect of different targets, calculate the total inforamtion gain using all cells irrespective of their lineage group origin
#The following function only calcuate the information gain for those mutated targets, no including non-mutated targets
#Note: LinTrans should be a dataframe for a specific sample, while linGrp should be a specific number representing a lineage group
IG_all__target_all_cells = function(LinTrans, cpf1_mutated=FALSE, cas9_mutated=FALSE) {
    #filter those allele and cell type information for a specific lineage group
    lineage_selected = LinTrans %>% mutate(alltarget=paste("LG", lineageGrp, intBC, r1, r2, r3, r4, sep="_")) %>% 
          mutate(cas9 = paste("LG", lineageGrp, intBC, r3, r4, sep="_"), cpf1=paste("LG", lineageGrp, intBC, r1, r2, sep="_")) %>% 
          mutate(r1 = paste("LG", lineageGrp, intBC, r1, sep="_"), r2 = paste("LG", lineageGrp,intBC, r2, sep="_"), 
          r3=paste("LG", lineageGrp,intBC, r3, sep="_"), r4=paste("LG", lineageGrp,intBC, r4, sep="_")) 
    if (cas9_mutated) {
        lineage_selected = lineage_selected[grep("[ID]",lineage_selected$cas9),]
    }
    #if (cpf1_mutated=="mutated") {
    if (cpf1_mutated) {
        lineage_selected = lineage_selected[grep("[ID]", lineage_selected$cpf1),]
    }

    #calculate father shannon
    target_barcode = lineage_selected %>% select(cellBC, r1, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r1), type = unique(cell_type))
    overall_barcode_cells = as.data.frame(table(target_barcode$type))
    father_shannon = shannon(overall_barcode_cells$Freq[overall_barcode_cells$Freq != 0])

    #calculate target r1, only considering the information_gain in all cells without requiring mutated cpf1 or mutated cas9
    target_barcode = lineage_selected %>% select(cellBC, r1, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r1), type = unique(cell_type))
    #checking the barcode according to single integration barcode
    #only selected those barcodes with mutation to calculate information gain
    barcodes = unique(grep("[ID]", lineage_selected$r1, value=TRUE))
    info_gain_r1 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_r1 = rbind(info_gain_r1, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_r1) = c('barcode', 'info_gain')

    #calculate target r2
    target_barcode = lineage_selected %>% select(cellBC, r2, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r2), type = unique(cell_type))

    barcodes = unique(grep("[ID]", lineage_selected$r2, value=TRUE))
    info_gain_r2 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_r2 = rbind(info_gain_r2, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_r2) = c('barcode', 'info_gain')

    #calculate target r3
    target_barcode = lineage_selected %>% select(cellBC, r3, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r3), type = unique(cell_type))

    barcodes = unique(grep("[ID]", lineage_selected$r3, value=TRUE))
    info_gain_r3 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_r3 = rbind(info_gain_r3, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_r3) = c('barcode', 'info_gain')

    #calculate target r4
    target_barcode = lineage_selected %>% select(cellBC, r4, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(r4), type = unique(cell_type))

    barcodes = unique(grep("[ID]", lineage_selected$r4, value=TRUE))
    info_gain_r4 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_r4 = rbind(info_gain_r4, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_r4) = c('barcode', 'info_gain')

    #calculate target cas9
    target_barcode = lineage_selected %>% select(cellBC, cas9, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(cas9), type = unique(cell_type))

    barcodes = unique(grep("[ID]", lineage_selected$cas9, value=TRUE))
    info_gain_cas9 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_cas9 = rbind(info_gain_cas9, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_cas9) = c('barcode', 'info_gain')

    #calculate target cpf1
    target_barcode = lineage_selected %>% select(cellBC, cpf1, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(cpf1), type = unique(cell_type))

    barcodes = unique(grep("[ID]", lineage_selected$cpf1, value=TRUE))
    info_gain_cpf1 = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_cpf1 = rbind(info_gain_cpf1, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_cpf1) = c('barcode', 'info_gain')

    #calculate all target
    target_barcode = lineage_selected %>% select(cellBC, alltarget, cell_type) %>% group_by(cellBC) %>% 
         summarise(barcodeset = scar_merge(alltarget), type = unique(cell_type))

    barcodes = unique(grep("[ID]", lineage_selected$alltarget, value=TRUE))
    info_gain_alltarget = data.frame('barcode' = character(), 'info_gain'=numeric())
    for (select_barcode in barcodes) {
         info_gain_alltarget = rbind(info_gain_alltarget, info_gain_barcode(target_barcode ,select_barcode))   
    }
    colnames(info_gain_alltarget) = c('barcode', 'info_gain')

    #add target label for separte plot
    if (nrow(info_gain_r1)>0) {
        info_gain_r1$target = "r1"
    }
    if (nrow(info_gain_r2)>0) {
        info_gain_r2$target = "r2"
    }
    if (nrow(info_gain_r3)>0) {
        info_gain_r3$target = "r3"
    }
    if (nrow(info_gain_r4)>0) {
        info_gain_r4$target = "r4"
    }
    if (nrow(info_gain_cpf1)>0) {
        info_gain_cpf1$target = "cpf1"
    }
    if (nrow(info_gain_cas9)>0) {
        info_gain_cas9$target = "cas9"
    }
    if (nrow(info_gain_alltarget)>0) {
        info_gain_alltarget$target = "alltarget"
    }

    #build a datafram to store all infomation gain for four targets
    info_gain_r1 = info_gain_r1 %>% arrange(info_gain)
    info_gain_r2 = info_gain_r2 %>% arrange(info_gain)
    info_gain_r3 = info_gain_r3 %>% arrange(info_gain)
    info_gain_r4 = info_gain_r4 %>% arrange(info_gain)
    
    info_gain_all = rbind(info_gain_r1, info_gain_r2, info_gain_r3, info_gain_r4)
    info_gain_all = info_gain_all[!duplicated(info_gain_all$barcode),]
    info_gain_all$target = factor(info_gain_all$target, levels=c("r1", "r2", "r3", "r4"))
    info_gain_all$barcode = factor(info_gain_all$barcode, levels=info_gain_all$barcode)
    #build a data.frame to store all information gain for either cas9 or cpf1
    info_gain_cas9 = info_gain_cas9 %>% arrange(info_gain)
    info_gain_cpf1 = info_gain_cpf1 %>% arrange(info_gain)
    info_gain_cas9cpf1 = rbind(info_gain_cas9, info_gain_cpf1)
    info_gain_cas9cpf1 = info_gain_cas9cpf1[!duplicated(info_gain_cas9cpf1$barcode),]
    info_gain_cas9cpf1$target = factor(info_gain_cas9cpf1$target, levels=c("cpf1", "cas9"))
    info_gain_cas9cpf1$barcode = factor(info_gain_cas9cpf1$barcode, levels=info_gain_cas9cpf1$barcode)
    
    #build a data.frame to store information gain for all targets
    info_gain_alltarget = info_gain_alltarget %>% arrange(info_gain)
    info_gain_alltarget$barcode = factor(info_gain_alltarget$barcode, levels=info_gain_alltarget$barcode)

    return(list(info_gain_all, info_gain_cas9cpf1, father_shannon, info_gain_alltarget))
}
#output: a list, the first data.frame contains the information gain for each target (r1, r2, r3, r4);
#the second dataframe contains IG for each cas9 or cpf1 combination
#the third one is father shannon index before splitting
#the fourth one contains information gain for each target combinations

#--------------------------------------------------------------
#data structure tranformation
#default order is cpf1 mutated first, cas9 mutated second
lineageD3data = function(EBmerge, sample, lineage,  cl_cols, cpf1_first) {
    lineage$cpf1 = paste(lineage$r1, lineage$r2, sep="|")
    lineage$cas9 = paste(lineage$r3, lineage$r4, sep="|") #be careful, cas9 and cpf1 is correlated by cells.
    LinTrans = left_join(lineage,EBmerge@meta.data, by="cellBC")
    LinTrans = LinTrans %>% select(cellBC, intBC, cpf1, cas9, lineageGrp, cell_type)
    LinTrans = LinTrans[!is.na(LinTrans$cell_type),]
    LinTrans$cpf1 = gsub("[ATCG]", "", LinTrans$cpf1) 
    LinTrans$cas9 = gsub("[ATCG]", "", LinTrans$cas9) 
    LinTrans = LinTrans %>% select(cellBC, lineageGrp, cpf1, cas9, cell_type)
    LinTrans1 = LinTrans %>% group_by(cellBC) %>% 
                                summarise(GeneBarcodeMerge = unique(lineageGrp), Cpf1scarMerge = scar_merge(cpf1), 
                                          Cas9scarMerge = scar_merge(cas9), seurat_cluster = unique(cell_type))


    #org_cells = subset(EBmerge, subset = orig.ident==sample)
    cl_orders = names(cl_cols)

    org_cells = as.data.frame(LinTrans1)
    org_cells$seurat_cluster = factor(org_cells$seurat_cluster)
    org_cells$GeneBarcodeMerge = paste0("linGrp", org_cells$GeneBarcodeMerge)
    
    org_cells <- org_cells[base::order(org_cells$GeneBarcodeMerge, org_cells$Cpf1scarMerge,  org_cells$Cas9scarMerge, setNames(1:length(cl_orders),cl_orders)[as.character(org_cells$seurat_cluster)]),]
    if (cpf1_first==TRUE) {
        org_cells$pathString <- paste(sample,  org_cells$GeneBarcodeMerge, org_cells$Cpf1scarMerge,org_cells$Cas9scarMerge,org_cells$cellBC, sep = "/")
        } else {
        org_cells$pathString <- paste(sample,  org_cells$GeneBarcodeMerge, org_cells$Cas9scarMerge, org_cells$Cpf1scarMerge,org_cells$cellBC, sep = "/")
    }
    
    #org_cells$pathString <- paste(sample,  org_cells$GeneBarcodeMerge, org_cells$Cpf1scarMerge,org_cells$Cas9scarMerge,org_cells$cellBC, sep = "/")
 
    #need to add cellBC to the rownames, as extract colors from lineage_cells_network is neccessary
    rownames(org_cells) = org_cells$cellBC
    return(org_cells)
    }


#define a function to plot lineage trees. Note: already contains file path information
#the following function do not highlight those subclones with mutated targets
#build a directory to save datasets
#EBmerge, single cell meta.data; sample, allele datasets for a specific sample;
#cl_cols, color group to label cell types; cas9, cas9 mutation?; cpf1, cpf1 mutation; 
#cpf1_first, whether cpf1 mutated earlier than cas9? 
#Note: provide input_dir to read files is better
#call_lineage = function(EBmerge, sample, cl_cols, cas9=TRUE, cpf1=TRUE, cpf1_first=TRUE, input_dir, dir_name) {
call_lineage = function(EBmerge, sample, cl_cols, cas9=TRUE, cpf1=TRUE, cpf1_first=TRUE, input_dir, output_dir) {

    #dir_name = "fig5lineagesV3"
    filename1 = paste0(input_dir, sample, 'allCassiopeiaAllele.csv')

    #filename1 = paste0('~/dualproject/lineageData/mannualSharedCellV4/', sample, 'allCassiopeiaAllele.csv')
    lineage = read.csv(filename1)
    #lineage = read.csv(sample_file)
    lineage$cpf1 = paste0(lineage$r1, lineage$r2)
    if (cpf1==TRUE) {
        lineage = lineage[grep("[DI]", lineage$cpf1), ]
    }
    
    lineage$cas9 = paste0(lineage$r3, lineage$r4)
    if (cas9==TRUE) {
        lineage = lineage[grep("[DI]", lineage$cas9) , ]
    }
    
    
    org_cells = lineageD3data(EBmerge, sample, lineage, cl_cols, cpf1_first)
    cat(c('The number of cell is', n_distinct(org_cells$cellBC), ', The number of lineage group is', n_distinct(org_cells$GeneBarcodeMerge), '.\n'))

    lineage_cells <- as.Node(org_cells)
    lineage_cells_network <- ToDataFrameNetwork(lineage_cells, "name")
    lineage_cells_list <- ToListExplicit(lineage_cells, unname = TRUE)
    cols <- c("#303030", cl_cols[as.character(org_cells[lineage_cells_network[,3],"seurat_cluster"])])
    cols[is.na(cols)] <- "#636363"
    jsarray <- paste0('["', paste(cols, collapse = '", "'), '"]')
    nodeStrokeJS <- JS(paste0('function(d, i) { return ', jsarray, '[i]; }'))
    radialNetwork(lineage_cells_list,  nodeColour = nodeStrokeJS, nodeStroke = NA, fontSize = 0)
    BJ333D60indels = radialNetwork(lineage_cells_list,  nodeColour = nodeStrokeJS, nodeStroke = NA, fontSize = 0)

    #create a path to store related lineage groups
    #save_path=paste0("~/dualproject/lineageData/",dir_name, "/")
    #output_dir <- file.path(paste0("~/dualproject/lineageData/",dir_name))
    #cat(c('The dir path is ', save_path))
    #cat(c('The dir path is ', output_dir, '.\n'))
    if (!dir.exists(output_dir)){
        dir.create(output_dir, recursive = T) 
    } 

    if (cas9==TRUE & cpf1==TRUE) {
        #filename2 = paste0('./fig5lineagesV3/', sample, 'SharedClone.html')
        #filename2 = paste0(save_path, sample, 'SharedClone.html')
        filename2 = paste0(output_dir, sample, 'SharedClone.html')
    } else if (cas9==TRUE & cpf1==FALSE) {
         #filename2 = paste0(save_path, sample, 'cas9SharedClone.html')
         filename2 = paste0(output_dir, sample, 'cas9SharedClone.html')
    } else if (cas9==FALSE & cpf1==TRUE) {
         #filename2 = paste0(save_path, sample, 'cpf1SharedClone.html')
         filename2 = paste0(output_dir, sample, 'cpf1SharedClone.html')
    } else {
        #filename2 = paste0(save_path, sample, 'allSharedClone.html')
        filename2 = paste0(output_dir, sample, 'allSharedClone.html')
    }
    
    saveNetwork(BJ333D60indels,file=filename2, selfcontained=T)
}

#optional check the iTracer plot, this itracer plot is related to information gain calculation
#waiting for further improvement
#linIG dataframe from EBlineages_3rd_zhum.ipynb
itracer_plot = function(amplicon_name,linIG,cl_cols, cpf1_first) {
    cl_orders = names(cl_cols)
    #linIG <- linIG[base::order(linIG$grp, linIG$cpf1,  linIG$cas9, setNames(1:length(cl_orders),cl_orders)[as.character(linIG$cell_type)]),]
    if (cpf1_first==TRUE) {
        linIG$pathString <- paste(amplicon_name,  linIG$grp, linIG$cpf1,linIG$cas9,linIG$cellBC, sep = "/")
        } else {
        linIG$pathString <- paste(amplicon_name,  linIG$grp, linIG$cas9, linIG$cpf1,linIG$cellBC, sep = "/")
    }
    linIG = as.data.frame(linIG)
    rownames(linIG) = linIG$cellBC
    lineage_cells <- as.Node(linIG)
    lineage_cells_network <- ToDataFrameNetwork(lineage_cells, "name")
    lineage_cells_list <- ToListExplicit(lineage_cells, unname = TRUE)
    cols <- c("#303030", cl_cols[as.character(linIG[lineage_cells_network[,3],"cell_type"])])
    cols[is.na(cols)] <- "#636363"
    jsarray <- paste0('["', paste(cols, collapse = '", "'), '"]')
    nodeStrokeJS <- JS(paste0('function(d, i) { return ', jsarray, '[i]; }'))
    itracerplot = radialNetwork(lineage_cells_list,  nodeColour = nodeStrokeJS, nodeStroke = NA, fontSize = 0)
    filename = paste0("~/dualproject/lineageData/testFig/",amplicon_name, "_test_clone.html")
    #filename2 = "~/dualproject/lineageData/testFig/D8a_test_clone.html"
    saveNetwork(itracerplot,file=filename, selfcontained=T)
}


#provide a dataframe target_barcode, which contains both selected barcodes and cell types
#return the relative information instead of real information gain, as the sign is more important
#Note: father shannon is distribution index of cell types
info_gain_bc = function(target_barcode, father_shannon) {    
        #this select all the barcodes instead of those barcodes with mutation    
        barcodes = unique(target_barcode$barcode)
        #barcodes = unique(target_barcode[, barcode])
        #print(target_barcode)
        #print(barcodes)
        info_gains = vector()
        for (select_barcode in barcodes) {
                #calculate those shannons with selected barcodes, 
                #careful filter function must refer to the column names instead of a variable
                selected_barcode_cells = target_barcode %>% filter(barcode==select_barcode) %>% 
                        select(cell_type) %>% table() %>% as.data.frame()
                positive = sum(selected_barcode_cells$Freq)
                positive_shannon = shannon(selected_barcode_cells$Freq[selected_barcode_cells$Freq != 0])
                
                #calculate those shannons without selected barcodes
                selected_barcode_cells = target_barcode%>% filter(!(barcode==select_barcode))  %>% 
                        select(cell_type) %>% table() %>% as.data.frame()
                negative = sum(selected_barcode_cells$Freq)
                negative_shannon = shannon(selected_barcode_cells$Freq[selected_barcode_cells$Freq != 0])
                son_shannon = positive_shannon*positive/nrow(target_barcode)+negative_shannon*negative/nrow(target_barcode)
                
                info_gain = round(father_shannon-son_shannon, 5)
                info_gains = append(info_gains, info_gain)
        }
        #print(info_gains)

        #info_gain_alltarget = data.frame(barcodes, info_gains) #no neccessary to create a frame
        #sum_infogain= sum(info_gain_alltarget$info_gain)
        sum_infogain= sum(info_gains)
        if (sum_infogain!=0) {
            relative_infogain = (2^sum_infogain)/(2^father_shannon)     
        } else {
            relative_infogain = 0 #should not contain value if sum_infogain equals to 0
        }
        
        return(relative_infogain)
        #return(c(relative_infogain, sum_infogain, father_shannon))
}

##calculation
#define a function select barcode type and cell types for individual calculation
#the select barcode type includes cas9, cpf1
#calcute a vector of information gain. Currently, each barcode represent the barcodeset of the same level
relative_ig = function(linIG, barcode, selected_cell_type="no_selection") {
    #barcode = "cas9"
    target_barcode = linIG[,c("cellBC", barcode, "cell_type")]
    colnames(target_barcode) = c("cellBC", "barcode", "cell_type")
    #target_barcode %>% filter(grp=="grp19")
    #selected_cell_type = "Hindgut"
    if (selected_cell_type!="no_selection") {
            target_barcode$cell_type = sapply(target_barcode$cell_type, function(x) {
        ifelse(x==selected_cell_type, selected_cell_type, "not")
        })
    }

    #define a function to calculate information gain for dataframe linIG
    overall_barcode_cells = as.data.frame(table(target_barcode$cell_type))
    father_shannon = shannon(overall_barcode_cells$Freq[overall_barcode_cells$Freq != 0])

    #calculate the relative contribution for the information gain of all barcodes
    #relative_ig_bctype= info_gain_bc(target_barcode,barcode, father_shannon)
    #remove barcode parameter function to avoid conflicts with filter function
    relative_ig_bctype= info_gain_bc(target_barcode, father_shannon)
    return(relative_ig_bctype)
    #return(c(relative_ig_bctype, father_shannon))
}

##calculation
#calculate the father shannon
relative_ig_father_shannon = function(linIG, selected_cell_type="no_selection") {
    #barcode = "cas9"
    target_barcode = linIG[,c("cellBC", "cell_type")]
    if (selected_cell_type!="no_selection") {
                target_barcode$cell_type = sapply(target_barcode$cell_type, function(x) {
            ifelse(x==selected_cell_type, selected_cell_type, "not")
            })
        }
    #define a function to calculate information gain for dataframe linIG
    overall_barcode_cells = as.data.frame(table(target_barcode$cell_type))
    father_shannon = shannon(overall_barcode_cells$Freq[overall_barcode_cells$Freq != 0])
    return(father_shannon)
}

##calculation
#define a function to calculate the relative_information gain change for a specific clone or samples
#if a clone, use clone_ig in function, if a sample, use linIG in functions
#create a void dataframe to store relative igs 
#barcodes as columns, cell type as rows
#may remove the calculation of lineageGrps while using clone_ig, waiting for further refinement

clone_ig_dataframe = function(clone_ig, linGrp, cpf1_first) {
    if (cpf1_first) {
        columns = c("grp", "cpf1", "cas9")
    } else {
        columns = c("grp", "cas9", "cpf1")
    }

    rows = unique(clone_ig$cell_type)
    relative_igs = matrix(0, nrow=length(rows), ncol=length(columns))
    relative_igs = as.data.frame(relative_igs)
    colnames(relative_igs) = columns
    rownames(relative_igs) = rows

    #calculate father shannon
    father_shannons = vector()
    for (selected_cell_type in rows) {
            father_shannons = append(father_shannons,relative_ig_father_shannon(clone_ig, selected_cell_type)) 
    }

    for (barcode in columns) {
        for (selected_cell_type in rows) {
            relative_igs[selected_cell_type, barcode] = relative_ig(clone_ig, barcode, selected_cell_type)     
        }
    }
    relative_igs = cbind(relative_igs, father_shannons)
    relative_igs$linGrp = linGrp
    relative_igs$clone_size = nrow(clone_ig)
    relative_igs$cell_type = rownames(relative_igs)
    return(relative_igs)
}


##preprocess
#preprocessing, default, keep only muatated cas9 and cpf1 scars only
#provide a LinTrans dataframe, which is obtained from allele_to_lineage conversion
#should be extremely careful with the naming system, hard to find errors
#lintrans_preprocess = function(LinTrans, cpf1_mutated="mutated", cas9_mutated="mutated") {
lintrans_preprocess = function(LinTrans, cpf1_mutated=TRUE, cas9_mutated=TRUE) {
    LinTrans[LinTrans=="[None]"] = "" 
    LinTrans$cpf1 = paste(LinTrans$r1, LinTrans$r2, sep="") #only for filtering
    LinTrans$cas9 = paste(LinTrans$r3, LinTrans$r4, sep="") #confilcts between parameters and columns
    #if (cas9_mutated=="mutated") {
    if (cas9_mutated) {
        LinTrans = LinTrans[grep("[ID]",LinTrans$cas9),]
    }
    #if (cpf1_mutated=="mutated") {
    if (cpf1_mutated) {
        LinTrans = LinTrans[grep("[ID]",LinTrans$cpf1),]
    }
    
    #LinTrans = LinTrans %>% select(cellBC, intBC, r1, r2, r3, r4, lineageGrp, cell_type)
    LinTrans$r1 = gsub("[", "", LinTrans$r1, fixed=T)
    LinTrans$r1 = gsub("]", "", LinTrans$r1, fixed=T)
    LinTrans$r2 = gsub("[", "", LinTrans$r2, fixed=T)
    LinTrans$r2 = gsub("]", "", LinTrans$r2, fixed=T)
    LinTrans$r3 = gsub("[", "", LinTrans$r3, fixed=T)
    LinTrans$r3 = gsub("]", "", LinTrans$r3, fixed=T)
    LinTrans$r4 = gsub("[", "", LinTrans$r4, fixed=T)
    LinTrans$r4 = gsub("]", "", LinTrans$r4, fixed=T)

    LinTrans$lineageGrp = paste0("grp", LinTrans$lineageGrp)
    LinTrans$cpf1 = paste(LinTrans$r1, LinTrans$r2, sep="") #for real check
    LinTrans$cas9 = paste(LinTrans$r3, LinTrans$r4, sep="") 

    LinTrans = LinTrans %>% select(cellBC, intBC, cas9, cpf1, lineageGrp, cell_type)
    return(LinTrans)
}


##preprocess
##try to calculate the change of information gain for time-sequential barcodes. 
#only care about barset, still keeps the coupling informatin between different targets?
#remove those groups with only one cell types, as they don't contain differentiation information
#actually those groups containing only one cell type may show the preference of clone differentiation biases

#lintrans_to_linig = function(linTrans, cl_cols, cpf1_first=FALSE,remove_single_type="yes") {
lintrans_to_linig = function(linTrans, cl_cols, cpf1_first=FALSE,remove_single_type=TRUE) {
    linIG = linTrans %>% group_by(cellBC) %>% 
            summarise(grp = unique(lineageGrp), cpf1 = scar_ig_merge(cpf1), 
            cas9 = scar_ig_merge(cas9), seurat_cluster = unique(cell_type))

    #link barcode subsets
    #cpf1_first = FALSE #Note: this value is specific to sample D8a
    if (cpf1_first) {
        #linIG$cpf1new = paste0(linIG$grp,linIG$cpf1, sep="")
        #linIG$cas9new = paste0(linIG$grp,linIG$cpf1, linIG$cas9, sep="")
        linIG$cpf1new = paste0(linIG$grp, "|",linIG$cpf1)
        linIG$cas9new = paste0(linIG$cpf1new, sep="l", linIG$cas9)
    } else {
        linIG$cas9new = paste0(linIG$grp, sep="|",linIG$cas9)
        linIG$cpf1new = paste0(linIG$cas9new, sep="l", linIG$cpf1)
    }
    linIG = linIG %>% select(cellBC, grp, cpf1new, cas9new, seurat_cluster)
    colnames(linIG) = c("cellBC", "grp", "cpf1", "cas9","cell_type")

    #to reduce calculation, remove those lineageGrp with only one cell type.
    #if (remove_single_type=="yes") {
    if (remove_single_type) {
        grp_ncelltype = linIG %>% group_by(grp) %>% summarise(n_celltype = n_distinct(cell_type)) 
        grps = grp_ncelltype$grp[grp_ncelltype$n_celltype>1]
        linIG = linIG[linIG$grp %in% grps, ]
    }

    #re-order linIG asendingly
    cl_orders = names(cl_cols)
    linIG <- linIG[base::order(linIG$grp, linIG$cpf1,  linIG$cas9, setNames(1:length(cl_orders),cl_orders)[as.character(linIG$cell_type)]),]
    return(linIG)
}


#this function is similar to iTracer_plot, but deals with linIG without specifying mutations or specifying with mutations
itracer_full_plot = function(amplicon_name,linIG,cl_cols, cpf1_first, output_dir) {
    #preprocessing
    cl_orders = names(cl_cols)
    if (cpf1_first==TRUE) {
        #linIG$cpf1 = paste0(linIG$grp, "|", linIG$cpf1) #add | as seperator
        #linIG$cas9 = paste0(linIG$cpf1, "l", linIG$cas9) #add l as separator
        #linIG <- linIG[base::order(linIG$grp, linIG$cpf1,  linIG$cas9, setNames(1:length(cl_orders),cl_orders)[as.character(linIG$cell_type)]),]
        linIG$pathString <- paste(amplicon_name,  linIG$grp, linIG$cpf1,linIG$cas9,linIG$cellBC, sep = "/")
    } else {
        #linIG$cas9 = paste0(linIG$grp, "|", linIG$cas9)
        #linIG$cpf1 = paste0(linIG$cas9, "l", linIG$cpf1)
        #linIG <- linIG[base::order(linIG$grp, linIG$cpf1,  linIG$cas9, setNames(1:length(cl_orders),cl_orders)[as.character(linIG$cell_type)]),]
        linIG$pathString <- paste(amplicon_name,  linIG$grp, linIG$cas9, linIG$cpf1,linIG$cellBC, sep = "/")
    }
    linIG = as.data.frame(linIG)
    rownames(linIG) = linIG$cellBC

    #create node object for lineageD3 plot
    lineage_cells <- as.Node(linIG)
    lineage_cells_network <- ToDataFrameNetwork(lineage_cells, "name")
    lineage_cells_list <- ToListExplicit(lineage_cells, unname = TRUE)

    #extract those nodes contains both the lineageGrp and scar1 connnection
    idx_node_first_scars <- grep("\\|", lineage_cells_network[,3])
    #select those nodes with mutation
    idx_node_first_scarred_scars <- idx_node_first_scars[sapply(unlist(lapply(strsplit(grep("\\|", lineage_cells_network[,3], value=T),"\\|"), "[", 2)), function(x) grepl("[ID]", x))]
    #those idx less than index_node_scars are lineageGrp nodes
    #idx_node_clones <- seq(min(idx_node_first_scars)-1) 

    #second level
    idx_node_second_scars = grep("l", lineage_cells_network[,3])
    idx_node_second_scarred_scars <- idx_node_second_scars[sapply(unlist(lapply(strsplit(grep("l", lineage_cells_network[,3], value=T),"l"), "[", 2)), function(x) grepl("[ID]", x))]



    #create a list of colors to color the nodes, the first node is sample name
    #both the first node and the cell type node will be color in the following step.
    cols <- c("#303030", cl_cols[as.character(linIG[lineage_cells_network[,3],"cell_type"])])
    cols[is.na(cols)] <- "#636363"

    #change the non-scarred nodes to color #dedede, lightgray, no need to label the upper layer clones
    cols[setdiff(idx_node_first_scars,idx_node_first_scarred_scars)+1] <- "#dedede"
    cols[setdiff(idx_node_second_scars,idx_node_second_scarred_scars)+1] <- "#dedede"
    #cols[idx_node_clones+1] <- "#dedede"

    #JS function to capture the color and names
    jsarray <- paste0('["', paste(cols, collapse = '", "'), '"]')
    nodeStrokeJS <- JS(paste0('function(d, i) { return ', jsarray, '[i]; }'))
    itracerplot = radialNetwork(lineage_cells_list,  nodeColour = nodeStrokeJS, nodeStroke = NA, fontSize = 0)
    #filename = paste0("~/dualproject/lineageData/testFig/",amplicon_name, "_test_clone.html")
    filename = paste0(output_dir,amplicon_name, "_test_clone.html")
    saveNetwork(itracerplot,file=filename, selfcontained=T)
}


#NOTE: This following script block are the main function to caculate the relative information gain
##the whole process for calculating the relative information gain for each barcode level
#calculate relative information gain for linGrp, barcoding first level, barcoding second level
#calculate_relativeig_singlesample = function(amplicon_name, singlecell, cl_cols, cpf1_mutated=TRUE, cas9_mutated=TRUE, cpf1_first=FALSE, remove_single_type=TRUE) {
calculate_relativeig_singlesample = function(amplicon_name, singlecell, cl_cols, cpf1_mutated=TRUE, cas9_mutated=TRUE, cpf1_first=FALSE, remove_single_type=TRUE, input_dir, output_dir) {
    filename1 = paste0(input_dir, amplicon_name, 'allCassiopeiaAllele.csv')
    #sample="D8a"
    #filename1 = paste0('~/dualproject/lineageData/mannualSharedCellV4/', amplicon_name, 'allCassiopeiaAllele.csv')
    allele_table = read.csv(filename1, row.names=1)
    ##preprocessing
    #both lineage and cell type information
    LinTrans = allele_to_lineage(allele_table, singlecell)
    #print(head(LinTrans))
    #preprocessing, default, keep only muatated cas9 and cpf1 scars only
    LinTrans_processed = lintrans_preprocess(LinTrans, cpf1_mutated, cas9_mutated)
    #head(LinTrans)
    #obtatin the final dataframe for information gain calculation
    linIG = lintrans_to_linig(LinTrans_processed, cl_cols, cpf1_first, remove_single_type)
    #print(head(linIG))
    ##optional check the iTracer plot, function defined in common_lineage_function
    if (nrow(linIG)==0) {
        #print("This sample has only one sample cell type after filtering non-cas9/cpf1 mutated alleles, which it will be filter out in the last step")
        #linIG_onetype=lintrans_to_linig(LinTrans_processed, cl_cols, cpf1_first, "no")
        print("iTracer for mutation, keep lineageGrp with only one cell type")
        linIG_onetype=lintrans_to_linig(LinTrans_processed, cl_cols, cpf1_first, FALSE)
        #itracer_full_plot(amplicon_name,linIG_onetype,cl_cols,cpf1_first)
        itracer_full_plot(amplicon_name,linIG_onetype,cl_cols,cpf1_first, output_dir)
    } else {
        itracer_full_plot(amplicon_name, linIG, cl_cols, cpf1_first, output_dir)
    }

    ##calculation
    #check the information gain for the whole samples without considering specific cell types
    #take sample D8a as sample
    columns = c("grp", "cas9", "cpf1")

    relative_igs = vector()
    #calculate relative information gains
    for (barcode in columns) {
        relative_igs = append(relative_igs, relative_ig(linIG, barcode))     
    }

    #calculate father shannon
    father_shannons = relative_ig_father_shannon(linIG)

    #combine the father shannon and relative_igs for further evaluation
    relative_igs = append(relative_igs, father_shannons)
    return(relative_igs)
}




#NOTE: this optional function hasn't been tested
##the whole process for calculating the relative information gain for each barcode level
#read the allele_table
#calculate_relativeig_singlesample=function(amplicon_name, singlecell, cl_cols, cpf1_mutated="mutated", cas9_mutated="mutated", cpf1_first=FALSE, remove_single_type="yes") {
calculate_relativeig_singlesample_optional = function(amplicon_name, singlecell, cl_cols, cpf1_mutated, cas9_mutated, cpf1_first, remove_single_type) {
    #sample="D8a"
    filename1 = paste0('~/dualproject/lineageData/mannualSharedCellV4/', amplicon_name, 'allCassiopeiaAllele.csv')
    allele_table = read.csv(filename1, row.names=1)
    ##preprocessing
    #both lineage and cell type information
    LinTrans = allele_to_lineage(allele_table, singlecell)
    #print(head(LinTrans))
    #preprocessing, default, keep only muatated cas9 and cpf1 scars only
    LinTrans_processed = lintrans_preprocess(LinTrans, cpf1_mutated, cas9_mutated)
    #head(LinTrans)
    #obtatin the final dataframe for information gain calculation
    linIG = lintrans_to_linig(LinTrans_processed, cl_cols, cpf1_first, remove_single_type)
    print(head(linIG))
    ##optional check the iTracer plot, function defined in common_lineage_function
    #if (nrow(linIG)>0 & cpf1_mutated=="mutated" & cas9_mutated=="mutated") {
    if (nrow(linIG)>0 & cpf1_mutated==TRUE & cas9_mutated==TRUE) {
        print("iTracer for mutation, remove lineageGrp with only one cell type")
        itracer_plot(amplicon_name,linIG,cl_cols,cpf1_first)
    } else if (nrow(linIG)==0 & cpf1_mutated==TRUE & cas9_mutated==TRUE) {
        #print("This sample has only one sample cell type after filtering non-cas9/cpf1 mutated alleles, which it will be filter out in the last step")
        #linIG_onetype=lintrans_to_linig(LinTrans_processed, cl_cols, cpf1_first, "no")
        print("iTracer for mutation, keep lineageGrp with only one cell type")
        linIG_onetype=lintrans_to_linig(LinTrans_processed, cl_cols, cpf1_first, FALSE)
        itracer_plot(amplicon_name,linIG_onetype,cl_cols,cpf1_first)
    } else if (nrow(linIG)>0 & cpf1_mutated==FALSE & cas9_mutated==FALSE) {
        print("iTracer for no mutation, remove lineageGrp with only one cell type")
        itracer_full_plot(amplicon_name, linIG, cl_cols, cpf1_first)
    } else {
        print("iTracer for no mutation, keep lineageGrp with only one cell type")
        linIG_onetype=lintrans_to_linig(LinTrans_processed, cl_cols, cpf1_first, FALSE)
        itracer_full_plot(amplicon_name,linIG_onetype,cl_cols,cpf1_first)
    }
    

    ##calculation
    #check the information gain for the whole samples without considering specific cell types
    #take sample D8a as sample
    columns = c("grp", "cas9", "cpf1")

    relative_igs = vector()
    #calculate relative information gains
    for (barcode in columns) {
        relative_igs = append(relative_igs, relative_ig(linIG, barcode))     
    }

    #calculate father shannon
    father_shannons = relative_ig_father_shannon(linIG)

    #combine the father shannon and relative_igs for further evaluation
    relative_igs = append(relative_igs, father_shannons)
    return(relative_igs)
}
