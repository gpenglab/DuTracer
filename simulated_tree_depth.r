#This script provide basic functions to calculate the simulated tree depth
library(dplyr)
library(reticulate) #for running python scripts
library(ggtree)
library(treeio)
library(tidytree)
library(tidyverse)
library(ggsignif)
library(ggsci)
library(scales)
library(ggplot2)
library(ungeviz)
#define functions to simulation table
simulation_table = function(allele_table) {
#simulation,sample haft intBCs for tree comparison
    #set.seed(123)
    #select half integration barcodes that contains all four targets
    #we don't need to care about sampling with replacement or without replacement, since this simulation doesn't care about probilities.
    #this simulation focus on selection only.
    alltarget_table = allele_table %>% group_by(cellBC) %>% summarise(intBC = sample(intBC, (n_distinct(intBC)-n_distinct(intBC)%%2)/2))
    alltarget_table = left_join(alltarget_table, allele_table, by=c("cellBC", "intBC"))
    alltarget_table = as.data.frame(alltarget_table)

    #print(dim(alltarget_table))

    #simulation,sample even number of intBCs  tree comparison, only select Cas9 target or Cpf1 target
    #set.seed(123)
    #select all the integration barcodes, then only selected cas9 targets or cpf1 targets
    half_target_table = allele_table %>% group_by(cellBC) %>% summarise(intBC = sample(intBC, (n_distinct(intBC)-n_distinct(intBC)%%2)))
    half_target_allele_table = left_join(half_target_table, allele_table, by=c("cellBC", "intBC"))
    half_target_allele_table = as.data.frame(half_target_allele_table)

    #print(dim(half_target_allele_table))

    cpf1_target_table = half_target_allele_table %>% select(cellBC, intBC, allele, r1, r2, lineageGrp, UMI, readCount)
    cas9_target_table = half_target_allele_table %>% select(cellBC, intBC, allele, r3, r4, lineageGrp, UMI, readCount)
    tables = list(alltarget_table, cpf1_target_table, cas9_target_table)
    return(tables)
}

#define a function to calculate mean tree deptn for the same number of targets for single/dual nuclease
#please check the paper for detailed simulation
call_mean_depth = function(allele_table,repetition=100) {
    alltarget_mean = vector()
    cpf1_target_mean = vector()
    cas9_target_mean = vector()
    cell_number = n_distinct(allele_table$cellBC)
    for (i in 1:repetition) {
        simulation_target_table = simulation_table(allele_table)
        alltarget_table = simulation_target_table[[1]]
        cpf1_target_table = simulation_target_table[[2]]
        cas9_target_table = simulation_target_table[[3]]
        write.csv(alltarget_table, "~/dualproject/sample293Tv2/T293_all_target_simulation_allele_table.csv")
        write.csv(cpf1_target_table, "~/dualproject/sample293Tv2/T293_cpf1_target_simulation_allele_table.csv")
        write.csv(cas9_target_table, "~/dualproject/sample293Tv2/T293_cas9_target_simulation_allele_table.csv")

        #test the running of python script to obtain the tree using cassiopeia packages
        #BE CAREFUL, source_python function will absorb the variables and functions from python, making R variable  strange
        #source_python("./T293_simulation_all_target.py")
        py_run_file("~/dualproject/sample293Tv2/T293_simulation_all_target.py")
        tree<-read.newick("~/dualproject/sample293Tv2/newick_trees/tree_all_target.newick",node.label = "label")
        #calculate the mean tree depth, the tree is reconstructed using Cas9 targets only
        tree_depth = castor::get_all_distances_to_root(tree)
        #summary(tree_depth[1:cell_number])
        #mean_all_tree_depth = mean(tree_depth[1:cell_number]) #mean value is integer and avoid outliers
        #mean seem to be the major cause for the apperance of biomodal distribution, change to mean
        mean_all_tree_depth = mean(tree_depth[1:cell_number])
        alltarget_mean = append(alltarget_mean, mean_all_tree_depth)

        #test the running of python script to obtain the tree using cassiopeia packages
        #source_python("./T293_simulation_cpf1_target.py")
        py_run_file("~/dualproject/sample293Tv2/T293_simulation_cpf1_target.py")
        #tree<-read.newick("~/dualproject/sample293Tv2/newick_trees/tree_cas9target.newick",node.label = "label") #one error here
        tree<-read.newick("~/dualproject/sample293Tv2/newick_trees/tree_cpf1target.newick",node.label = "label")
        #calculate the mean tree depth, the tree is reconstructed using Cas9 targets only
        tree_depth = castor::get_all_distances_to_root(tree)
        #summary(tree_depth[1:cell_number])
        mean_cpf1_tree_depth = mean(tree_depth[1:cell_number]) #mean value is integer and avoid outliers, but may induce multi-modal distribution
        #mean_cpf1_tree_depth = mean(tree_depth[1:cell_number]) #change to using "mean" instead
        cpf1_target_mean = append(cpf1_target_mean, mean_cpf1_tree_depth)

        #test the running of python script to obtain the tree using cassiopeia packages
        #source_python("./T293_simulation_cas9_target.py")
        py_run_file("~/dualproject/sample293Tv2/T293_simulation_cas9_target.py")
        tree<-read.newick("~/dualproject/sample293Tv2/newick_trees/tree_cas9target.newick",node.label = "label")
        #calculate the mean tree depth, the tree is reconstructed using Cas9 targets only
        tree_depth = castor::get_all_distances_to_root(tree)
        #summary(tree_depth[1:cell_number])
        mean_cas9_tree_depth = mean(tree_depth[1:cell_number]) #mean value is integer and avoid outliers 
        cas9_target_mean = append(cas9_target_mean, mean_cas9_tree_depth)
        }
    simulation_mean = data.frame(alltarget_mean,  cpf1_target_mean, cas9_target_mean)
    return(simulation_mean)
           
}

options(dplyr.summarise.inform = FALSE)
mean_tree_depth = call_mean_depth(allele_table, repetition=100)
tree_depth = mean_tree_depth %>% pivot_longer(cols=1:3, names_to = "target", values_to = "tree_depth_mean") 
