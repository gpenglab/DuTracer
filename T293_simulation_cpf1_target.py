#!bin/python

#this script only plot the simulated Cas12a/Cpf1 targets
from IPython.display import Image

import numpy as np
import pandas as pd
#import newick
import cassiopeia as cas
import matplotlib.pyplot as plt
import re

#need to write a if else function, avoid the same variable name in python and R when using r-reticulate
allele_table_cpf1 = pd.read_csv('/data1/home/gdpeng/chengchen/dualproject/sample293Tv2/T293_cpf1_target_simulation_allele_table.csv',
                           usecols = ['cellBC', 'intBC', 'r1', 'r2', 'allele', 'lineageGrp', 'readCount', 'UMI'])
                           

# #function may be required for next step
# def get_default_cut_site_columns(allele_table_cpf1):
#     """Retrieves the default cut-sites columns of an AlleleTable.
#     A basic helper function that will retrieve the cut-sites from an AlleleTable
#     if the AlleleTable was created using the Cassiopeia pipeline. In this case,
#     each cut-site is denoted by an integer preceded by the character "r", for
#     example "r1" or "r2".
#     Args:
#         allele_table: AlleleTable
#     Return:
#         Columns in the AlleleTable corresponding to the cut sites.
#     """
#     cut_sites = [
#         column
#         for column in allele_table_cpf1.columns
#         if bool(re.search(r"r\d", column))
#     ]

#     return cut_sites

#estimate indel priors
#indel_priors = cas.pp.compute_empirical_indel_priors(allele_table_cpf1, grouping_variables=['intBC', 'lineageGrp'])
indel_priors = pd.read_csv("/data1/home/gdpeng/chengchen/dualproject/sample293Tv2/sc293T_indel_priors.csv", index_col=0)

character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(allele_table_cpf1,
                                                                                         allele_rep_thresh = 0.9,
                                                                                         mutation_priors = indel_priors)

#build the tree structure
cas_tree = cas.data.CassiopeiaTree(character_matrix=character_matrix, priors=priors)

# create a basic vanilla greedy solver
vanilla_greedy = cas.solver.VanillaGreedySolver()

# reconstruct the tree
vanilla_greedy.solve(cas_tree, collapse_mutationless_edges=True)

#save newick files
f = open("/data1/home/gdpeng/chengchen/dualproject/sample293Tv2/newick_trees/tree_cpf1target.newick", 'w')
f.write(cas_tree.get_newick())
f.close()

#save tree figures
figure = cas.pl.plot_matplotlib(cas_tree, add_root=True)[0]
figure.savefig("/data1/home/gdpeng/chengchen/dualproject/sample293Tv2/figures/tree_cpf1target.pdf")
