import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from tqdm import tqdm
import networkx as nx
from collections import defaultdict
import pylab

def get_max_depth(G, root):
    """
    Gives the maximum depth of the graph G from some node (typically the root).

    :param G:
        input tree
    :param root:
        Node from which to compute the maximum depth:

    :returns:
        The max depth.
    """

    md = 0

    for n in nx.descendants(G, root):

        if G.nodes[n]["depth"] > md:

            md = G.nodes[n]["depth"]

    return md

def set_progeny_size(G, root):
    """
    Store as a node attribute the number of leaves below a given node. Traverses from the root provided.

    :param G:
        input tree
    :param root:
        Root of tree
    
    :return:
        Graph with progeny size stored as a node attribute (for some node n you can get this value like so:
        G.nodes[n]['prog_size'])
    """

    s = get_progeny_size(G, root)

    G.nodes[root]["prog_size"] = s

    for d in tqdm(G.nodes(), desc="Computing progeny size for each internal node"):

        s = get_progeny_size(G, d)
        G.nodes[d]["prog_size"] = s

    return G

def get_progeny_size(G, node):
    """
    Get the progeny size of a node in the graph. 

    :param G:
        input tree
    :param root:
        internal node.
    
    :return:
        Integer value denoting the number of leaves below that node.

    """

    all_prog = [node for node in nx.dfs_preorder_nodes(G, node)]

    return len([n for n in all_prog if G.out_degree(n) == 0 and G.in_degree(n) == 1])

def get_children_of_clade(G, node):
    """
    Get the names of the leaves below the node in G.

    :param G:
        input tree
    :param root:
        internal node.
    
    :return:
        List of all leaves below the node.

    """

    all_prog = [node for node in nx.dfs_preorder_nodes(G, node)]
    return [n for n in all_prog if G.out_degree(n) == 0 and G.in_degree(n) == 1]

def get_meta_counts(G, node, metavals):
    """
    For all the leaves below the node in G, count how many occurences of each value in metavals.  

    :param G:
        input tree
    :param root:
        internal node.
    :param metavals:
        List of possible metavals.
    
    :return:
        Dictionary mapping each metaval to the number of occurences.
    """

    meta_counts = defaultdict(dict)
    children_vals = [G.nodes[n]["meta"] for n in get_children_of_clade(G, node)]
    for m in metavals:
        meta_counts[m] = children_vals.count(m)

    return meta_counts


def set_depth(G, root):
    """
    Store the depth of each node as an attribute in the graph. 

    :param G:
        input tree
    :param root:
        root of tree.
    
    :return:
        Graph with depth as an attribute for each node. You can access this by using G.nodes[n]['depth'] for any node n.
    """

    depth = nx.shortest_path_length(G, root)

    for d in depth.keys():

        G.nodes[d]["depth"] = depth[d]

    return G

def cut_tree(G, depth):
    """
    Gets the internal nodes at the depth specified in the tree.
    
    :param G:
        input tree
    :param depth:
        Depth at which to cut the tree.

    :return:
        List of internal nodes at the depth specified.
    """

    nodes = []
    leaves = [n for n in G.nodes if G.out_degree(n) == 0]

    for leaf in leaves:
        if G.nodes[leaf]['depth'] <= depth:
            if G.nodes[leaf]['depth'] == 1: 
                nodes.append(leaf)
            else:
                nodes.append(list(G.predecessors(leaf))[0])

    for n in G.nodes:
        if (G.nodes[n]["depth"] == depth) & (n not in leaves):
            nodes.append(n)
    
    return set(nodes)

def assign_meta(G, meta):
    """
    Assign meta items to all leaves in G and store as a node attribute.

    :param G:
        Input graph.
    :param meta:
        pandas Series of meta items, where the index are sample labels.

    :return:
        Graph with meta items assigned to leaves. 
    """


    root = [node for node in G.nodes() if G.in_degree(node) == 0][0]

    leaves = [x for x in G.nodes() if G.out_degree(x) == 0 and G.in_degree(x) == 1]
    metadict = {}

    for l in leaves:
        G.nodes[l]['meta'] = meta.loc[l]

    return G

def add_redundant_leaves(G, cm):
    """
    To fairly take into account sample purity, we'll add back in 'redundant' leaves (i.e.
    leaves that were removed because of non-unique character strings).
    """

    # create lookup value for duplicates
    cm["lookup"] = cm.astype('str').apply(''.join, axis=1)
    net_nodes = np.intersect1d(cm.index, [n for n in G])

    uniq = cm.loc[net_nodes]

    # find all non-unique character states in cm
    nonuniq = np.setdiff1d(cm.index, np.array([n for n in G]))

    for n in nonuniq:

        new_node = str(n)
        _leaf = uniq.index[uniq["lookup"] == cm.loc[n]["lookup"]][0]
        parents = list(G.predecessors(_leaf))
        for p in parents:
            G.add_edge(p, new_node)

        G.nodes[new_node]["depth"] = G.nodes[_leaf]["depth"]

    return G


def nearest_neighbor_dist(G):
    """
    Compute the distance for each leaf to the nearest leaf with the same meta value. 

    :param G:
        Input graph with meta values already mapped.
    
    :return:
        A list consisting of a vector of all nearest neighbor distances and the max distance of the tree to normalize
        by.
    """

    _leaves = [n for n in G if G.out_degree(n) == 0]
    n = len(_leaves)
    tree_dists = np.zeros((n, n))

    g = G.to_undirected()

    for i in tqdm(range(n), desc="computing distances between leaves"):
        l1 = _leaves[i]
        for j in range(i+1, n):
            l2 = _leaves[j]

            if g.nodes[l1]["meta"] == g.nodes[l2]["meta"]:
                tree_dists[i, j] = nx.shortest_path_length(g, l1, l2)

    tree_dists = tree_dists + tree_dists.T - np.diag(np.diag(tree_dists)) # reflect array
    tree_dists[tree_dists == 0] = np.inf # set all 0 values to inf b/c these are not the same meta values

    min_dists = []
    for i in range(n):
        min_dists.append(np.min(tree_dists[i, :]))

    max_dist = nx.diameter(g)

    return min_dists, max_dist

def sample_chisq_test(G, metavals, depth=0):
    """
    Calculates the association between clades and meta variables with a Chi-Squared test every possible depth of the graph. 
    At some given depth, we first  cut the tree (i.e. get the internal nodes at that depth) and compute the number of
    occurences for each meta value under that node. We form a contingency table .. math::(T) of size :math: (C x M) where :math:(C) is the number
    of clades and :math:(M) is the number of meta values. The elements in the table, :math: (m_{i,j}), are the frequencies of meta item 
    :math(m_j) in clade :math:(c_i) . We can then use a Chi-Squared Test to compute the association.

    :param G:
        input tree
    ;param metavals:
        Possible meta values.
    :param depth:
        Depth at which to compute the chisq test.
    
    :return:
        A list consisting of the test statistic, the number of clades at the
        depth.
    """

    nodes = cut_tree(G, depth)

    if len(nodes) == 0:
        return 0, 1

    # metacounts is a list of dictionaries, each tallying the number of
    # occurrences of each meta value in the subclade below the node n
    metacounts = dict(zip(nodes, [get_meta_counts(G, n, metavals) for n in nodes]))

    num_leaves = pylab.sum([G.nodes[m]["prog_size"] for m in metacounts.keys()])

    # make chisq pivot table for test -- M rows (# of subclades) x K cols (# of possible meta items)

    csq_table = np.zeros((len(metacounts.keys()), len(metavals)))

    clade_ids = list(metacounts.keys())
    for i in range(len(clade_ids)):
        k = clade_ids[i]
        clade = metacounts[k]

        for j in range(len(metavals)):
            meta_item = metavals[j]
            csq_table[i, j] = clade[meta_item]

    # drop cols where all 0, an infrequent occurence but can happen when the clades are really unbalanced
    #good_rows = (np.sum(csq_table, axis=1) > 5)
    #csq_table = csq_table[good_rows, :]
    #good_cols = (np.sum(csq_table, axis=0) > 5)
    #csq_table = csq_table[:, good_cols]
    #print(csq_table)

    # screen table before passing it to the test - make sure all variables passed the zero filter
    if np.any(np.sum(csq_table, axis=1) == 0) or np.any(np.sum(csq_table, axis=0) == 0) or len(csq_table) == 0:
        return 0, 0, 1, csq_table.shape[0]

    chisq = stats.chi2_contingency(csq_table)
    tstat, pval = chisq[0], chisq[1]

    n = np.sum(csq_table, axis=None)
    V = np.sqrt(tstat / (n * min(csq_table.shape[0]-1, csq_table.shape[1]-1)))

    return csq_table.shape[0]

def calc_info_gain(G, metavals, depth=0):
    nodes = cut_tree(G, depth)

    leaves = [n for n in G if G.out_degree(n) == 0]

    if len(nodes) == 0:
        return 0

    metacounts = dict(zip(nodes, [get_meta_counts(G, n, metavals) for n in nodes]))

    num_leaves = dict(zip(nodes, [G.nodes[m]["prog_size"] for m in metacounts.keys()]))

    if len(num_leaves) == 1:
        return 0

    entropy = defaultdict(list)

    for k, v in metacounts.items():
        tot = num_leaves[k]
        # print(tot)
        if tot == 1:
            entropy[k].append(0)
            continue
        # print(metacounts[k].values())
        probs = [ s / float(tot) for s in metacounts[k].values()]
        probs = [i for i in probs if i != 0]
        # print(probs)

        entropy[k].append(-np.sum(probs * np.log2(probs)))

    info_gain = np.sum([(num_leaves[k] / float(len(leaves))) * entropy[k][0] for k, v in entropy.items()])
    
    return entropy, info_gain