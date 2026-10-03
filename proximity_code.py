#!/usr/bin/env python
# coding: utf-8

# In[12]:


import csv,os,time
import random
from collections import defaultdict
import networkx, random
import numpy
import pandas as pd
import matplotlib.pyplot as plt

def get_network(network_file, only_lcc):
    network = create_network_from_sif_file(network_file, use_edge_data = False, delim =',', include_unconnected=True)
    #print len(network.nodes()), len(network.edges())
    if only_lcc and not network_file.endswith(".lcc"):
        print ("Shrinking network to its LCC", len(network.nodes()), len(network.edges()))
        components = get_connected_components(network, False)
        network = get_subgraph(network, components[0])
        print ("Final shape:", len(network.nodes()), len(network.edges()))
        #print len(network.nodes()), len(network.edges())
        network_lcc_file = network_file + ".lcc"
        if not os.path.exists(network_lcc_file ):
            f = open(network_lcc_file, 'w')
            for u,v in network.edges():
                f.write("%s 1 %s\n" % (u, v))
            f.close()
    return network

# Network creation from connected network file
def create_network_from_sif_file(network_file_in_sif, use_edge_data = False, delim =',', include_unconnected=True):
    setNode, setEdge, dictDummy, dictEdge = get_nodes_and_edges_from_sif_file(network_file_in_sif, store_edge_type = use_edge_data, delim = delim)
    g = create_graph()
    if include_unconnected:
        g.add_nodes_from(setNode)
    if use_edge_data:
        for e,w in dictEdge.items():
            u,v = e
            g.add_edge(u,v,w=w) #,{'w':w})
    else:
        g.add_edges_from(setEdge)
    return g

def get_nodes_and_edges_from_sif_file(file_name, store_edge_type = False, delim=',', data_to_float=True):
    """
    Parse sif file into node and edge sets and dictionaries
    returns setNode, setEdge, dictNode, dictEdge
    store_edge_type: if True, dictEdge[(u,v)] = edge_value
    delim: delimiter between elements in sif file, if None all whitespaces between letters are considered as delim
    """
    setNode = set()
    setEdge = set()
    dictNode = {}
    dictEdge = {}
    flag = False
    f=open(file_name)
    for line in f:
        if delim == ',':
            words = line.rstrip("\n").split(',')
        else:
            words = line.rstrip("\n").split(delim)
        id1 = words[0]
        id2 = words[1]
        
        setNode.add(id1)
        setNode.add(id2)
        setEdge.add((id1, id2))
        
    f.close()
    if len(setEdge) == 0:
        setEdge = None
    if len(dictNode) == 0:
        dictNode = None
    if len(dictEdge) == 0:
        dictEdge = None
    if flag:
        print( "Warning: Ignored extra columns in the file!")
    return setNode, setEdge, dictNode, dictEdge


def create_graph(directed=False):
    """
        Creates & returns a graph
    """
    if directed:
        g = networkx.DiGraph()
    else:
        g = networkx.Graph()
    return g

def get_connected_components(G, return_as_graph_list=True):
    """
        Finds (strongly in the case of directed network) connected components of graph
        returnAsGraphList: returns list of graph objects corresponding to connected components (from larger to smaller)
        otherwise returns list of node list corresponding nodes in connected components
    """
    result_list = []

    if return_as_graph_list:
        result_list = networkx.connected_component_subgraphs(G)
    else:
        result_list = [c for c in sorted(networkx.connected_components(G), key=len, reverse=True)]

    return result_list

def get_subgraph(G, nodes):
    """
    NetworkX subgraph method wrapper
    """
    
    return G.subgraph(nodes)

def calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=None, n_random=1000, min_bin_size=2, seed=452456, lengths=None, out_file="output1.txt"):
    """
    Calculate proximity from nodes_from to nodes_to
    If degree binning or random nodes are not given, they are generated
    lengths: precalculated shortest path length dictionary
    """
    #distance = "closest"
    #lengths = network_utilities.get_shortest_path_lengths(network, "../data/toy.sif.pcl")
    #d = network_utilities.get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
    nodes_network = set(network.nodes())
    nodes_from = set(nodes_from) & nodes_network 
    nodes_to = set(nodes_to) & nodes_network
    if len(nodes_from) == 0 or len(nodes_to) == 0:
        return None # At least one of the node group not in network
    d = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
    print(d)
    if bins is None and (nodes_from_random is None or nodes_to_random is None):
        bins = get_degree_binning(network, min_bin_size, lengths) # if lengths is given, it will only use those nodes
    if nodes_from_random is None:
        nodes_from_random = get_random_nodes(nodes_from, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    if nodes_to_random is None:
        nodes_to_random = get_random_nodes(nodes_to, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
    random_values_list = zip(nodes_from_random, nodes_to_random)
    values = numpy.empty(len(nodes_from_random)) #n_random
    for i, values_random in enumerate(random_values_list):
        nodes_from, nodes_to = values_random
        #values[i] = network_utilities.get_separation(network, lengths, nodes_from, nodes_to, distance, parameters = {})
        values[i] = calculate_closest_distance(network, nodes_from, nodes_to, lengths)
        
    #pval = float(sum(values <= d)) / len(values) # needs high number of n_random
    
    m, s = numpy.mean(values), numpy.std(values)
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
        print(m)
        print(s)
    return d, z, (m, s) #(z, pval)

def calculate_closest_distance(network, nodes_from, nodes_to, lengths=None):
    values_outer = []
    if lengths is None:
        for node_from in nodes_from:
            values = []
            for node_to in nodes_to:
                val = get_shortest_path_length_between(network, node_from, node_to)
                values.append(val)
            d = min(values)
            #print d,
            values_outer.append(d)
    else:
        for node_from in nodes_from:
            values = []
            vals = lengths[node_from]
            for node_to in nodes_to:
                val = vals[node_to]
                values.append(val)
            d = min(values)
            values_outer.append(d)
    d = numpy.mean(values_outer)
    #print d
    return d

def get_shortest_path_length_between(G, source_id, target_id):
    return networkx.shortest_path_length(G, source_id, target_id)
def get_degree_binning(g, bin_size, lengths=None):
    degree_to_nodes = {}
    for node, degree in g.degree(): #.iteritems(): # iterator in networkx 2.0
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    value = degree_to_nodes.keys()
    values= list(sorted(value))
    print(type(values))
    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1 
        #print i, low, high, len(val) 
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins

def get_random_nodes(nodes, network, bins=None, n_random=1000, min_bin_size=2, degree_aware=True, seed=None):
    if bins is None:
        # Get degree bins of the network
        bins = get_degree_binning(network, min_bin_size) 
    nodes_random = pick_random_nodes_matching_selected(network, bins, nodes, n_random, degree_aware, seed=seed) 
    return nodes_random

def pick_random_nodes_matching_selected(network, bins, nodes_selected, n_random, degree_aware=True, connected=False, seed=None):
    """
    Use get_degree_binning to get bins
    """
    if seed is not None:
        random.seed(seed)
    values = []
    nodes = network.nodes()
    for i in range(n_random):
        if degree_aware:
            if connected:
                raise ValueError("Not implemented!")
            nodes_random = set()
            node_to_equivalent_nodes = get_degree_equivalents(nodes_selected, bins, network)
            for node, equivalent_nodes in node_to_equivalent_nodes.items():
                #nodes_random.append(random.choice(equivalent_nodes))
                chosen = random.choice(equivalent_nodes)
                for k in range(20): # Try to find a distinct node (at most 20 times)
                    if chosen in nodes_random:
                        chosen = random.choice(equivalent_nodes)
                nodes_random.add(chosen)
            nodes_random = list(nodes_random)
        else:
            if connected:
                nodes_random = [ random.choice(nodes) ]
                k = 1
                while True:
                    if k == len(nodes_selected):
                        break
                    node_random = random.choice(nodes_random)
                    node_selected = random.choice(network.neighbors(node_random))
                    if node_selected in nodes_random:
                        continue
                    nodes_random.append(node_selected)
                    k += 1
            else:
                nodes_random = random.sample(nodes, len(nodes_selected))
        values.append(nodes_random)
    return values

def get_degree_equivalents(seeds, bins, g):
    seed_to_nodes = {}
    for seed in seeds:
        d = g.degree(seed)
        for l, h, nodes in bins:
            if l <= d and h >= d:
                mod_nodes = list(nodes)
                mod_nodes.remove(seed)
                seed_to_nodes[seed] = mod_nodes
                break
    return seed_to_nodes



def calculate_proximity_multiple(network, from_file=None, to_file=None, n_random=1000, min_bin_size=2, seed=452456, lengths=None, out_file="output2.txt"):
    """
    Run proximity on each entries of from and to files in a pairwise manner
    output is saved in out_file (e.g., output.txt)
    """
    nodes = set(network.nodes())
    pathwaygenes= get_pathway_genes(from_file, nodes = nodes)
    #drug_to_targets = dict((drug, nodes & targets) for drug, targets in drug_to_targets.iteritems())
    drug_target_genes= get_drug_target_genes(to_file, nodes = nodes)
    # Calculate proximity values
    print (len(drug_target_genes), len(pathwaygenes))
    # Get degree binning 
    bins = get_degree_binning(network, min_bin_size)
    f = open(out_file, 'w')
    for drug, nodes_from in drug_target_genes.items():
        for pathways, nodes_to in pathwaygenes.items():
            print (pathways,drug)
            d, z, (m, s) = calculate_proximity(network, nodes_from, nodes_to, nodes_from_random=None, nodes_to_random=None, bins=bins, n_random=n_random, min_bin_size=min_bin_size, seed=seed, lengths=lengths)
            f.write("%s\t%s\t%f\n" % (pathways,drug,z))
        
    f.close()
    return 

def get_drug_target_genes(from_file, nodes=None, network=None):
    """
    If nodes is not None, keep only nodes in the network
    If network is not None, keep only LCC
    """
    drug_t0_genes = {}
    for line in open(from_file):
        words = line.strip("\n").split("\t")
        drug = words[0].strip('"')
        genes =set(words[1:])
        if nodes is not None:
            genes &= nodes
            if len(genes) == 0:
                continue
        if network is not None:
            network_sub = network.subgraph(genes)
            genes =get_connected_components(network_sub, False)[0]
        drug_t0_genes[drug] = genes
    return drug_t0_genes

def get_pathway_genes(pathway_file, nodes=None, network=None):
    """
    If nodes is not None, keep only nodes in the network
    If network is not None, keep only LCC
    """
    pathway_t0_genes = {}
    for line in open(pathway_file):
        line = line.strip("\n").split("\t")
        reactome = line[0].strip('"')
        genes = set(line[1:])
        if nodes is not None:
            genes &= nodes
            if len(genes) == 0:
                continue
        if network is not None:
            network_sub = network.subgraph(genes)
            genes =get_connected_components(network_sub, False)[0]      
        pathway_t0_genes[reactome]= genes
        
    return pathway_t0_genes
print( 'starting time:, ', time.ctime())
file_name="toy.csv"
from_file="from.txt"
to_file="to.txt"
network = get_network(file_name, only_lcc = True)
#calculate_proximity(network, from_file, to_file, n_random=1000, min_bin_size=2, seed=452456, lengths=None, out_file="output_p.txt")
calculate_proximity_multiple(network, from_file, to_file, n_random=1000, min_bin_size=2, seed=452456, lengths=None, out_file="output_Gemcitabine_pm.txt")
print( 'Ending time:, ', time.ctime())


# In[ ]:





# In[ ]:





# In[ ]:




