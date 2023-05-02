import numpy as np
import networkx as nx
import itertools as it
import random as rd
import pickle as pk
import pandas as pd


def random_selection(G,size):
    while True:
        rnd_sample=rd.sample(list(G.nodes()),size)
        if 'nan' in rnd_sample:
            rnd_sample=rd.sample(list(G.nodes()),size)
        else:
            break
    return rnd_sample

def calculate_closest_distance(G_ppi_lcc,spl, nodes_from, nodes_to):
    values_outer = []
    for node_from in nodes_from:
        values = []
        for node_to in nodes_to:
            if node_from==node_to:
                val =0
            else:
                try:
                    try:
                        val = spl[node_from,node_to]
                    except:
                        val = spl[node_to,node_from]
                except:
                    val=len(nx.shortest_path(G_ppi_lcc,source=node_from, target=node_to))
            values.append(val)
        d = min(values)
        #print d,
        values_outer.append(d)
    d = np.mean(values_outer)
    #print d
    return d

def calculate_proximity(network,spl, nodes_from, nodes_to, exp_random):   #let's calculate the proximity between the two sets
    d=calculate_closest_distance(network,spl,nodes_from,nodes_to)         #exposure and disease genesets
    m=exp_random[len(nodes_from),len(nodes_to)][0]
    s=exp_random[len(nodes_from),len(nodes_to)][1]
    if s == 0:
        z = 0.0
    else:
        z = (d - m) / s
    return d, z, (m, s) #(z, pval)

#Here, we define the PPI
ppi = pd.read_csv("input/ppi_symbol_lcc.csv",delimiter= ',',
           skipinitialspace=True)
G_ppi = nx.from_pandas_edgelist(ppi, 'symbol1', 'symbol2')
print("The total PPI is read and it is big: %s" %(len(G_ppi)))

#Here, we define the set of degs for each condition

with open('input/heart_degs_ppi_and_core_components_dict.pickle', 'rb') as handle:
    heart_degs_ppi_and_core_components_dict = pk.load(handle)

#Here we define the set of diseases

with open('input/disgenet_gda_curated_score.pickle', 'rb') as handle:
    disgenet_gda_curated_score = pk.load(handle)
print("This is the number of all diseases from disgenet: %s" %(len(disgenet_gda_curated_score)))

#Let's consider only diseases with at least 3 associated genes
disgenet_gda_curated_score_multigenes={}
for dis,geneset in disgenet_gda_curated_score.items():
    new_geneset=set(geneset) & G_ppi.nodes()
    if len(new_geneset)>2:
        disgenet_gda_curated_score_multigenes[dis]=new_geneset

#Let's import the pre-calculated spl dictionary
with open('input/ppi_spl.pickle', 'rb') as handle: #This file is too large for Github, but it can be created by running PPI_spl.py
    spl = pk.load(handle)

disease_geneset_distribution=[]
for v in disgenet_gda_curated_score_multigenes.values():
    disease_geneset_distribution.append(len(v))

heart_degs_ppi_and_core_components_proximity={}
for condition,genelist in heart_degs_ppi_and_core_components_dict.items():
    bpa_spec_size=len(genelist)
    bpa_random_distance={}

    for dis_size in set(disease_geneset_distribution):
        S = 10000
        l_rnd_dist = []
        for s in range(S):
            rnd_exp_set=random_selection(G_ppi,bpa_spec_size)
            rnd_dis_set=random_selection(G_ppi,dis_size)
            rnd_dist=calculate_closest_distance(G_ppi,spl,rnd_exp_set,rnd_dis_set)
            l_rnd_dist.append(rnd_dist)
        mu = np.mean(l_rnd_dist)
        std = np.std(l_rnd_dist)

        bpa_random_distance[bpa_spec_size,dis_size]=[mu,std]

    for dis in disgenet_gda_curated_score_multigenes.keys():
        heart_degs_ppi_and_core_components_proximity[condition,dis]=calculate_proximity(G_ppi,spl,genelist,disgenet_gda_curated_score_multigenes[dis],bpa_random_distance)


#Let's save it
with open('output/heart_degs_ppi_and_core_components_proximity.pickle', 'wb') as handle:
    pk.dump(heart_degs_ppi_and_core_components_proximity, handle, protocol=pk.HIGHEST_PROTOCOL)
