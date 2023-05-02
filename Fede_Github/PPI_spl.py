import os
import sys
import numpy as np
import networkx as nx
import itertools as it
import random as rd
import umap
import pickle as pk
import pandas as pd



ppi = pd.read_csv("input/ppi_symbol_lcc.csv"delimiter= ',',
           skipinitialspace=True)

G_ppi = nx.from_pandas_edgelist(ppi, 'symbol1', 'symbol2')

G_ppi_lcc = G_ppi.subgraph(max(nx.connected_components(G_ppi), key=len))  # extract lcc graph
print(G_ppi_lcc.number_of_nodes())
print(G_ppi_lcc.number_of_edges())
spl={}

pairwise=it.combinations(G_ppi_lcc.nodes(), 2)
for pair in pairwise:
    spl[pair]= nx.shortest_path_length(G_ppi_lcc, pair[0], pair[1])

with open('input/ppi_spl.pickle', 'wb') as handle:
    pickle.dump(spl, handle, protocol=pickle.HIGHEST_PROTOCOL)
