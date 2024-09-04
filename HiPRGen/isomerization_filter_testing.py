# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 08:54:47 2024

@author: jacob
"""

import copy
import networkx as nx
from networkx.exception import NetworkXNoCycle
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
import glob
from monty.serialization import dumpfn
import os
import pickle
import sys

pickle_directory = \
    "C:/Users/jacob/Kinetiscope_Utils/name_molecules"
    
os.chdir(pickle_directory)

#Note: the pickle file needs to be in the same folder as HiPRGen for the 
#mol_entry objects to be loaded

print('Loading mol_entries.pickle...')

with open('mol_entries.pickle', 'rb') as f: #loads HiPRGen mol_entry objs
        
    mol_entries = pickle.load(f)
    
print('Done!')

for species in mol_entries:
   species_ugraph = nx.Graph(species.mol_graph.graph)
   cycle_basis = nx.cycle_basis(species_ugraph)
   if cycle_basis:
       print(cycle_basis)
   # try:
   #     print(nx.find_cycle(species_ugraph))
   # except NetworkXNoCycle:
   #     pass