# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:26:19 2023

@author: JRMilton
"""

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
import glob
import pickle
from monty.serialization import loadfn, dumpfn
import copy
import networkx as nx
import time

#want to write the json we've created as a "key" showing names associated with
#molecule graphs

#want to write to an excel file reactions of the form: name1 + name2 -> name3+ name4

def find_molecules_from_mpculeids(participants, mol_entries):
    """
    Takes a list of lists of reactant and product mpculeids and converts them to
    pymatgen Molecule objects

    participants: a dictionary of the form: {"reatants": "reactantmpculeid", "products": etc.}
    mol_entries: a pickled representation of an object from the open file object file with the reconstituted object
    hierarchy--contains pymatgen Molecule objects representing species in our reaction network which are associated
    with a given mpculeid

    Returns: new_dict, a new dictionary with the same keys but a list of Molecule objects as the values
    """
    new_dict = {}
    for key, side in participants.items():
        molecule_list = []
        for mpculeid in side:
            for mol_entry in mol_entries:
                m_id = mol_entry.entry_id
                if m_id == mpculeid:
                    molecule_list.append(mol_entry)
        new_dict[key] = molecule_list
    return new_dict
    
def name_molecule(mol, func_group_dict):
    """
    Takes a pymatgen Molecule object and generates a "name" for that object based
    on the functional groups present. Any atoms that are not associated with a 
    functional group will be be added to the name; this also means that species
    with no functional groups present will simply return their composition.
    
    Parameters
    ----------
    mol : pymatgen Molecule object
        The species we wish to name
    func_group_dict : dictionary
        A dictionary whos keys are the names of functional groups and values are
        the Molecule objects associated with those groups

    Returns
    -------
    name : string
        The generated "name" of the species

    """
    mol_graph = nx.Graph(mol.graph)
    mol_copy = copy.deepcopy(mol_graph)
    name = ""
    for n, group in func_group_dict.items(): 
        group_undirected_graph = nx.Graph(MoleculeGraph.with_local_env_strategy(group, OpenBabelNN(order = False)).graph) #Creates a Networkx (undirected) Graph object from func Molecule object
        while functional_group_present(mol_copy, group_undirected_graph)[0]: #want to repeat until func group is no longer present
            if n not in name:
                name = name + n + '_' #if present, add func group name to name
            else: #if the functional group is already in the name, we just increase the number present
                if name[-2].isnumeric():
                    new_number = int(name[-2]) + 1
                    name = name[0:len(name)-2]+ str(new_number) + '_'
                else:
                    name = name[0:len(name)-1] + '2_'
            mappings = list(functional_group_present(mol_copy, group_undirected_graph)[1]) #points to location of functional group in molecule
            for atom_index in mappings[0].keys():
                mol_copy.remove_node(atom_index)  #and remove the group
    if name and mol_copy: #i.e. if we have added some functional groups to the name but there remain atoms not in functional group
        comp = ""
        for node in mol_copy.nodes(data = True): #each node is a 2-tuple containing the index of that node and a dictionary with data related to that node
            species = node[1]['specie']
            if species in comp:
                species_index = comp.index(species)
                new_number = int(comp[species_index+1]) + 1
                comp = comp[0:len(comp)-1]+str(new_number)
            else:
                comp = comp + species + "1"
        name = name + comp + '_'
    elif name and not mol_copy: #i.e. our functional groups encompass all atoms in the molecule
        pass
    else: #i.e. the molecule contains no functional groups
        name = str(mol.molecule.composition).replace(" ", "")
        name = name + '_'
    if mol.molecule.charge == 1: #add charges to the end of the names
        name = name + "+" + str(mol.molecule.charge)
    else:
        name = name + str(mol.molecule.charge)
    return name

def functional_group_present(mol_graph, func):
    """
    Tests whether or not a given functional group is present in a molecule via
    testing if the Networkx graphical representation of a molecule contains a 
    subgraph that is isomorphic to the functional group.
    
    Parameters
    ----------
    mol_graph : Networkx Undirected Graph
        Graph of the molecule you want to name
    func : Networkx Undirected Graph
        Graph of the functional group

    Returns
    -------
    Boolean
        True if a subgraph of the molecule is isomorphic to the functional group
        graph, False otherwise.
    Generator over isomorphisms between a subgraph of G1 and G2.
        This is a generator of mappings between a subgraph of the molecule and
        the functional group--this helps us deterimine where the functional
        group is in the original molecule

    """
    nm = nx.isomorphism.categorical_node_match("specie", None) #ensures isomorphic graphs must have the same atoms
    isomorphism = nx.isomorphism.GraphMatcher(mol_graph, func, node_match = nm)
    return isomorphism.subgraph_is_isomorphic(), isomorphism.subgraph_isomorphisms_iter()

def stereoisomer_test(test_dict, mol_entry_id, mol_name):
    """
    Structural isomers--species with the same formula but different graphs--will return
    the same name under this paradigm. This function numbers each isomer such that the
    names remain unique

    Parameters
    ----------
    test_dict : dictionary
        dictionary whose keys are names and values are entry_ids of the species
        associated with that name
    mol_entry_id : string
        the unique entry id associated with a species
    mol_name : string
        the name generated for a species from the name_molecule function, which
        may not yet be unique

    Returns
    -------
    stereos : dictionary
        a dictionary of length 1 or 2, whose keys are now unique names of stereoisomers
        and values are entryids

    """
    stereos = {}
    if mol_name in test_dict and mol_entry_id != test_dict[mol_name]: #i.e. if name already in dictionary but graphs of molecules are different
        old_id = test_dict[mol_name]
        name_1 = mol_name + '_#1'
        name_2 = mol_name + '_#2'
        stereos.update({name_1: old_id, name_2: mol_entry_id})
    elif mol_name + '_#1' in test_dict and mol_entry_id != test_dict[mol_name + '_#1']: #i.e. if multiple stereoisomers are already in the dictionary
        current_max = 3
        while mol_name + '_#' + str(current_max) in test_dict and mol_entry_id != test_dict[mol_name + '_#' + str(current_max)]:
            current_max += 1
        new_name = mol_name + '_#' + str(current_max)
        stereos[new_name] = mol_entry_id
    return stereos

start = time.time()
print('Associating functional groups with their Molecule objects...')
func_group_dict = {} #take a group of xyz files associated with our functional groups and generate a dictionary associating the molecule graph of that functional group with its name
for filename in glob.glob('*.xyz'):
   mol = Molecule.from_file(filename)
   name = filename.replace('.xyz', '')
   func_group_dict[name] = mol
print('Done!')
        
print('Loading mol_entries.pickle...')
with open('mol_entries.pickle', 'rb') as f: #loads pymatgen Molecule objects from pickle file
    mol_entries = pickle.load(f)
print('Done!')

print('Opening json...')  
third_name = 'reaction_tally'
third_entries = loadfn(third_name + ".json")
print('Done!')  

n = 0
test_dict = {}
entry_ids = set()
print('Naming Molecules...')
for reaction in third_entries["pathways"].keys():
    if third_entries["pathways"][reaction] > 500: #only add network products found >500 times
        for rxn in third_entries["reactions"].keys():
            if str(reaction) == rxn:
                molecule_dict = find_molecules_from_mpculeids(third_entries["reactions"][rxn], mol_entries) #want to associate an mpculeid with a name
                for l in molecule_dict.values():
                    for molecule in l:
                        if molecule.entry_id not in entry_ids: #i.e. we haven't named this molecule yet
                            name = name_molecule(molecule, func_group_dict) 
                            if stereoisomer_test(test_dict, molecule.entry_id, name): #returns empty dict if no stereoisomers already in dict
                                test_dict.update(stereoisomer_test(test_dict, molecule.entry_id, name))
                                if len(stereoisomer_test(test_dict, molecule.entry_id, name)) == 2: #If two entries are returned, this means we've named two stereoisomers #1 and #2, and need to remove the duplicate name from the dict
                                    test_dict.pop(name)
                            else:
                                test_dict[name] = molecule.entry_id
                            entry_ids.add(molecule.entry_id)
                                

print('Done naming molecules!')
end = time.time()
total = end - start
time_min = total/60
time_min = round(time_min, 2)
print('named ', len(test_dict), 'molecules and took', time_min, ' minutes')

dict_set = set(test_dict.values())
if entry_ids.difference(dict_set):
    print('Unnamed molecule hashes: ', entry_ids.difference(dict_set))
        
l_list = []
    
print('Testing for duplicates...')
for name, h in test_dict.items(): 
    for na, ha in test_dict.items():
        if na == name and ha != h:
            l_list.append((ha, h))

print('# of duplicate entries found: ', len(l_list))
if l_list:
    print(l_list)
print('Done!') 

dumpfn(test_dict, 'named_molecules.json')

# print('Done! ', len(mpcule_ids), ' reactions total')

# kinetiscope_reaction_list = []
