# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 13:13:41 2023

@author: JRMilton
"""
from monty.serialization import loadfn
import operator
import pickle
from species_filter import really_covalent_isomorphic

# #Provide Evan with a list of the form:

# [{"reactants":[MPCuleID_A, MPCuleID_B], "products":[MPCuleID_C, MPCuleID_D]}, etc]


# 1. all non-electron-containing from phase 1 tally
# 2. for each network product, 10 most abundant one-step pathways
# 3. all non-charge-transfer from phase 2 tally with >= 500 occurences

# from all of those, ensure no duplicates + no duplicates modulo species charge

# (A + B- -> C- + D
# A- + B -> C- + D
# A- + B -> C + D-
# A + B- -> C + D- 

# Only one reaction!)


# total should be less than 7000
# JSONs contain reaction IDs, rn.sqlite contains mapping from reaction ID to species IDs, mol_entries.pickle contains mapping from species IDs to 
mpcule_ids = []
added =[]
added_hashes = {}
first_name = 'reaction_tally_p1'
first_entries = loadfn(first_name + ".json") #loads json as a dictionary whose keys are mol ids and values are
                                                #dictionaries whose keys are labels of values

print('Loading mol_entries.pickle...')
with open('mol_entries.pickle', 'rb') as f:
    mol_entries = pickle.load(f)
    print(mol_entries)
print('Done!')

def resonant_reaction(reaction_dict, added_hashes):
    #take in mpculeids for reaction we're testing--a dictionary containing a tuple of reactants as the key and products as the values
    hash_dict = {}
    reactants = []
    products = []

    for reactant in mpcule_dict.keys(): #find graph hashes corresponding to those mpculeids
        r_hash = mol_entries[reactant].covalent_hash
        reactants.append(r_hash)
    reactants.sort()
    for product in mpculte_dict.values:
        p_hash = mol_entries[product].covalent_hash
        reactants.append(r_hash)
    products.sort()

    reaction_dict[reactants] = products

    for reactant_pair in added_hashes.keys():  #for each reaction currently in mpculids,
        if reactant_pair == hash_dict.keys()[0]: #compare the reactant hashes of the new reaction and the old one
            product_pair = saved_hashes[reactant_pair] #if they're the same, compare the product hashes 
            if product_pair == hash_dict.values()[0]
                return {} #if those are the same, return true, otherwise return false

    return hash_dict

print('Adding reactions from phase 1...')   
for reaction in first_entries["pathways"].keys():
    for rxn in first_entries["reactions"].keys():
        if str(reaction) == rxn:
            electron_test = False
            for entry in first_entries["reactions"][rxn].values():
                for species in entry:
                    if species == None:
                        electron_test = True
            if not electron_test:
                if reaction not in added:
                    reaction_dict = {}
                    reactant_mpcule_list = list(first_entries["reactions"][rxn]["reactants"].values())
                    product_mpcule_list = list(first_entries["reactions"][rxn]["products"].values())
                    reaction_dict[reactants] = products
                    resonance = resonant_reaction(reaction_dict, added_hashes)
                    if resonance != {}:
                        added.append(reaction)
                        mpcule_ids.append(first_entries["reactions"][rxn])
                        added_hashes[list(resonance.keys())] = list(resonance.values())
                        

print('Done! ', len(mpcule_ids), ' added')


# print('Adding reactions for phase 2 network products...')                                             
# second_name = 'sink_report'
# second_entries = loadfn(second_name + ".json") #loads json as a dictionary whose keys are mol ids and values are
#                                                 #dictionaries whose keys are labels of values
# print(len(second_entries.keys()))
# for dictionary in second_entries.values():
#     species_index = dictionary["species_index"]
#     reaction_json = loadfn(str(species_index) + "_pathway.json") #dictionary containing two keys: pathways and reactions. 
#     pathways = reaction_json["pathways"]
#     pathways.sort(key = operator.itemgetter('frequency'), reverse = True)
#     top_pathways = []
#     n = 1
#     ten_saved = False
    
#     while not ten_saved:
#         for dictionary in pathways:
#             if int(dictionary['weight']) == n:
#                 top_pathways.append(dictionary['pathway'])
#                 if len(top_pathways) == 10:
#                     ten_saved = True
#                     break
#         n += 1
    
#     for rxn in top_pathways:
#         for reaction in rxn:
#             for num in reaction_json["reactions"].keys():
#                 if num == str(reaction):
#                     mpcule_ids.append(reaction_json["reactions"][num])
#     n = 1
                
# print('Done! ', len(mpcule_ids), ' reactions total')

# print('Adding reactions for phase 2 tally...')  
# third_name = 'reaction_tally'
# third_entries = loadfn(third_name + ".json")

# for reaction in third_entries["pathways"].keys():
#     if third_entries["pathways"][reaction] > 500:
#         for rxn in third_entries["reactions"].keys():
#             if str(reaction) == rxn:
#                 mpcule_ids.append(third_entries["reactions"][rxn])

# print('Done! ', len(mpcule_ids), ' reactions total')