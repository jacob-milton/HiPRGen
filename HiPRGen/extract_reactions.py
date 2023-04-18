# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 13:13:41 2023

@author: JRMilton
"""
from monty.serialization import loadfn
import operator
import pickle

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
print('Done!')

def resonant_reaction(reaction_dict, added_hashes):
    #take in reaction_dict with {[reactant_hashes],[product_hashes]:reaction_charge}--a dictionary containing a tuple of reactants as the key and products as the values
    for reaction in added_hashes.keys():  #for each reaction currently in mpculids,
        if list(reaction_dict.keys()) == reaction: #compare the reactant hashes of the new reaction and the old one
            if added_hashes[reaction] == list(reaction_dict.hashes()): #if they're the same, compare the product hashes 
                return True

    return False

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
                    reactants = first_entries["reactions"][rxn]["reactants"]
                    products = first_entries["reactions"][rxn]["products"]
                    participants = [reactants, products]
                    reaction_charges = []
                    for side in participants:
                        for species in side: #takes the list of ids, finds their corresponding mol_entries, which have their charges and hashes
                            for mol_entry in mol_entries:
                                m_id = mol_entry.entry_id
                                if m_id == species:
                                    species_charge = mol_entry.charge
                                    print(species_charge)
                                    print(type(species_charge))
                                    reaction_charges.append(species_charge)
                                    species_hash = mol_entry.covalent_hash
                                    species_index = side.index(species)
                                    side[species_index] = species_hash
                    for side in participants:
                        side.sort()
                        side = tuple(side)
                    reaction_dict = {}
                    participants = tuple(participants)
                    reaction_dict[participants] = sum(reaction_charges)
                    if not resonant_reaction(reaction_dict, added_hashes):
                        added.append(reaction)
                        mpcule_ids.append(first_entries["reactions"][rxn])
                        added_hashes[participants] = sum(reaction_charges)
                        
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