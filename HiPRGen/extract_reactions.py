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
# mpcule_ids = []
# first_name = input("Name of input json file: ")
# first_entries = loadfn(first_name + ".json") #loads json as a dictionary whose keys are mol ids and values are
#                                                 #dictionaries whose keys are labels of values

# for reaction in first_entries["pathways"].keys():
#     for rxn in first_entries["reactions"].keys():
#         if str(reaction) == rxn:
#             electron_test = False
#             for entry in first_entries["reactions"][rxn].values():
#                 for species in entry:
#                     if species == None:
#                         electron_test = True
#             if not electron_test:
#                 mpcule_ids.append(first_entries["reactions"][rxn])
# print(len(mpcule_ids))
                                               
# second_name = input("Name of input json file: ")
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
                
# print(len(mpcule_ids))

# third_name = input("Name of input json file: ")
# third_entries = loadfn(third_name + ".json")

# for reaction in third_entries["pathways"].keys():
#     if third_entries["pathways"][reaction] > 500:
#         for rxn in third_entries["reactions"].keys():
#             if str(reaction) == rxn:
#                 mpcule_ids.append(third_entries["reactions"][rxn])

# print(len(mpcule_ids))

with open('mol_entries.pickle', 'rb') as f:
    mol_entries = pickle.load(f)
    print(mol_entries)
