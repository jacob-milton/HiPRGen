# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 13:13:41 2023

@author: JRMilton
"""
from monty.serialization import loadfn, dumpfn
import copy
import operator
import pickle
from openpyxl import workbook
from openpyxl import load_workbook

mpcule_ids = [] #list of reactions we want to save of the form: [{'reactants': [mpculid, mpculeid], 'products': etc.}]
added =[] #keeps track of reactions we have already added to prevent redundancy
added_hashes = {}
first_name = 'reaction_tally_p1'
first_entries = loadfn(first_name + ".json") #loads json as a dictionary whose keys are mol ids and values are
                                                #dictionaries whose keys are labels of values

print('Loading mol_entries.pickle...')
with open('mol_entries.pickle', 'rb') as f: #loads mol_entries from pickle file
    mol_entries = pickle.load(f)
print('Done!')

def create_reaction_dict(participants):
    """
    Takes a list of lists of reactant and product mpculeids and saves them to a dictionary
    whose key is a 2-tuple of the reactant and product graph hashes and whose value is the 
    sum of the charges of all species in the reaction. 

    participants: a list of the form: [[reactant1mpculid, reactant2mpculid][product1mpculid, product2mpculid]]

    Returns: reaction_dict
    """
    reaction_charges = []
    participants_copy = copy.deepcopy(participants)
    for side in participants_copy:
        for species in side: #takes the list of mpculeids, finds their corresponding mol_entries, which have their charges and hashes
            for mol_entry in mol_entries:
                m_id = mol_entry.entry_id
                if m_id == species:
                    species_charge = mol_entry.charge
                    reaction_charges.append(species_charge)
                    species_hash = mol_entry.covalent_hash
                    species_index = side.index(species)
                    side[species_index] = species_hash
    for side in participants_copy: #we want these sorted so we can directly compare reaction tuples later
        side.sort()
        side_tuple = tuple(side)
        side_index = participants_copy.index(side)
        participants_copy[side_index] = side_tuple
    reaction_dict = {}
    participants_copy = tuple(participants_copy) #we save these as tuples because dictionary keys must be hashable and lists are not
    reaction_dict[participants_copy] = sum(reaction_charges)
    return reaction_dict

def resonant_reaction(reaction_dict, added_hashes):
    """
    If a reaction is of the type:

    A + B -> C + D

    We assume that the transition state for this reaction will be the same regardless
    of the charges of the reactants or products, and thus do not add a reaction
    if a 'resonant reaction' is already in our list.

    reaction_dict: dictionary whose key is a 2-tuple of sorted graph hashes for a 
    reaction and whose value is the sum of all of the charges of the reaction.

    added_hashes: a dictionary containing reaction_dicts that have already been
    addded to our list to save

    Returns: True if the reaction we're testing resonantes with one already added
    to our list, and False otherwise.
    """
    for reaction in added_hashes.keys(): #compare the graph hashes of the reaction we're testing to all those already added
        for rxn in reaction_dict.keys():
            if rxn == reaction: #if the graphs of both products and reactants are the same,
                if added_hashes[reaction] == reaction_dict[rxn]: #and the sum of the charges is the same, return True.
                    return True

    return False

def reverse_reaction(reaction_dict, added_hashes):
    """
    If a reaction is of the form A+B -> C+D and we have already saved the reverse
    reaction, the reactions will have the same transition state and therefore
    only one needs to be saved.

    reaction_dict: dictionary of the form: {((reactant_hashes), (product_hashes)): sum of charges of all species in the reaction}

    added_hashes: a dictionary containing all reaction_dicts that have already been added to our reaction list

    Returns: True if the reverse reaction is already present in our list, False otherwise.
    """
    for reaction in reaction_dict.keys():
        reverse = (reaction[1], reaction[0])
        for rxn in added_hashes:
            if reverse == rxn and reaction_dict[reaction] == added_hashes[rxn]:
                return True

    return False

def charge_transfer_reaction(reaction_dict):
    """
    Removes reactions where no bonds are broken or formed, only electrons are transferred
    between reactants. 

    reaction_dict: dictionary of the form: {((reactant_hashes), (product_hashes)): sum of charges of all species in the reaction}

    Returns True if a reaction is electron transfer, False otherwise. 
    """
    for reaction in reaction_dict.keys():
        reactant_total_hashes = set(reaction[0])
        product_total_hashes = set(reaction[1])
        reactant_hash_list = list(reaction[0])
        product_hash_list = list(reaction[1])

    if len(reactant_total_hashes.intersection(product_total_hashes)) == 2: #if the products and reactants have the same graph, the reaction is electron transfer
        return True
    elif len(reactant_total_hashes.intersection(product_total_hashes)) == 1 and reactant_hash_list[0] == reactant_hash_list[1] and product_hash_list[0] == product_hash_list[1]:
        return True #finds special cases for electron transfer between two of the same reactant
    else:
        return False
        
    return False

# print('Adding reactions from phase 1...')   
# for reaction in first_entries["pathways"].keys():
#     for rxn in first_entries["reactions"].keys():
#         if str(reaction) == rxn:
#             electron_test = False #removes reactions occuring between a reactant/product and an "electron"
#             for entry in first_entries["reactions"][rxn].values():
#                 for species in entry:
#                     if species == None:
#                         electron_test = True
#             if not electron_test:
#                 if reaction not in added: #don't add reactions twice
#                     reactants = first_entries["reactions"][rxn]["reactants"]
#                     products = first_entries["reactions"][rxn]["products"]
#                     participants = [reactants, products] 
#                     reaction_dict = create_reaction_dict(participants)
#                     if not resonant_reaction(reaction_dict, added_hashes):
#                         added.append(reaction)
#                         mpcule_ids.append(first_entries["reactions"][rxn])
#                         added_hashes.update(reaction_dict.items())
                        
# print('Done! ', len(mpcule_ids), ' added')

print('Adding reactions for phase 2 network products...')                                             
second_name = 'sink_report'
second_entries = loadfn(second_name + ".json") 

for dictionary in second_entries.values(): #iterates through each network product
    species_index = dictionary["species_index"]
    reaction_json = loadfn(str(species_index) + "_pathway.json") #loads each network_product.json file as dictionary containing  
    pathways = reaction_json["pathways"]                         #two keys: pathways and reactions.
    pathways.sort(key = operator.itemgetter('frequency'), reverse = True) #sorts by frequency
    top_pathways = []
    n = 1
    ten_saved = False
    
    while not ten_saved: #makes sure we save top ten reactions even for network products that don't have 10 1-step reactions forming them
        for dictionary in pathways:
            if int(dictionary['weight']) == n:
                top_pathways.append(dictionary['pathway'])
                if len(top_pathways) == 10:
                    ten_saved = True
                    break
        n += 1
    
    for rxn in top_pathways:
        for reaction in rxn:
            for num in reaction_json["reactions"].keys():
                if num == str(reaction):
                    reactants = reaction_json["reactions"][num]["reactants"]
                    products = reaction_json["reactions"][num]["products"]
                    participants = [reactants, products] 
                    reaction_dict = create_reaction_dict(participants)
                    if not resonant_reaction(reaction_dict, added_hashes):
                        if not reverse_reaction(reaction_dict, added_hashes):
                            if not charge_transfer_reaction(reaction_dict):
                                mpcule_ids.append(reaction_json["reactions"][num])
                                added.append(num)
                                added_hashes.update(reaction_dict.items())
    n = 1
                
print('Done! ', len(mpcule_ids), ' reactions total')

print('Adding reactions for phase 2 tally...')  
third_name = 'reaction_tally'
third_entries = loadfn(third_name + ".json")

for reaction in third_entries["pathways"].keys():
    if third_entries["pathways"][reaction] > 500: #only add network products found >500 times
        for rxn in third_entries["reactions"].keys():
            if str(reaction) == rxn:
                reactants = third_entries["reactions"][reaction]["reactants"]
                products = third_entries["reactions"][reaction]["products"]
                participants = [reactants, products] 
                reaction_dict = create_reaction_dict(participants)
                if not resonant_reaction(reaction_dict, added_hashes):
                    if not reverse_reaction(reaction_dict, added_hashes):
                        if not charge_transfer_reaction(reaction_dict):
                            mpcule_ids.append(third_entries["reactions"][reaction])
                            added.append(reaction)
                            added_hashes.update(reaction_dict.items())

print('Done! ', len(mpcule_ids), ' reactions total')

kinetiscope_reaction_list = []

print('Converting mpculeids to formulas and entry_indcies...')
for reaction in mpcule_ids: #mpculeids is a list of dictionarys containing two keys: reactants and products with mpcule ids as the values
    new_reaction = {}
    reactant_list = []
    reactants = list(reaction['reactants'])
    for reactant in reactants:
        for mol_entry in mol_entries:
            m_id = mol_entry.entry_id
            if m_id == reactant:
                m_index = str(mol_entry.ind)
                hyphen_index = reactant.find('-')
                formula = reactant[hyphen_index+1:len(reactant)]
                reactant_list.append(formula + '(' + m_index + ')')
    new_reaction['reactants'] = tuple(reactant_list)
    product_list = []
    products = list(reaction['products'])
    for product in products:
        for mol_entry in mol_entries:
            m_id = mol_entry.entry_id
            if m_id == product:
                m_index = str(mol_entry.ind)
                formula = product[hyphen_index+1:len(product)]
                product_list.append(formula + '(' + m_index + ')')          
    new_reaction['products'] = tuple(product_list)
    kinetiscope_reaction_list.append(new_reaction)

template = load_workbook(filename = 'Kinetiscope_rxn_template.xlsx')
ts = template.active
cell_range = ts['A1':'R1']
new_file = openpyxl.workbook()

ws = new_sheet.active
                
                
# dumpfn(mpcule_ids, 'euvl_TSreactions_041823.json')