# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 13:19:00 2024

@author: jacob
"""

from Rxn_classes import HiPRGen_Reaction

def find_charge(mpculeid):
    """
    This function simply returns the charge of a given species as an interger
    based on its mpculeid.

    Parameters
    ----------
    mpculeid : string
        a string of the form: "graph_hash-formula-charge-spin" 

    Returns
    -------
    charge : int
        charge of the species
    """
    
    charge_str = mpculeid.split("-")[2]
    if "m" in charge_str: #m stands for minus in the string
        charge = -int(charge_str.replace("m", ""))
        
    else:
        charge = int(charge_str)
        
    return charge

def create_reaction_dict(reaction_components, molecule_entry_dict):
    # Initialize total_charge and hash conversion
    total_charge = 0
    sorted_species_hashes = []

    for component_side in reaction_components:
        component_hashes = []
        for mpculeid in component_side:
            # molecule_entry = molecule_entry_dict.get(molecule_id)
            # if molecule_entry:
            total_charge += find_charge(mpculeid)
            graph_hash = mpculeid.split('-')[0]
            component_hashes.append(graph_hash)
        sorted_species_hashes.append(sorted(component_hashes))

    # Create the reaction dictionary with the sorted tuples as keys
    reaction_dict = {tuple(map(tuple, sorted_species_hashes)): total_charge}

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

reaction_list = [] #list of reactions we want to save of the form: [{'reactants': [mpculid, mpculeid], 'products': etc.}]
rxns_added = set() #keeps track of reactions we have already added to prevent redundancy
added_hashes = {}
first_name = 'reaction_tally_p1'
first_entries = loadfn(first_name + ".json") #loads json as a dictionary whose keys are mol ids and values are
                                                #dictionaries whose keys are labels of values

print('Loading mol_entries.pickle...')
with open('mol_entries.pickle', 'rb') as f: #loads mol_entries from pickle file
    mol_entries = pickle.load(f)
print('Done!')

# Create a mapping from molecule id to its properties for quick lookup
mol_entry_dict = {entry.entry_id: entry for entry in mol_entries}
reaction_indicies = first_entries["pathways"].keys()
reaction_index_dict = first_entries["reactions"]
    
print('Adding reactions from phase 1...')   
for reaction_index in reaction_indicies:
    reaction = reaction_index_dict.get(reaction_index, False)
    assert bool(reaction) == True
    is_ionization_reaction = False
    for species in reaction.values():
        if species == None:
            is_ionization_reaction = True
    if not is_ionization_reaction:
        new_rxn = HiPRGen_Reaction(reaction,phase=1)
        if new_rxn.reaction_hash not in rxns_added:
            reaction_dict = create_reaction_dict(reaction, mol_entry_dict)
            if not resonant_reaction(reaction_dict, added_hashes):
                rxns_added.add(new_rxn.reaction_hash)
                reaction_list.append(reaction)
                # added_hashes.update(reaction_dict.items())
            
    # for rxn in first_entries["reactions"].keys():
    #     if str(reaction) == rxn:
             #removes reactions occuring between a reactant/product and an "electron"
            # for entry in first_entries["reactions"][rxn].values():
            #     for species in entry:
            #         if species == None:
            #             electron_test = True
            # if not electron_test:
                # if reaction not in added: #don't add reactions twice
                #     reactants = first_entries["reactions"][rxn]["reactants"]
                #     products = first_entries["reactions"][rxn]["products"]
                #     participants = [reactants, products] 
                    
                    
                        
                        
print('Done! ', len(mpcule_ids), ' added')


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
                    reaction_dict = create_reaction_dict(participants, mol_entry_dict)
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
                reaction_dict = create_reaction_dict(participants, mol_entry_dict)
                if not resonant_reaction(reaction_dict, added_hashes):
                    if not reverse_reaction(reaction_dict, added_hashes):
                        if not charge_transfer_reaction(reaction_dict):
                            mpcule_ids.append(third_entries["reactions"][reaction])
                            added.append(reaction)
                            added_hashes.update(reaction_dict.items())

print('Done! ', len(mpcule_ids), ' reactions total')
dumpfn(mpcule_ids, 'euvl_TSreactions_041823.json')