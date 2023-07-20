# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import os
import copy
from itertools import combinations_with_replacement
import time

from monty.serialization import dumpfn, loadfn
import numpy as np

from openbabel import openbabel as ob
from openbabel import pybel as pb

from pymatgen.analysis.graphs import MoleculeGraph, MolGraphSplitError
from pymatgen.core.structure import Molecule, IMolecule
from pymatgen.core.sites import Site
from pymatgen.analysis.local_env import OpenBabelNN, metal_edge_extender

from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
import networkx as nx

from atomate.qchem.database import QChemCalcDb
import scine_molassembler #Molassembler is a C++ library that aims to facilitate crossings between Cartesian and graph representations of molecules. 

def add_hot_atom_tags(mol_graph: MoleculeGraph, nbo_spins): #identifies radical rich sites
    """
    Adds metadata to a given MoleculeGraph object describing the location of radical
    sites present within a given fragment if they are present

    Parameters
    ----------
    mol_graph : MoleculeGraph
        The fragment to be labelled
    nbo_spins : List
        A list describing the natural bond orbital spins of different sites
        in a molecule

    Returns
    -------
    mol_graph : MoleculeGraph
        The molecule graph with the added metadata denoting what sites in this fragment
        are "radical rich"
        
    """
    hot_atoms = []
    for ii in range(len(mol_graph.molecule)): #adds len(mol_graph.molecule) 0's to the hot_atoms list, where len(mol_graph.molecule) is the number of atoms in
        hot_atoms.append(0)                   #the molecule
    for atom_ind in nbo_spins:
        if nbo_spins[atom_ind] > 0.09 and str(mol_graph.molecule[int(atom_ind)].specie) != "H": #how did we choose 0.09 as the cutoff?
            if str(mol_graph.molecule[int(atom_ind)].specie) == "C": #Prevents tetrahedral carbons from reacting
                connected_sites = mol_graph.get_connected_sites(atom_ind)
                if len(connected_sites) > 3:
                    print("previously hot 4-neighbor carbon caught!")
                    continue
            hot_atoms[int(atom_ind)] = 1 #if there is a radical and it is not tertiary, this labels it as a "hot atom"
    mol_graph.molecule.add_site_property("hot", hot_atoms)
    return mol_graph


def generate_obmol(mol_graph: MoleculeGraph): #Creates an openbabel OBMol object from out given pymatgen Molecule object
    """
    Creates an openbabel OBMol object from a given pymatgen MoleculeGraph

    Parameters
    ----------
    mol_graph : MoleculeGraph  

    Returns
    -------
    obmol: OBMol
    
    """
    obmol = ob.OBMol() 
    obmol.BeginModify() 
    for site in mol_graph.molecule:
        coords = [c for c in site.coords]
        atomno = site.specie.Z
        obatom = ob.OBAtom()
        obatom.thisown = 0
        obatom.SetAtomicNum(atomno)
        obatom.SetVector(*coords)
        obmol.AddAtom(obatom)
        del obatom

    obmol.SetTotalSpinMultiplicity(mol_graph.molecule.spin_multiplicity)
    obmol.SetTotalCharge(int(mol_graph.molecule.charge))
    obmol.Center()

    for edge in list(mol_graph.graph.edges()):
        # Always assume single bond
        obmol.AddBond(edge[0] + 1, edge[1] + 1, 1)

    obmol.EndModify()
    return obmol

def generate_scinemol(mol_graph: MoleculeGraph, directory):
    """
    Creates a scine_molassembler Molecule instance from a given pymatgen MoleculeGraph
    Parameters
    ----------
    mol_graph : MoleculeGraph
    directory : String

    Returns
    -------
    scine_mol : scine_molassembler Molecule instance

    """
    mol_graph.molecule.to("xyz", os.path.join(directory,"tmp.xyz")) #in pymatgen, the .to function for Istructure objects (the parent class of molecules) outputs the molecule to a file or string. In this case the name of the file will be "xyz" and it will be written as an XYZ file
    scine_mol = scine_molassembler.io.read(os.path.join(directory,"tmp.xyz")) #reads our .xyz file and generates a scine_molassmbler Molecule instance from it 
    os.remove(os.path.join(directory,"tmp.xyz"))
    return scine_mol

def combine_molecules_molassembler(molgraph_1: MoleculeGraph, molgraph_2: MoleculeGraph, #this is recombination for the scine_molassember Molecule instances
                                   index_1: int, index_2: int):
    """
    Takes the MoleculeGraph objects of two fragments and combines them into a single new 
    MoleculeGraph, ensuring that the fragments do not overlap with each other via a conformational
    search
    
    Parameters
    ----------
    molgraph_1 : MoleculeGraph
    molgraph_2 : MoleculeGraph
    index_1 : int
    index_2 : int

    Returns
    -------
    new_mg : MoleculeGraph
        A MoleculeGraph resulting from the combination of the two input fragments

    """
    scinemol1 = generate_scinemol(molgraph_1, "/global/home/groups/lr_mp/smblau/EUVL/dec_recomb/")
    scinemol2 = generate_scinemol(molgraph_2, "/global/home/groups/lr_mp/smblau/EUVL/dec_recomb/")
    combined_scinemol = scine_molassembler.editing.connect(scinemol1, scinemol2, index_1, index_2, scine_molassembler.BondType(1)) #Pretty sure BondType 1 in scine_molassembler is a double bond and 0 would be single...
    conformations = scine_molassembler.dg.generate_ensemble(combined_scinemol, 50, 1010) #Generate a set of 50 3D positions for a molecule...somehow           
    conformation = None
    for confo in conformations:
        if not isinstance(confo, scine_molassembler.dg.Error): #discards conformations that produce an error
            conformation = confo
            break
    if conformation is not None: #If we find a conformation, generate a new mol_graph for it with a bond between the two "reacting" atoms
        new_mg = combine_mol_graphs(molgraph_1, molgraph_2)
        new_mg.add_edge(index_1, index_2 + len(molgraph_1.molecule))
        
        for ii, pos in enumerate(conformation):
            new_mg.molecule[ii].coords = pos*0.529177 #0.529177 is just used to convert units

        return new_mg
    else:
        print("FAILED")
        return None
    
def combine_mol_graphs(molgraph_1: MoleculeGraph, molgraph_2: MoleculeGraph):
    """
    Combines two mol_graphs, ensuring that they don't overlap in the new graph
    
    Parameters
    ----------
    molgraph_1 : MoleculeGraph
    molgraph_2 : MoleculeGraph

    Returns
    -------
    copy_1 : MoleculeGraph
        A new molecule graph resulting from the recombination of the two fragments

    """
    radius_1 = np.amax(molgraph_1.molecule.distance_matrix) #finds the maximum distance between two atoms in each mol_graph and defines that as the radius
    radius_2 = np.amax(molgraph_2.molecule.distance_matrix) #of a sphere containing the molecule/ion

    copy_1 = copy.deepcopy(molgraph_1)
    if "hot" in copy_1.molecule.site_properties: #removes the hot atom properties we assigned earlier because they aren't true for this new molecule
        copy_1.molecule.remove_site_property("hot")
    copy_1.molecule.translate_sites(list(range(len(molgraph_1.molecule))), #Rotate all sites about the center of mass of the molecule
                                    -1 * molgraph_1.molecule.center_of_mass)

    copy_2 = copy.deepcopy(molgraph_2)
    if "hot" in copy_2.molecule.site_properties:
        copy_2.molecule.remove_site_property("hot")
    copy_2.molecule.translate_sites(list(range(len(copy_2.molecule))), #Translate all sites by the inverse of the center of mass of the molecule added to the sum of both radii and 1.0?...
                                    -1 * copy_2.molecule.center_of_mass + np.array([radius_1 + radius_2 + 1.0, 0.0, 0.0]))

    for site in copy_2.molecule: #inserts all sites in copy_2.molecule into copy_1
        copy_1.insert_node(len(copy_1.molecule), site.specie, site.coords)

    for edge in copy_2.graph.edges(): #i.e. iterating over bonds in the graph
        side_1 = edge[0] + len(molgraph_1.molecule) #edge is a 2-tuple containing two indicies that in this case represent atoms
        side_2 = edge[1] + len(molgraph_1.molecule) #so here we're increasing the indicies so that they don't overlap with molgraph_1
        copy_1.add_edge(side_1, side_2)

    copy_1.molecule.set_charge_and_spin(molgraph_1.molecule.charge + molgraph_2.molecule.charge)

    return copy_1

# def combine_molecules(molgraph_1: MoleculeGraph, molgraph_2: MoleculeGraph,
#                       index_1: int, index_2: int):
#     """

#     Parameters
#     ----------
#     molgraph_1 : MoleculeGraph
#     molgraph_2 : MoleculeGraph
#     index_1 : int
#         Index of one reactive atom to be bonded to another in the recombined
#         MoleculeGraph
#     index_2 : int

#     Returns
#     -------
#     new_mg : MoleculeGraph
#         DESCRIPTION.

#     """
#     new_mg = combine_mol_graphs(molgraph_1, molgraph_2)

#     obmol = generate_obmol(new_mg)
#     builder = ob.OBBuilder()
#     builder.Connect(obmol, index_1 + 1, index_2 + len(molgraph_1.molecule) + 1) #Atoms index_1+1 and index_2+ 1 are part of two fragments that are not connected in mol. 
#                                                                                 #Connect will translate and rotate the fragment that contains b so that a and b are seperated by a bond. This bond is also added.
#     new_mg.add_edge(index_1, index_2 + len(molgraph_1.molecule)) #somehow adding new bonds to the new mol_graph?

#     for ii, atom in enumerate(ob.OBMolAtomIter(obmol)): #iterates over our new mol_graph and adds all of its atoms' coordinates
#         new_mg.molecule[ii].coords = [atom.GetX(), atom.GetY(), atom.GetZ()]

#     return new_mg

def identify_connectable_atoms(mol_graphs): 
    """
    Generates two lists associated with each MoleculeGraph objects, describing
    the indicies of "hot atoms" and "underbonded atoms" within that MoleculeGraph.
    "Hot atoms" are radical-rich sites, found as described in add_hot_atom_tags, while "
    underbonded" sites are those that have less than the ideal number of bonds for a given atom.
    Note that sulfur will always be classified as underbonded under this paradigm. 
    
    :param mol_graphs: list of MoleculeGraph objects
    :return: a list of lists of indicies of underbonded and hot atoms for each mol   
    """

    bond_max = {"C": 4,
                "P": 5,
                "S": 6,
                "O": 2,
                "N": 3,
                "Cl": 1,
                "F": 1
    }

    underbonded_atoms_index_list = list() 
    hot_atoms_index_list = list()
    for i, mol_graph in enumerate(mol_graphs): 
        if not nx.is_connected(mol_graph.graph.to_undirected()): #from wolfram: in a connected graph, there is a path from any point to any other point in the graph--i.e. in chemical terms all atoms are bonded in a connected graph. If that is not the case we skip it.
            continue
        underbonded_atoms_in_mol = list() #this will collect a list of indicies of underbonded atoms in a fragment
        hot_atoms_in_mol = list()
        num_atoms = len(mol_graph.molecule)
        hot_atoms = mol_graph.molecule.site_properties["hot"]

        two_neighbor_carbon_found = False #searches for carbons with only two neighborns 
        for j in range(num_atoms):
            if str(mol_graph.molecule[int(j)].specie) == "C":
                connected_sites = mol_graph.get_connected_sites(j)
                if len(connected_sites) == 2:
                    nitrogen_found = False #that are not in cyano groups
                    for site in connected_sites:
                        if str(site.site.specie) == "N":
                            nitrogen_found = True
                    if not nitrogen_found:
                        two_neighbor_carbon_found = True
                        underbonded_atoms_in_mol.append(j)
                        break #if we have a carbon with only two neighbors, immediately return that as the only hot atom
                elif len(connected_sites) == 1: 
                    two_neighbor_carbon_found = True
                    underbonded_atoms_in_mol.append(j)
                    break

        if not two_neighbor_carbon_found:
            for j in range(num_atoms):
                connected_sites = mol_graph.get_connected_sites(j)
                num_connected_sites = len(connected_sites)
                element = str(mol_graph.molecule[j].specie)
                if element == "C": #looks for rings
                    total_weight = 0
                    for site in connected_sites: 
                        total_weight += site.weight #weight is always one
                    if total_weight == 3 and num_connected_sites == 3:
                        rings = mol_graph.find_rings(including=[j])
                        if len(rings) > 0:
                            total_weight += 1
                    num_connected_sites = total_weight

            
                if element in ["Li", "Mg", "H"]: #these species are underbonded if they aren't bonded to anything
                    if num_connected_sites == 0 and num_atoms == 1:
                        underbonded_atoms_in_mol.append(j)

                else:
                    metal_count = 0

                    for k, site in enumerate(connected_sites):
                        if str(site.site.specie) in ["Li", "Mg"]:
                            metal_count += 1

                    if num_connected_sites - metal_count < bond_max[element]:
                        underbonded_atoms_in_mol.append(j)
                    elif hot_atoms[j] == 1:
                        hot_atoms_in_mol.append(j)
                    
        underbonded_atoms_index_list.append(underbonded_atoms_in_mol)
        hot_atoms_index_list.append(hot_atoms_in_mol)
        
    return underbonded_atoms_index_list, hot_atoms_index_list


def generate_combinations(db_entries, mol_graphs, directory, age_list=None):
    """
    Generate all combination of molecule/atom indices that can participate in recombination,
    by looping through all molecule pairs(including a mol and itself) and all connectable heavy atoms in each molecule.
    :param db_entries: [MoleculeGraph] already calculated, checked against to avoid duplicates. Typically set equivalent to mol_graphs.
    :param mol_graphs: [MoleculeGraph] to participate in recombination
    :param directory: String location to print combinations.txt and mol_graphs_recombination.json
    :param age_list: ["old" OR "new"] where old species cannot recombine with old species, to avoid duplicates from a previous recombination campaign
    :return: list of string [ 'mol1_index'+'_'+'mol2_index'+'_'+'atom1_index'+'_'+'atom2_index'].
    """

    combinations_file = os.path.join(directory, "dec_combinations.txt")
    mol_graphs_file = os.path.join(directory, "dec_mol_graphs_recombination.json")

    db_hashes = [] #add hashes for our new MoleculeGraphs
    for mg in db_entries:
        if nx.is_connected(mg.graph.to_undirected()):
            mg_hash = weisfeiler_lehman_graph_hash(mg.graph.to_undirected(), node_attr="specie")
            db_hashes.append({"hash":mg_hash, "charge":mg.molecule.charge})
    
    with open(combinations_file, "w") as combos:
        combos.write("mol_1\tatom_1\tmol_2\tatom_2\n") #This atom in this molecule can recombine with that atom in that molecule
        final_list = []
        final_hashes_list = []
        underbonded_atoms_index_list, hot_atoms_index_list = identify_connectable_atoms(mol_graphs) #returns two lists
        num_mols = len(mol_graphs)
        all_mol_pair_index = list(combinations_with_replacement(range(num_mols), 2)) #this is a list of tuples containing all possible combinations of indicies meaning any two molecules can recombine
        for pair_index in all_mol_pair_index:
            mol_graph1 = mol_graphs[pair_index[0]] #this pair index is the molecule index
            mol_graph2 = mol_graphs[pair_index[1]]

            total_charge = mol_graph1.molecule.charge + mol_graph2.molecule.charge
            # total_electrons = mol_graph1.molecule._nelectrons + mol_graph2.molecule._nelectrons

            if int(total_charge) not in {-1, 0, 1}:
                continue
            #elif total_electrons % 2 != 0:
            #    continue
            elif age_list:
                if age_list[pair_index[0]] == "old" and age_list[pair_index[1]] == "old":#we've already combined these fragments so we ignore them
                    continue
            
            underbonded_atoms_1 = underbonded_atoms_index_list[pair_index[0]] 
            underbonded_atoms_2 = underbonded_atoms_index_list[pair_index[1]]
            hot_atoms_1 = hot_atoms_index_list[pair_index[0]]
            hot_atoms_2 = hot_atoms_index_list[pair_index[1]]
            connectable_atoms_1 = underbonded_atoms_1 + hot_atoms_1
            connectable_atoms_2 = underbonded_atoms_2 + hot_atoms_2
            if len(connectable_atoms_1) == 0 or len(connectable_atoms_2) == 0:
                continue
            else:
                for i, atom1 in enumerate(connectable_atoms_1): #recombines radical rich sites with underbonded sites and enumerates over all possiblities
                    for j, atom2 in enumerate(connectable_atoms_2):
                        if atom1 in hot_atoms_1 and atom2 in hot_atoms_2: #forbid radical rich recombine with radical rich
                            continue
                        specie1 = str(mol_graph1.molecule[atom1].specie)
                        specie2 = str(mol_graph2.molecule[atom2].specie)

                        if specie1 in ["Li", "Mg"] and specie2 in ["Li", "Mg"]:
                            continue

                        combined_mol_graph = combine_mol_graphs(mol_graph1, mol_graph2) #generates a new MoleculeGraph object for the recombinant
                        combined_mol_graph.add_edge(atom1, atom2 + len(mol_graph1.molecule)) #adds the new bond in that recombinant

                        # Remove species with a bare aromatic ring carbon
                        #bare_aromatic_carbon_found = False
                        #rings = combined_mol_graph.find_rings()
                        #if len(rings) > 0:
                        #    for ring in rings:
                        #        for edge in ring:
                        #            k = edge[0]
                        #            if str(combined_mol_graph.molecule[k].specie) == "C":
                        #                connected_sites = combined_mol_graph.get_connected_sites(k)
                        #                if len(connected_sites) == 2:
                        #                    bare_aromatic_carbon_found = True
                        #                    break
                        #if bare_aromatic_carbon_found:
                        #    continue
                        
                        if not nx.is_connected(combined_mol_graph.graph.to_undirected()):
                            raise RuntimeError("Disconnected combined mol graph found! Exiting...")
                        
                        cmg_hash = weisfeiler_lehman_graph_hash(combined_mol_graph.graph.to_undirected(), node_attr="specie") #makes a new graph hash for our recombinant graph
                        cmg_charge = combined_mol_graph.molecule.charge
                        
                        match = False #this test makes sure we're not adding any redundant recombinants
                        for entry in db_hashes: 
                            if cmg_hash == entry["hash"] and cmg_charge == entry["charge"]:
                                match = True
                                break
                        if match:
                            continue
                        else:
                            index = None
                            for ii, entry in enumerate(final_hashes_list):
                                if cmg_hash == entry["hash"] and cmg_charge == entry["charge"]:
                                    index = ii
                                    break
                            if index is None: #generates a new index for our new recombinant
                                index = len(final_list)
                                final_list.append(combined_mol_graph)
                                final_hashes_list.append({"hash":cmg_hash, "charge":cmg_charge})
                                
                                combos.write("{}\t{}\t{}\t{}\t{}\n".format(pair_index[0],
                                                                           atom1,
                                                                           pair_index[1],
                                                                           atom2,
                                                                           index))

    dumpfn(final_list, mol_graphs_file)

    return final_list

def parse_combinations_file(filename):
    """
    Reads through a text file (typically) generated from the generate_cominbations function
    and returns reactions as a list 

    Parameters
    ----------
    filename : String

    Returns
    -------
    reactions : List
        A list of reactions of the form [mol1_index, mol2_index, atom1_index, atom2_index]

    """
    with open(filename) as combo_file:
        lines = combo_file.readlines()

        reactions = list()
        # Skip first line - header
        for line in lines[1:]:
            line_parsed = [int(x) for x in line.strip().split("\t")] #The .strip() method for a string removes all spaces from the beginning and end, while split returns a list
            reactions.append(tuple(line_parsed))                     #Split according to the seperator ("\t")

        return reactions


# def find_underbonded_atoms(mol_graphs):
#     bond_max = {"C": 4,
#                 "P": 5,
#                 "S": 6,
#                 "O": 2,
#                 "N": 3,
#                 "Cl": 1,
#                 "F": 1
#     }
#     underbonded_atoms_index_list = []
#     for i, mol_graph in enumerate(mol_graphs): #what is the point of this index tracking if we don't use it in the for loop?
#         if not nx.is_connected(mol_graph.graph.to_undirected()): #from wolfram: in a connected graph, there is a path from any point to any other point in the graph--i.e. in chemical terms all atoms are bonded in a connected graph. If that is not the case we skip it.
#             continue
#         underbonded_atoms_in_mol = [] #this will collect a list of indicies of underbonded atoms in a fragment
#         num_atoms = len(mol_graph.molecule)

#         # two_neighbor_carbon_found = False #searches for carbons with only two neighborns 
#         # for j in range(num_atoms):
#         #     if str(mol_graph.molecule[int(j)].specie) == "C":
#         #         connected_sites = mol_graph.get_connected_sites(j)
#         #         if len(connected_sites) == 2:
#         #             nitrogen_found = False #that are not in cyano groups
#         #             for site in connected_sites:
#         #                 if str(site.site.specie) == "N":
#         #                     nitrogen_found = True
#         #             if not nitrogen_found:
#         #                 two_neighbor_carbon_found = True
#         #                 underbonded_atoms_in_mol.append(j)
#         #                 break #if we have a carbon with only two neighbors, immediately return that as the only hot atom
#         #         elif len(connected_sites) == 1: #when would a carbon only have one neighbor?...
#         #             two_neighbor_carbon_found = True
#         #             underbonded_atoms_in_mol.append(j)
#         #             break

#         # if not two_neighbor_carbon_found:
#             for j in range(num_atoms):
#                 connected_sites = mol_graph.get_connected_sites(j)
#                 num_connected_sites = len(connected_sites)
#                 element = str(mol_graph.molecule[j].specie)
#                 if element == "C": #looks for rings
#                     total_weight = 0
#                     for site in connected_sites: 
#                         total_weight += site.weight #weight is always one
#                     if total_weight == 3 and num_connected_sites == 3:
#                         rings = mol_graph.find_rings(including=[j])
#                         if len(rings) > 0:
#                             total_weight += 1
#                     num_connected_sites = total_weight

            
#                 if element in ["Li", "Mg", "H"]: #these species are underbonded if they aren't bonded to anything
#                     if num_connected_sites == 0 and num_atoms == 1:
#                         underbonded_atoms_in_mol.append(j)

#                 else:
#                     metal_count = 0

#                     for k, site in enumerate(connected_sites):
#                         if str(site.site.specie) in ["Li", "Mg"]:
#                             metal_count += 1

#                     if num_connected_sites - metal_count < bond_max[element]:
#                         underbonded_atoms_in_mol.append(j)
#                     # elif hot_atoms[j] == 1:
#                     #     hot_atoms_in_mol.append(j)
                    
#         underbonded_atoms_index_list.append(underbonded_atoms_in_mol)
#         # hot_atoms_index_list.append(hot_atoms_in_mol)

def gen_double_bonded_products(mol_graphs):
    '''
    Takes fragments that have already been generated and generates double-bonded
    species from them. If, for example, a fragment is a radical, breaking any bond
    adjacent to the radical site will generate a species containing a double bond, 
    and this function will generate those "recombinants" from fragment MoleculeGraph
    objects.
    
    Parameters
    ----------
    mol_graphs : list of MoleculeGraph objects

    Returns
    -------
    fragments : dictionary
        A dictionary mapping the IMolecule object of a fragment to its possible
        "recombinants" containing double bonds
    '''
    fragments = {}
    underbonded_atoms_index_list, hot_atoms_index_list = identify_connectable_atoms(mol_graphs) #find all fragments that contain underbonded atoms
    for i, atom_index_list in enumerate(underbonded_atoms_index_list):
        if atom_index_list: #fires only if a fragment contains "underbonded" atoms
            mol_graph = mol_graphs[i]
            for index in atom_index_list: #generate recombinants by breaking only bonds adjacent to the underbonded atoms
                pi_bond_found = False
                if mol_graph.find_rings(): #making "recombinants" in aromatic rings would generate arynes which we do not want
                    if str(mol_graph.molecule[int(index)].specie) == "C" and len(mol_graph.get_connected_sites(index)) == 2:
                        pi_bond_found = True
                elif str(mol_graph.molecule[int(index)].specie) == "C" and len(mol_graph.get_connected_sites(index)) == 3:
                    pi_bond_found = False
                    neighbors = mol_graph.get_connected_sites(index)
                    for neighbor in neighbors:
                        if str(mol_graph.molecule[int(neighbor[2])].specie) == "C" and len(mol_graph.get_connected_sites(neighbor[2])) == 3:
                            pi_bond_found = True
                        elif str(mol_graph.molecule[int(neighbor[2])].specie) == "N" and len(mol_graph.get_connected_sites(neighbor[2])) != 3:
                           pi_bond_found = True
                        elif str(mol_graph.molecule[int(neighbor[2])].specie) == "O" and len(mol_graph.get_connected_sites(neighbor[2])) == 1:
                           pi_bond_found = True
                if not pi_bond_found:
                    underbonded_atom_neighbors = mol_graph.get_connected_sites(index) #connectedsites returns a named tuple of neighbors of site index: periodic_site, jimage, index, weight
                    for neighbor in underbonded_atom_neighbors:
                        to_disconnect = mol_graph.get_connected_sites(neighbor[2]) #finds neighbors of atom bonded to the underbonded site
                        mol_fragments = []
                        for atom in to_disconnect:
                            bond = [(neighbor[2], atom[2])]  
                            try:
                                double_bonded_fragment = mol_graph.split_molecule_subgraphs(bond, allow_reverse=True) #returns a list of MoleculeGraphobjects
                                mol_fragments.append(double_bonded_fragment)
                            except:
                                continue
                        if mol_fragments:
                            mol_graph_dict = mol_graph.molecule.as_dict() #convert our Molecule object to an IMolecule object, as Imolecule objects are hashable
                            Imol = IMolecule.from_dict(mol_graph_dict)
                            fragments[Imol] = mol_fragments
    return fragments
                      
print("Loading NBO tasks from DB") #loading NBO info from fragment DFT outputs

db_file = os.path.join(os.environ["HOME"],"/global/home/users/jrmilton/db.json") #points to the directory of the NBO json
mmdb = QChemCalcDb.from_db_file(db_file, admin=True) #converts that db_file to a QChemCalcDb
old_target_entries = list(mmdb.collection.find({"tags.class":"EUVL_prod"})) #all of these create lists of fragments to iterate through
new_target_entries = list(mmdb.collection.find({"tags.class":"EUVL_Nov_SP"}))
pfbs_entries = list(mmdb.collection.find({"tags.class":"EUVL_pfbs_SP"}))
cnbz_entries = list(mmdb.collection.find({"tags.class":"EUVL_4cnbz_SP"}))
NBO_entries = []
age_list = []
for entry in old_target_entries:
    if "NBO" in entry["task_label"]:
        mg = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(entry["output"]["initial_molecule"]), OpenBabelNN())
        if nx.is_connected(mg.graph.to_undirected()):
            NBO_entries.append(entry)
            age_list.append("old")
print(len(NBO_entries), "old entries")
nov_sp = 0
for entry in new_target_entries:
    mg = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(entry["output"]["initial_molecule"]), OpenBabelNN()) #creates a MoleculeGraph object 
    if nx.is_connected(mg.graph.to_undirected()):
        if "F" in entry["formula_alphabetical"]:
            if entry["formula_alphabetical"] == "C4 F9" or entry["formula_alphabetical"] == "C4 F9 O3 S1":
                NBO_entries.append(entry)
                #age_list.append("new")
                age_list.append("old")
                nov_sp += 1
        else:
            NBO_entries.append(entry)
            #age_list.append("new")
            age_list.append("old")
            nov_sp += 1
    else:
        print("Not connected!")
        print(entry["formula_alphabetical"])
        print(entry["output"]["initial_molecule"]["charge"])
        #print(Molecule.from_dict(entry["output"]["initial_molecule"]))
print(nov_sp, "Nov_SP entries")
pfbs = 0
for entry in pfbs_entries:
    mg = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(entry["output"]["initial_molecule"]), OpenBabelNN())
    if nx.is_connected(mg.graph.to_undirected()):
        NBO_entries.append(entry)
        #age_list.append("new")
        age_list.append("old")
        pfbs += 1
    else:
        print("Not connected!")
        print(entry["formula_alphabetical"])
        print(entry["output"]["initial_molecule"]["charge"])
print(pfbs, "PFBS entries")
cnbz = 0
for entry in cnbz_entries:
    mg = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(entry["output"]["initial_molecule"]), OpenBabelNN())
    if mg.molecule.charge > 0:
        continue
    if entry["task_label"] == "FFopt_frag4_C0_SP_NBO_tzvp":
        continue
    #mg.molecule.to("xyz", "/global/home/groups/lr_mp/smblau/EUVL/dec_recomb/check_cnbz_xyzs/"+str(cnbz)+".xyz")
    #print(entry["task_label"])
    #print(mg.molecule)
    if nx.is_connected(mg.graph.to_undirected()):
        NBO_entries.append(entry)
        age_list.append("new")
        cnbz += 1
    else:
        print("Not connected!")
        print(entry["formula_alphabetical"])
        print(entry["output"]["initial_molecule"]["charge"])
print(cnbz, "4cnbz entries")
#print(huh)
print("Done", len(NBO_entries), "NBO docs")

print()
print("=====")
print()

print("Building a molecule graph for each doc and tagging hot atoms via NBO spins for open shell fragments")

mg_list = []
for entry in NBO_entries:
    mg = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(entry["output"]["initial_molecule"]), OpenBabelNN())
    if entry["output"]["initial_molecule"]["spin_multiplicity"] == 2:
        nbo_spins = entry["output"]["nbo"]["natural_populations"][0]["Density"]
    else:
        nbo_spins = []
    tagged_mg = add_hot_atom_tags(mg, nbo_spins) #if nbo spins is empty and this function is called, just adds an empty list as the "hot atoms"
    mg_list.append(tagged_mg)
print("Done", len(mg_list), "tagged mol_graphs") #all fragments in this list are radicals with "radical rich" sites

print()
print("=====")
print()

print("Building and appending hydrogen molecule graphs")

H_site = Site("H", [0.0, 0.0, 0.0])
H0_mol = Molecule.from_sites([H_site])
H0_mg = MoleculeGraph.with_local_env_strategy(H0_mol, OpenBabelNN())
Hminus_mol = Molecule.from_sites([H_site])
Hminus_mol.set_charge_and_spin(-1)
Hminus_mg = MoleculeGraph.with_local_env_strategy(Hminus_mol, OpenBabelNN())
Hplus_mol = Molecule.from_sites([H_site])
Hplus_mol.set_charge_and_spin(1)
Hplus_mg = MoleculeGraph.with_local_env_strategy(Hplus_mol, OpenBabelNN())
mg_list.append(add_hot_atom_tags(H0_mg, []))
mg_list.append(add_hot_atom_tags(Hminus_mg, []))
mg_list.append(add_hot_atom_tags(Hplus_mg, []))
age_list.append("old")
age_list.append("old")
age_list.append("old")
print("Done")

# print()
# print("=====")
# print()

# # print("Generating recombinations")
# # recomb = generate_combinations(mg_list, mg_list, "/global/home/groups/lr_mp/smblau/EUVL/dec_recomb/", age_list)
# # print("Done", len(recomb), "recombinations")

# print()
# print("=====")
# print()

# print("Generating recombinant structures")
# recomb_mgs = []
# combos = parse_combinations_file("/global/home/groups/lr_mp/smblau/EUVL/dec_recomb/dec_combinations.txt") #This function returns a list of tuples
# for ii,reaction in enumerate(combos): #pairs each tuple in combos with every other tuple in combos
#     print(ii)
#     mg1 = mg_list[reaction[0]] #mg_list contains all species valid for recombination?
#     atom1 = reaction[1] #considering the call for the combine_molecules_molassembler function I can only assume this "atom" is an index corresponding to a reactive atom
#     mg2 = mg_list[reaction[2]]
#     atom2 = reaction[3]
#     new_mg = combine_molecules_molassembler(mg1, mg2, atom1, atom2)
#     recomb_mgs.append(new_mg)
# num_recomb_mgs = 0
# for ii, mg in enumerate(recomb_mgs):
#     if mg is not None:
#         print(mg.molecule.charge)
#         num_recomb_mgs += 1
#         mg.molecule.to("json", "/global/home/groups/lr_mp/smblau/EUVL/dec_recomb/dec_recomb_jsons/"+str(ii)+".json") #outputs new recombinants to a json file for later DFT
#         mg.molecule.to("xyz", "/global/home/groups/lr_mp/smblau/EUVL/dec_recomb/dec_recomb_xyzs/"+str(ii)+".xyz")
# print("Done", num_recomb_mgs, "recombinant jsons written")

hash_list = [] #we'll use these hashes later to ensure we aren't adding any duplicate fragments
for entry in mg_list:
    mg_hash = weisfeiler_lehman_graph_hash(entry.graph.to_undirected(), node_attr="specie") 
    hash_list.append(mg_hash)
    
num_db_recombs = 0
double_bond_recombs = gen_double_bonded_products(mg_list)
 
for i, mol_fragment_list in enumerate(double_bond_recombs.values()):
    for recomb in mol_fragment_list:
        for fragment in recomb:
            print("number added: ", num_db_recombs)
            frag_hash = weisfeiler_lehman_graph_hash(fragment.graph.to_undirected(), node_attr="specie")
            match = False #this test makes sure we're not adding any redundant recombinants
            for entry in hash_list: 
                if entry == frag_hash:
                    match = True
                    break
            if match:
                continue
            elif not match: #can probably get rid of the above if statement
                fragment.molecule.set_charge_and_spin(0.0)
                two_neighbor_carbon_found = False #searches for carbons with only two neighborns
                lone_fragment_test = False #searches for single atom fragments
                num_atoms = len(fragment.molecule)
                if num_atoms == 1:
                    lone_fragment_test = True
                for j, atom in enumerate(fragment.molecule):
                    if str(fragment.molecule[int(j)].specie) == "C":
                        connected_sites = fragment.get_connected_sites(j)
                        if len(connected_sites) <= 2:
                             two_neighbor_carbon_found = True
                if not two_neighbor_carbon_found and not lone_fragment_test:
                    if fragment.molecule.spin_multiplicity == 1: #removes open-shell recombinants--we are looking for closed-shell species
                        num_db_recombs += 1
                        fragment.molecule.to("/global/home/users/jrmilton/"+str(num_db_recombs)+".json", "json") #outputs new recombinants to a json file for later DFT
                        fragment.molecule.to("/global/home/users/jrmilton/"+str(num_db_recombs)+".xyz", "xyz")
                        hash_list.append(frag_hash)
    
print("Done", num_db_recombs, "recombinant jsons written")
