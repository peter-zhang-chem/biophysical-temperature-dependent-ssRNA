# written by Peter Zhang and Hung Nguyen
# This program builds coarse-grained RNA structure with three interaction sites (ribose sugar, phosphate, and nucleobase) either by input sequence or input pdb file.

from openmm import app #simtk.openmm is deprecated
from openmm import unit
import openmm as omm
from Bio.PDB import PDBParser
import numpy as np
import random


atoms = []                  # Initializing the list storing all the atoms
topo = omm.app.Topology()   # Initializing topology <Topology: 0 chains, 0 residues, 0 atoms, 0 bonds>

# ** build coarse-grained model of RNA sequence ** #
def build_by_seq(seq, N_Mg, box_size, forcefield):
    
    # ** Defining geometry of each nucleotide on the x,y plane ** #
    # ! openMM has unit in nanometer! * #
    # Ribose sugar is placed at the origin (0, 0, 0)
    # geometry = {'Nucleotide': {'Ribose Sugar (S)': [x, y, z], 'Phosphate (P)': [x, y, z], 'Base (A, G, C, or U)': [x, y, z]}
    # S-P distance was obtained from line 819 in RNA-cg-master (4.6010 Å), x-value of P was calulated with 4.6010/sec(π - θ), y-value of P was calculated with 4.6010/csc(π - θ) 
    
    center = (box_size/2)*0.1
    geometry = {'ADE': {'S': [center, center, center], 'P': [center-0.15997, center+0.43139, center], 'A': [center+0.48515, center, center]}, # 'A' x-value = 4.8515 (line 727), θ(A-S-P) = 1.9259 (line 735) obtained from RNA-cg-master(force.h).
                'GUA': {'S': [0, 0, 0], 'P': [-0.15526, 0.43311, 0], 'G': [0.49659, 0, 0]}, # 'G' x-value = 4.9659 (line 759), θ(G-S-P) = 1.9150 (line 765) obtained from RNA-cg-master(force.h).
                'CYT': {'S': [0, 0, 0], 'P': [-0.17692, 0.42472, 0], 'C': [0.42738, 0, 0]}, # 'C' x-value = 4.2738 (line 744), θ(C-S-P) = 1.9655(line 750) obtained from RNA-cg-master(force.h).
                'URA': {'S': [0, 0, 0], 'P': [-0.17726, 0.42458, 0], 'U': [0.42733, 0, 0]}, # 'U' x-value = 4.2733 (line 774), θ(U-S-P) = 1.9663(line 780) obtained from RNA-cg-master(force.h).
                'ADE5': {'S': [center, center, center], 'P': [center-0.15997, center-0.43139, center], 'A': [center+0.48515, center, center]}, # 5' ADE
                'GUA5': {'S': [0, 0, 0], 'P': [-0.15526, -0.43311, 0], 'G': [0.49659, 0, 0]}, # 5' GUA
                'CYT5': {'S': [0, 0, 0], 'P': [-0.17692, -0.42472, 0], 'C': [0.42738, 0, 0]}, # 5' CYT
                'URA5': {'S': [0, 0, 0], 'P': [-0.17726, -0.42458, 0], 'U': [0.42733, 0, 0]}, # 5' URA
                'ADE3': {'S': [center, center, center], 'P': [center-0.15997, center+0.43139, center], 'P3': [center-0.15997, center-0.43139, center], 'A': [center+0.48515, center, center]}, # 3' ADE
                'GUA3': {'S': [0, 0, 0], 'P': [-0.15526, 0.43311, 0], 'P3': [-0.15526, -0.43311, 0], 'G': [0.49659, 0, 0]}, # 3' GUA
                'CYT3': {'S': [0, 0, 0], 'P': [-0.17692, 0.42472, 0], 'P3': [-0.17692, -0.42472, 0], 'C': [0.42738, 0, 0]}, # 3' CYT
                'URA3': {'S': [0, 0, 0], 'P': [-0.17726, 0.42458, 0], 'P3': [-0.17726, -0.42458, 0], 'U': [0.42733, 0, 0]}} # 3' URA
    # Name map couples identifier to one base ('ADE' -> 'A'), ('GUA' -> 'G'), ('CYT' -> 'C'), ('URA' -> 'U')
    name_map = {'A': 'ADE', 'G': 'GUA', 'C': 'CYT', 'U': 'URA'} 
    step = 0.63581 # distance between adjacent sugar, the angle of S-P-S is 1.5256 rad (line 828) obtained from RNA-cg-master(force.h).
    
    # ** Constuct the initial RNA structure on the x,y-plane ** #
                    # P             #     
                    #   \           #
                    #    S --- Base #
                    #   /           #
                    # P             #
                    #   \           #    
                    #    S --- Base #
                    #   /           #    
                    # P             #
    
    # ** delcare initial variables ** #
    chain = topo.addChain('A')  # <Topology: 1 chains, 0 residues, 0 atoms, 0 bonds>
    ion_chain = topo.addChain('B')
    positions = []              # Initial the list to store all the atom positions (same as the ones listed in geometry above)
    idx_offset = 0              # idx_offset is updated with (len(forcefield._templates[symbol].atoms)) on line 80
    
    # ** adding all the atoms to atoms[] ** #
    for i, resSymbol in enumerate(seq):     # getting the index of nucleotide and its one letter code
        symbol = name_map[resSymbol]        # assign symbol to the three letter code from name_map ('ADE' --> 'A', 'GUA' --> 'G', 'CYT' --> 'C', or 'URA' --> 'U')
        if i == len(seq) - 1:               # if it is at the 3' terminal, the symbol will be ('ADE3', 'GUA3', 'CYT3', or 'URA3')
            symbol = symbol + "3"
        elif i == 0:
            symbol = symbol + "5"           # if it is at the 5' terminal, the symbol would be ('ADE5', 'GUA5', 'CYT5', or 'URA5')                
        res_geometry = geometry[symbol]     # getting the 'S', 'P', 'A', 'P3' (for 3' phosine) coordinates of the corresponding nucleotide
        res = topo.addResidue(symbol, chain) # output: <Residue # (symbol) of chain #)
        for atom in forcefield._templates[symbol].atoms:
            atoms.append(topo.addAtom(atom.name, forcefield._atomTypes[atom.type].element, res)) # addAtom(name, element, residue[, id]) -> to create a new atom and add it to the topology
            if atom.name in res_geometry:
                positions.append(res_geometry[atom.name])   # populate positions[] with all the atom's position, the list will be updated at the end with a new list "updated_positions"
            else:
                print ("Residue %s not found!!" % atom.name)
                exit()
        
        # ** Assigning all the internal bonds ** #
        for bond in forcefield._templates[symbol].bonds:        # refer to the bonds section in ForceField.xml
            #print ("Adding bond  %s - %s" % (atoms[bond[0] + idx_offset], atoms[bond[1] + idx_offset]))
            topo.addBond(atoms[bond[0] + idx_offset], atoms[bond[1] + idx_offset])  # forming internal bond

        # ** Assigning all the external bonds (bonds between nucleotides) ** #
        p_idx = None #p_idx tracks all the phosphates
        s_idx = None #s_idx tracks all the sugar
        for bond in forcefield._templates[symbol].externalBonds:
            p_idx = bond + idx_offset # update p_idx (make sure the forcefield is in order)
        if p_idx > 2:
            s_idx = p_idx - 4 #the external bond between sugar is always with phosphate 4 atoms ahead
            #print ("Adding external bond  to %s  %s" % (atoms[p_idx], atoms[s_idx]))
            topo.addBond(atoms[p_idx], atoms[s_idx])    # add to openMM topology, returned at the end
        idx_offset += len(forcefield._templates[symbol].atoms)  #idx_offset moves by the number of atoms assigned to a symbol (ex. idx_offset of GUA3 would be 4).
    
    # ** update all the coordinates along the xy-axis ** #
    updated_positions = []  # initializing the list that stores updated coordinate
    for i in range (len(atoms)):
        start_Base = positions[0][1] # starting coordinates of first base
        start_Sugar = positions[1][1]  # starting coordinates of first sugar
        start_Phos = positions[2][1]   # starting coordinates of first phosphate
        cur_pos = positions[i].copy()      # variable cur_pos (current position) is updated by increasing i
        if i >= 3:  # the first three atoms do not needs to be updated 
            if i%3 == 0 and i != len(atoms) - 1:    # index of all the base will give i%3 = 0, and the 3' phosphate is an exception; hence the and statement.
                cur_pos[1] = start_Base + (i//3*step)    # update the coordinate of base
            elif i%3 == 1:  # update all the sugar
                cur_pos[1] = start_Sugar + (i//3*step)
            elif i%3 == 2:  # update all the phosphate except the 3' phosphate
                cur_pos[1] = start_Phos + (i//3*step)
            elif i == len(atoms) - 1:   # update the 3' phosphate
                cur_pos[1] = start_Phos + (len(seq)*step)
        updated_positions.append(cur_pos) # populate the update_positions list
    
        # ** adding all the divalent ions to topology ** #
    for Mg in range (N_Mg): 
        symbol = "Mg"
        res = topo.addResidue(symbol, ion_chain)
        topo.addAtom("Mg", forcefield._templates["Mg"].atoms[0].element, res) # topo.addAtom(atom.name, forcefield._atomTypes[atom.type].element, res))        
        
        ion_position = [random.uniform(0, box_size * 0.1) for _ in range(3)]
        #print("ion_position = ", ion_position)
        updated_positions.append(ion_position)

    return topo, updated_positions

def build_by_pdb(pdb, N_Mg, box_size, forcefield):
    return_positions = []
    # Name map couples identifier to one base ('ADE' -> 'A'), ('GUA' -> 'G'), ('CYT' -> 'C'), ('URA' -> 'U')
    parser = PDBParser(QUIET=True)
    name_map = {'A': 'ADE', 'G': 'GUA', 'C': 'CYT', 'U': 'URA'}
    structure = parser.get_structure("RNA", pdb)
    seq = []
    for model in structure:
        for chain in model:
            for residue in chain:
                res_symbol = residue.get_resname()
                seq.append(res_symbol)
                #print("res_symbol = ", res_symbol)
                # get the coarse grained coordinates
                phos_coords = []
                sugar_coords = []
                base_coords = []
                phos3_coords = []
                for base in residue:
                    atom_name = base.get_name()
                    if atom_name in ["N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "C2", "C4", "C5", "C6", "C8", "O2", "O4", "O6"]:
                        base_coords.append(list(base.get_vector()))
                if(len(base_coords) != 0):
                    return_positions.append(list(np.mean(base_coords, axis = 0) * 0.1))
                
                for sugar in residue:
                    atom_name = sugar.get_name()
                    if atom_name in ["C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O3'", "O4'", "O5'"]:
                        #print("Atom name = %s, phos.get_vector() = %s " % (atom_name, sugar.get_vector()))
                        sugar_coords.append(list(sugar.get_vector()))
                if(len(sugar_coords) != 0):
                    return_positions.append(list(np.mean(sugar_coords, axis = 0) * 0.1))
                
                for phos in residue:
                    atom_name = phos.get_name()
                    if atom_name in ['P', "OP1", "OP2"]:
                        #print("Atom name = %s, phos.get_vector() = %s " % (atom_name, phos.get_vector()))
                        phos_coords.append(list(phos.get_vector()))
                if(len(phos_coords) != 0):
                    return_positions.append(list(np.mean(phos_coords, axis = 0) * 0.1 ))
                
                for phos3 in residue:
                    atom_name = phos3.get_name()
                    if atom_name in ["P3", "OP3", "OP3E"]:
                        phos3_coords.append(list(phos3.get_vector()))
                if(len(phos3_coords) != 0):
                    return_positions.append(list(np.mean(phos3_coords, axis = 0) * 0.1))
    
                #print("Base coords = %s, length = %s" % (base_coords, len(base_coords)))
                #print("Sugar coords = %s, length = %s" % (sugar_coords, len(sugar_coords)))
                #print("Phosphate coords = %s, length = %s" % (phos_coords, len(phos_coords)))
                #print("Phos 3 coords = %s, length = %s" % (phos3_coords, len(phos3_coords)))

    chain = topo.addChain('A')
    ion_chain = topo.addChain('B')
    idx_offset = 0

    for i, resSymbol in enumerate(seq):
        symbol = name_map[resSymbol]
        if i == len(seq) - 1:               # if it is at the 3' terminal, the symbol will be ('ADE3', 'GUA3', 'CYT3', or 'URA3')
            symbol = symbol + "3"
        elif i == 0:
            symbol = symbol + "5"           # if it is at the 5' terminal, the symbol would be ('ADE5', 'GUA5', 'CYT5', or 'URA5')                
        res = topo.addResidue(symbol, chain) # output: <Residue # (symbol) of chain #)
        for atom in forcefield._templates[symbol].atoms:
            atoms.append(topo.addAtom(atom.name, forcefield._atomTypes[atom.type].element, res)) # addAtom(name, element, residue[, id]) -> to create a new atom and add it to the topology

    # ** Assigning all the internal bonds ** #
        for bond in forcefield._templates[symbol].bonds:        # refer to the bonds section in ForceField.xml
            #print ("Adding bond  %s - %s" % (atoms[bond[0] + idx_offset], atoms[bond[1] + idx_offset]))
            topo.addBond(atoms[bond[0] + idx_offset], atoms[bond[1] + idx_offset])  # forming internal bond
        # ** Assigning all the external bonds (bonds between nucleotides) ** #
        p_idx = None #p_idx tracks all the phosphates
        s_idx = None #s_idx tracks all the sugar
        
        for bond in forcefield._templates[symbol].externalBonds:
            p_idx = bond + idx_offset # update p_idx (make sure the forcefield is in order)
        if p_idx > 2:
            s_idx = p_idx - 4 #the external bond between sugar is always with phosphate 4 atoms ahead
            #print ("Adding external bond  to %s  %s" % (atoms[p_idx], atoms[s_idx]))
            topo.addBond(atoms[p_idx], atoms[s_idx])    # add to openMM topology, returned at the end
        idx_offset += len(forcefield._templates[symbol].atoms)  #idx_offset moves by the number of atoms assigned to a symbol (ex. idx_offset of GUA3 would be 4).
    
    print(N_Mg)
    for Mg in range (N_Mg): 
        symbol = "Mg"
        res = topo.addResidue(symbol, ion_chain)
        topo.addAtom("Mg", forcefield._templates["Mg"].atoms[0].element, res) # topo.addAtom(atom.name, forcefield._atomTypes[atom.type].element, res))        
        
    
    for to_add in range (0, N_Mg):
        #print("cutoff = ", box_size)
        ion_position = [random.uniform(0, box_size * 0.1) for _ in range(3)] # * here we genearate random Mg coordinates, box_size * 0.1 for unit to be in nm (the input is in angstrom).
        return_positions.append(ion_position)        
    
    print("sequence = ", seq)
    print("return positions = ", len(return_positions))
    print(topo)
    return topo, return_positions

def build_by_coordinates(coordinates, N_Mg, box_size, forcefield):
    return_positions = []
    # Name map couples identifier to one base ('ADE' -> 'A'), ('GUA' -> 'G'), ('CYT' -> 'C'), ('URA' -> 'U')
    parser = PDBParser(QUIET=True)
    name_map = {'A': 'ADE', 'G': 'GUA', 'C': 'CYT', 'U': 'URA'}
    structure = parser.get_structure("RNA", coordinates)
    seq = []
    Mg_count = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                res_symbol = residue.get_resname()
                #print("res_symbol = ", res_symbol)
                if res_symbol != "M":
                    seq.append(res_symbol)
                for atom in residue:
                    atom_name = atom.get_name()
                    #print(atom_name)
                    if atom_name in ['A', 'C', 'U', 'G']:
                        return_positions.append(np.divide(list(atom.get_vector()),10))
                    if atom_name == 'S':
                        return_positions.append(np.divide(list(atom.get_vector()),10))
                    if atom_name == 'P':
                        return_positions.append(np.divide(list(atom.get_vector()),10))
                    if atom_name == 'P3':
                        return_positions.append(np.divide(list(atom.get_vector()),10))
                    #elif atom_name == 'Mg':
                        #Mg_count +=1
                        #return_positions.append(np.divide(list(atom.get_vector()),10))
    print("sequence = ", seq)
    #print("return positions = ", return_positions)
    chain = topo.addChain('A')
    ion_chain = topo.addChain('B')
    idx_offset = 0

    for i, resSymbol in enumerate(seq):
        symbol = name_map[resSymbol]
        if i == len(seq) - 1:               # if it is at the 3' terminal, the symbol will be ('ADE3', 'GUA3', 'CYT3', or 'URA3')
            symbol = symbol + "3"
        elif i == 0:
            symbol = symbol + "5"           # if it is at the 5' terminal, the symbol would be ('ADE5', 'GUA5', 'CYT5', or 'URA5')                
        res = topo.addResidue(symbol, chain) # output: <Residue # (symbol) of chain #)
        for atom in forcefield._templates[symbol].atoms:
            atoms.append(topo.addAtom(atom.name, forcefield._atomTypes[atom.type].element, res)) # addAtom(name, element, residue[, id]) -> to create a new atom and add it to the topology
    # ** Assigning all the internal bonds ** #
        for bond in forcefield._templates[symbol].bonds:        # refer to the bonds section in ForceField.xml
            #print ("Adding bond  %s - %s" % (atoms[bond[0] + idx_offset], atoms[bond[1] + idx_offset]))
            topo.addBond(atoms[bond[0] + idx_offset], atoms[bond[1] + idx_offset])  # forming internal bond
        # ** Assigning all the external bonds (bonds between nucleotides) ** #
        p_idx = None #p_idx tracks all the phosphates
        s_idx = None #s_idx tracks all the sugar
        
        for bond in forcefield._templates[symbol].externalBonds:
            p_idx = bond + idx_offset # update p_idx (make sure the forcefield is in order)
        if p_idx > 2:
            s_idx = p_idx - 4 #the external bond between sugar is always with phosphate 4 atoms ahead
            #print ("Adding external bond  to %s  %s" % (atoms[p_idx], atoms[s_idx]))
            topo.addBond(atoms[p_idx], atoms[s_idx])    # add to openMM topology, returned at the end
        idx_offset += len(forcefield._templates[symbol].atoms)  #idx_offset moves by the number of atoms assigned to a symbol (ex. idx_offset of GUA3 would be 4).

    print(N_Mg)
    for Mg in range (N_Mg): 
        symbol = "Mg"
        res = topo.addResidue(symbol, ion_chain)
        topo.addAtom("Mg", forcefield._templates["Mg"].atoms[0].element, res)
        ion_position = [random.uniform(0, box_size * 0.1) for _ in range(3)] # * here we genearate random Mg coordinates, box_size * 0.1 for unit to be in nm (the input is in angstrom).
        return_positions.append(ion_position)


    # for explicit ion positions
    #ion_position = [[1.6845, 3.6339, 4.0680], [4.4230, 5.2953, 1.1042]]
    #ion_position = [[2.1093, 3.3062, 1.9455], [4.6942, 1.4761, 5.5380], [5.0061, 4.7129, 3.9030], [4.0429, 5.0551, 0.4129], [2.0214, 1.1954, 3.8172], [4.6097, 2.4112, 1.4029],
    #                [0.9750, 0.0886, 0.6183], [5.8688, 3.4875, 5.5773], [2.8083, 3.8483, 4.4700], [1.8854, 2.2283, 5.8770], [3.9062, 0.0623, 2.2428], [0.1669, 3.7190, 2.5461],
    #                [5.9001, 0.5336, 4.1782]]

    #for to_add in range (0, N_Mg):
    #    return_positions.append(ion_position[to_add])
    
    print(topo)
    print(len(return_positions))

    return topo, return_positions

                    
