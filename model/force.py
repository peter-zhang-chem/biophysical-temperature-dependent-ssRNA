# written by Peter Zhang and Hung Nguyen
# This program add interaction forces of the structure built from build.py

import openmm as omm
import numpy as np
from openmm import unit
from openmm import app
from build import atoms, topo


# * bond force * #
# ** E = kappa(p-p0)**2 (force constant was multiplied by 2 as OpenMM harmonic energy is represented by (1/2)k(angle - equilibrium angle)^2#
# ** the bond force need to be considered are the following: ** #
# ** A --- S (l = 4.8515, k = 10.0) ** #
# ** C --- S (l = 4.2738, k = 10.0) ** #
# ** G --- S (l = 4.9659, k = 10.0) ** #
# ** U --- S (l = 4.2733, k = 10.0) ** #
# ** P --- S (start with P end with S (from the previous residue) (l = 3.8157, k = 64.0)) ** #
# ** S --- P (start with S end with P (from the same residue) (l = 4.6010, k = 23.0)) ** #
def add_bond_force(topology, system, forcegroup):
    bondforce = omm.HarmonicBondForce()
    kappa_B_S = 2*10.0000*(unit.kilocalorie_per_mole/unit.angstroms**2) # * kappa_B_S (Base to Sugar kappa value)
    kappa_S_P = 2*64.0000*(unit.kilocalorie_per_mole/unit.angstroms**2) # * Sugar to Phosphate downstream
    kappa_P_S = 2*23.0000*(unit.kilocalorie_per_mole/unit.angstroms**2) # * Phosphate to sugar downstream
    i = 0
    for bond in topo.bonds():
        bond_atom_1_name = atoms[bond[0].index].name
        bond_atom_2_name = atoms[bond[1].index].name
        bond_atom_1 = bond[0].index
        bond_atom_2 = bond[1].index
        if bond_atom_1_name == 'A' and bond_atom_2_name == 'S':
            bondforce.addBond(bond_atom_1, bond_atom_2, 4.8515*unit.angstroms, kappa_B_S)
            #print("Bond added between %s and %s" % (bond_atom_1, bond_atom_2))
        elif bond_atom_1_name == 'C' and bond_atom_2_name == 'S':
            bondforce.addBond(bond_atom_1, bond_atom_2, 4.2738*unit.angstroms, kappa_B_S)
        elif bond_atom_1_name == 'G' and bond_atom_2_name == 'S':
            bondforce.addBond(bond_atom_1, bond_atom_2, 4.9659*unit.angstroms, kappa_B_S)
            #print("Bond added between %s and %s" % (bond_atom_1, bond_atom_2))
        elif bond_atom_1_name == 'U' and bond_atom_2_name == 'S':
            bondforce.addBond(bond_atom_1, bond_atom_2, 4.2733*unit.angstroms, kappa_B_S)
        elif bond_atom_1_name == 'P' and bond_atom_2_name == 'S':
            bondforce.addBond(bond_atom_1, bond_atom_2, 3.8157*unit.angstroms, kappa_S_P)
            #print("Bond added between %s and %s" % (bond_atom_1, bond_atom_2))
        elif bond_atom_1_name == 'S' and (bond_atom_2_name == 'P'):
            bondforce.addBond(bond_atom_1, bond_atom_2, 4.6010*unit.angstroms, kappa_P_S)
            #print("Bond added between %s and %s" % (bond_atom_1, bond_atom_2))
        elif bond_atom_1_name == 'S' and bond_atom_2_name == 'P3':
            bondforce.addBond(bond_atom_1, bond_atom_2, 3.8157*unit.angstroms, kappa_S_P) # * from S -> P downstream
            #print("Bond added between %s and %s" % (bond_atom_1, bond_atom_2))
        i += 1
    bondforce.setUsesPeriodicBoundaryConditions(True)
    bondforce.setForceGroup(forcegroup)
    system.addForce(bondforce)

# * angle force * #
# ** angle needs to be considered are the following: ** #
    # * base-sugar-phosphate angle force
          # ** Type             Bond                            Parameters            **#
            # ** 1) A --- S --- P (of the next residue) (theta = 1.9259, kappa = 5.0) ** #
            # ** 2) C --- S --- P (of the next residue) (theta = 1.9655, kappa = 5.0) ** #
            # ** 3) G --- S --- P (of the next residue) (theta = 1.9150, kappa = 5.0) ** #
            # ** 4) U --- S --- P (of the next resiude) (theta = 1.9663, kappa = 5.0) ** #
            # ** 5) A --- S --- P (theta = 1.7029, kappa = 5.0) ** #
            # ** 6) C --- S --- P (theta = 1.5803, kappa = 5.0) ** #
            # ** 7) G --- S --- P (theta = 1.7690, kappa = 5.0) ** #
            # ** 8) U --- S --- P (theta = 1.5735, kappa = 5.0) ** #
            # ** 9) P --- S --- P (S and P of the previous residue) (theta = 1.4440, kappa = 20.0) ** #
            # ** 10) S --- P --- S (of the previous residue) (theta = 1.5256, kappa = 20.0) ** #

angleforce = omm.HarmonicAngleForce()

def add_angle_force(topology, system, forcegroup):
    kappa_B_S_P = 2*5.0*unit.kilocalorie_per_mole/(unit.radians**2) # ** kappa_B_S_P(Base-Sugar-Phosphate+1 kappa value)
    kappa_backBone = 2*20.0*unit.kilocalorie_per_mole/(unit.radians**2) # ** kappa_backBone (P-S-P or S-P-S kappa value)
    for base in range(len(atoms) - 1):  # ** loop through the first base to the phosphate before the 3' phosphate
        # ** adding all the angle force that involves base ** #
        if atoms[base].name != 'S' or 'P':      # * exclude index starting with S or P
            if atoms[base].name == 'A':         # * if the base is A
                if (base == len(atoms) - 4):    # * if the base is the 3' base
                    angleforce.addAngle(base, base + 1, base + 2, 1.7029*unit.radian, kappa_B_S_P) # * Type 5 bond (see chart above)
                    #print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[base + 2].name)
                    angleforce.addAngle(base, base + 1, len(atoms) - 1, 1.9259*unit.radian, kappa_B_S_P) # * Type 1 bond  with the 3' phosphate
                    #print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[len(atoms) - 1].name)
                else: # * if the base is not at the 3'
                    angleforce.addAngle(base, base + 1, base + 2, 1.7029*unit.radian, kappa_B_S_P) # * Type 5 bond
                    #print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[base + 2].name)
                    angleforce.addAngle(base, base + 1, base + 5, 1.9259*unit.radian, kappa_B_S_P) # * Type 3 bond
                    #print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[base + 5].name)
            elif atoms[base].name == 'C':   # * if the base is C
                if (base == len(atoms) - 4):    # * if the base is the 3' base
                    angleforce.addAngle(base, base + 1, base + 2, 1.5803*unit.radian, kappa_B_S_P) # * Type 2 bond
                    # print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[base + 2].name)
                    angleforce.addAngle(base, base + 1, len(atoms) - 1, 1.9655*unit.radian, kappa_B_S_P) # * Type 2 bond with the 3' phosphate
                    # print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[len(atoms) - 1].name)
                else: # * if the base is not at the 3'
                    angleforce.addAngle(base, base + 1, base + 2, 1.5803*unit.radian, kappa_B_S_P) 
                    # print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[base + 2].name)
                    angleforce.addAngle(base, base + 1, base + 5, 1.9655*unit.radian, kappa_B_S_P)
                    # print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[base + 5].name)
            elif atoms[base].name == 'G':   # * if the base is G
                if (base == len(atoms) - 4): # * if the base is the 3' base
                    angleforce.addAngle(base, base + 1, base + 2, 1.7690*unit.radian, kappa_B_S_P)
                    #print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[base + 2].name)
                    angleforce.addAngle(base, base + 1, len(atoms) - 1, 1.9150*unit.radian, kappa_B_S_P)
                    #print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[len(atoms) - 1].name)
                else: # * if the base is not at the 3'
                    angleforce.addAngle(base, base + 1, base + 2, 1.7690*unit.radian, kappa_B_S_P)
                    #print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[base + 2].name)
                    angleforce.addAngle(base, base + 1, base + 5, 1.9150*unit.radian, kappa_B_S_P)
                    #print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[base + 5].name)
            elif atoms[base].name == 'U':   # * if the base is U
                if (base == len(atoms) - 4):    # * if the base is the 3' base
                    angleforce.addAngle(base, base + 1, base + 2, 1.5735*unit.radian, kappa_B_S_P)
                    # print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[base + 2].name)
                    angleforce.addAngle(base, base + 1, len(atoms) - 1, 1.9663*unit.radian, kappa_B_S_P)
                    # print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[len(atoms) - 1].name)
                else: # * if the base is not at the 3'
                    angleforce.addAngle(base, base + 1, base + 2, 1.5735*unit.radian, kappa_B_S_P)
                    # print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[base + 2].name)
                    angleforce.addAngle(base, base + 1, base + 5, 1.9663*unit.radian, kappa_B_S_P)
                    # print("angleforce added between ", atoms[base].name, atoms[base + 1].name, atoms[base + 5].name)
        # ** adding all the angle force of backbones
        if atoms[base].name == 'P':     # * add the angle force of P --- S --- P
            if (base == len(atoms) - 2):    # * at the phosphate before 3' Phosphate
                angleforce.addAngle (base, base - 1, base + 1, 1.4440*unit.radian, kappa_backBone)
                #print("angleforce added between ", atoms[base].name, atoms[base - 1].name, atoms[base + 1].name)
            else:   # * not the 3' Phosphate
                angleforce.addAngle (base, base - 1, base + 3, 1.4440*unit.radian, kappa_backBone)
                #print("angleforce added between ", atoms[base].name, atoms[base - 1].name, atoms[base + 3].name)
        if atoms[base].name == 'S':      # ** add the angle force of S --- P --- S
            if (base != len(atoms) - 3):    # * none 3' Sugar (3' Sugar do not have S --- P --- S angle)
                angleforce.addAngle (base, base + 4, base + 3, 1.5256*unit.radian, kappa_backBone)
                #print("angleforce added between ", atoms[base].name, atoms[base + 4].name, atoms[base + 3].name)
    angleforce.setUsesPeriodicBoundaryConditions(True)
    angleforce.setForceGroup(forcegroup)                                                                                                                                                                                                                                                             
    system.addForce(angleforce)
    return angleforce

# ** Debye-Huckel force ** #
# ** P3_charge: the charge of 3' phosphate (RNA specific value) ** #
def add_DH_force(topology, system, P3_charge, simu, forcegroup):
    # angleforce = add_angle_force(topology, system, forcegroup)
    DHforce = omm.CustomNonbondedForce("scale * q1 * q2 * exp(-kappa*r)/r")
    DHforce.addGlobalParameter("scale", simu.l_Bjerrum * unit.kilocalorie_per_mole / unit.elementary_charge**2)
    DHforce.addGlobalParameter("kappa", simu.kappa)
    DHforce.addPerParticleParameter('q')
    print("Phosphate charge before updates = ", -simu.Q)
    b3 = simu.b_unitless ** 3
    V2 = 4 * np.pi * np.e * b3 * (1. + simu.ion_q) * (1 / simu.Q_unitless - 1./ simu.ion_q)
    V1 = 4 * np.pi * np.e * b3 * (1 + 1) * (1/simu.Q_unitless - 1)
    C1V1 = simu.Kconc * V1 * 6.022e-07
    #print("C1V1 = ", C1V1)
    C2V2 = simu.Mconc * V2 * 6.022e-07
    #print("C2V2 = ", C2V2)
    theta1 = 0
    if(simu.Mconc > 0):
        if (simu.ion_q == 2):
            theta1 = (np.sqrt(pow(C1V1, 4) + 8 * np.e * C1V1 * C1V1 * C2V2 * (1.0 - simu.Q_unitless)) - C1V1 * C1V1) / (4 * np.e * C2V2)
            simu.Q = (1.0 - theta1) * unit.elementary_charge
    print("Phosphate charge after updates = ", -simu.Q)

    # ** Add elemental charges (P = simu.Q, S = B = 0, P3 charge depends on RNA (P3_charge))
    for atom in topology.atoms():
        if atom.name == 'P':
            DHforce.addParticle([-simu.Q])
            #print("SIMU.Q = ", -simu.Q)
        elif atom.name == 'P3':
            DHforce.addParticle([P3_charge * unit.elementary_charge])
            #print("P3 charge = ", P3_charge)
        elif atom.name in ['A', 'C', 'U', 'G', 'S'] :
            DHforce.addParticle([0 * unit.elementary_charge])
            #print(0 * unit.elementary_charge)
        elif atom.name in ['Mg', 'Ca']:
            DHforce.addParticle([simu.ion_q * unit.elementary_charge])
        # print("particle parameter = ", )
    
    # ** Exclude bonded interactions
    for bond in topology.bonds():
        DHforce.addExclusion(bond[0].index, bond[1].index)
        #print ("Adding exclusion bond  to %s  %s" % (bond[0], bond[1]))
    

    Mg_set1 = set()
    Mg_set2 = set()
    Phos_set1 = set()
    Phos_set2 = set()

    for atom in topology.atoms():
        if(atom.name == 'Mg'):
            Mg_set1.add(atom.index)
            Mg_set2.add(atom.index)
        
        elif(atom.name in ['P', 'P3']):
            Phos_set1.add(atom.index)
            Phos_set2.add(atom.index)

    DHforce.addInteractionGroup(Mg_set1, Mg_set2)    
    DHforce.addInteractionGroup(Phos_set1, Phos_set2)

    DHforce.setCutoffDistance(simu.cutoff)
    DHforce.setForceGroup(forcegroup)
    DHforce.setNonbondedMethod(omm.CustomNonbondedForce.CutoffPeriodic)
    system.addForce(DHforce)

def add_MP_Potential(topology, system, P3_charge, simu, forcegroup):
    MP_potential = omm.CustomNonbondedForce("(W(r) + ((D_H - W(r)) * exp(-a^2/r^2))); D_H=scale * q1 * q2 * exp(-kappa*r)/r")
    MP_potential.addGlobalParameter("scale", simu.l_Bjerrum * unit.kilocalorie_per_mole / unit.elementary_charge**2)
    MP_potential.addGlobalParameter("kappa", simu.kappa)
    MP_potential.addGlobalParameter("a", 5 * unit.angstroms)
    MP_potential.addPerParticleParameter('q')
    
    temp_unitless = int(simu.temp.value_in_unit(unit.kelvin) - 273.15)
    with open(f'pmf_MgP.t{temp_unitless}', 'r') as file:
        tabulated_data = [float(line.strip()) * unit.kilocalorie_per_mole for line in file]

    W_r = omm.Continuous1DFunction(tabulated_data, 1.67500000 * unit.angstroms, 50.000000 * unit.angstroms)
    
    #with open('pmf_Mg_P_1264', 'r') as file:        # * tabulated data here is calculated from reference interaction site model (see Hung's PNAS 2019, file obtained from C code)
        #tabulated_data = [float(line.strip()) * unit.kilocalorie_per_mole for line in file]
        #print("tabulated data type = ", type(tabulated_data))
        #print(tabulated_data)
    #W_r = omm.Continuous1DFunction(tabulated_data, 1.77500000 * unit.angstroms, 50.000000 * unit.angstroms)

    MP_potential.addTabulatedFunction('W', W_r)

    # ** Add elemental charges (P = simu.Q, S = B = 0, P3 charge depends on RNA (P3_charge))
    for atom in topology.atoms():
        if atom.name == 'P':
            MP_potential.addParticle([-simu.Q])
            #print("SIMU.Q = ", -simu.Q)
        elif atom.name == 'P3':
            MP_potential.addParticle([P3_charge * unit.elementary_charge])
            #print("P3 charge = ", P3_charge)
        elif atom.name in ['A', 'C', 'U', 'G', 'S'] :
            MP_potential.addParticle([0 * unit.elementary_charge])
            #print(0 * unit.elementary_charge)
        elif atom.name in ['Mg', 'Ca']:
            MP_potential.addParticle([simu.ion_q * unit.elementary_charge])
        # print("particle parameter = ", )
    
    # ** Exclude bonded interactions
    for bond in topology.bonds():
        MP_potential.addExclusion(bond[0].index, bond[1].index)
        #print ("Adding exclusion bond  to %s  %s" % (bond[0], bond[1]))
    
    Mg_set1 = set()
    Phos_set1 = set()

    for atom in topology.atoms():
        if(atom.name == 'Mg'):
            Mg_set1.add(atom.index)
        elif (atom.name in ['P']):
            Phos_set1.add(atom.index)


    MP_potential.addInteractionGroup(Mg_set1, Phos_set1)
    MP_potential.setCutoffDistance(simu.cutoff)
    MP_potential.setForceGroup(forcegroup)
    MP_potential.setNonbondedMethod(omm.CustomNonbondedForce.CutoffPeriodic)
    system.addForce(MP_potential)

# ** WCA force ** #
def add_WCA_force(topology, system, forcegroup):
   ########################### D0 from PNAS 2019 Hung #####################################################
    D0 = [
  #    A        U        C        G       S        P=P3     Mg2+    Ca2+   
    0.32000, 0.32000, 0.32000, 0.32000, 0.51300, 0.44100, 0.45200, 0.53200,   # A  
    0.32000, 0.32000, 0.32000, 0.32000, 0.50400, 0.43200, 0.44200, 0.52300,   # U   
    0.32000, 0.32000, 0.32000, 0.32000, 0.50400, 0.43200, 0.44300, 0.52300,   # C
    0.32000, 0.32000, 0.32000, 0.32000, 0.53100, 0.45900, 0.47000, 0.55000,   # G
    0.51300, 0.50400, 0.50400, 0.53100, 0.52200, 0.45000, 0.46100, 0.54100,   # S
    0.44100, 0.43200, 0.43200, 0.45900, 0.45000, 0.37800, 0.38900, 0.46900,   # P=P3
    0.45200, 0.44200, 0.44300, 0.47000, 0.46100, 0.38900, 0.40000, 0.48000,   # Mg2+
    0.53200, 0.52300, 0.52300, 0.55000, 0.54100, 0.46900, 0.48000, 0.56000]   # Ca2+

    WCA_cutoff = 10.0 * unit.angstroms

############################ WCA Functional Form (A combination of PNAS 2019 and Nature Chemistry 2021) ############################################
    energy_function = "step(sig - r) * ep * ((R6 - 2) * R6 + 1); R6 = (sig / r)^6; sig = D0(type1, type2);"
    WCA_force = omm.CustomNonbondedForce(energy_function)
    WCA_force.addTabulatedFunction('D0', omm.Discrete2DFunction(8, 8, D0))  # D0 from Hung PNAS 2019
    WCA_force.addGlobalParameter('ep', 1.0 * unit.kilocalories_per_mole)
    WCA_force.addPerParticleParameter('type')
    
    Base_sugar_set = set() 
    Mg_set = set ()
    P_set = set() 
    Ca_set = set()

    for atom in topology.atoms():
       if atom.name == 'A':
           WCA_force.addParticle([0])
           Base_sugar_set.add(atom.index)
       elif atom.name == 'U':
           WCA_force.addParticle([1])
           Base_sugar_set.add(atom.index)
       elif atom.name == 'C':
           WCA_force.addParticle([2])
           Base_sugar_set.add(atom.index)
       elif atom.name == 'G':
           WCA_force.addParticle([3])
           Base_sugar_set.add(atom.index)
       elif atom.name == 'S':
           WCA_force.addParticle([4])
           Base_sugar_set.add(atom.index)
       elif atom.name == 'P' or atom.name == "P3":
           WCA_force.addParticle([5])
           P_set.add(atom.index)
       elif atom.name == 'Mg':
           WCA_force.addParticle([6])
           Mg_set.add(atom.index)
       elif atom.name == 'Ca':
           WCA_force.addParticle([7])
           Ca_set.add(atom.index)


########################## add exclusion forces ##########################
    for bond in topology.bonds():
        WCA_force.addExclusion(bond[0].index, bond[1].index)
    

################## exclude Mg and phosphate interactions ##################
    WCA_force.addInteractionGroup(Base_sugar_set, Mg_set)
    WCA_force.addInteractionGroup(Base_sugar_set, Ca_set)
    WCA_force.addInteractionGroup(Base_sugar_set, Base_sugar_set)
    WCA_force.addInteractionGroup(Base_sugar_set, P_set)
    WCA_force.addInteractionGroup(P_set, P_set)
    WCA_force.addInteractionGroup(Mg_set, Mg_set)
    WCA_force.addInteractionGroup(Mg_set, Ca_set)
    WCA_force.addInteractionGroup(Ca_set, Ca_set)


    WCA_force.setCutoffDistance(WCA_cutoff)
    WCA_force.setForceGroup(forcegroup)
    WCA_force.setNonbondedMethod(omm.CustomNonbondedForce.CutoffPeriodic)
    system.addForce(WCA_force)

# ** Stacking ** #
                #-----Secondary Stacking-----#
                    # P                    #     
                    #   \                  #
                    #    S --- Base---     #
                    #   /            |     #
                    # P              | r0  #
                    #   \            |     #    
                    #    S --- Base---     #
# ** @param ss_D: delta G for the dimer. It correts the h variable ()
def add_Stack_force(topology, system, ss_D, simu, forcegroup):
   
    # * Stacking energy function with correction for dihedral angle * #
    # energyfunction = "U0/(1.0 + kbond*(distance(p6,p7) - r0)^2  + kphi1*(dihedral(p1,p2,p3,p4) - phi10 + pi*(((pi - dihedral(p1,p2,p3,p4) + phi10)/abs(pi - dihedral(p1,p2,p3,p4) + phi10)) -  ((pi + dihedral(p1,p2,p3,p4) - phi10)/abs(pi + dihedral(p1,p2,p3,p4) - phi10)))  )^2  +  kphi2*(dihedral(p2,p3,p4,p5) - phi20 + pi*(((pi - dihedral(p2,p3,p4,p5) + phi20)/abs(pi - dihedral(p2,p3,p4,p5) + phi20)) -  ((pi + dihedral(p2,p3,p4,p5) - phi20)/abs(pi + dihedral(p2,p3,p4,p5) - phi20)))  )^2 )"
    
    # * Stacking energy function
    # * Here we take the absolute value around the maximum and correct by -pi, then to compensate the correction, we 
    energyfunction = "U0/(1.0 + kbond * (distance(p6, p7) - r0)^2 + kphi1 * (abs(dihedral(p1, p2, p3, p4) - pi - phi10) - pi)^2 + kphi2 * (abs(dihedral(p2, p3, p4, p5) + pi - phi20) - pi)^2)"

    # ** Create CustomCompoundBondForce object ** #
    stack_force = omm.CustomCompoundBondForce(7, energyfunction)    
    stack_force.addPerBondParameter("U0")
    stack_force.addPerBondParameter("r0")
    stack_force.addPerBondParameter("phi10")
    stack_force.addPerBondParameter("phi20")

    # ** Add global parameters ** #
    stack_force.addGlobalParameter("kbond", 1.4/unit.angstroms**2)
    stack_force.addGlobalParameter("kphi1", 4.0/unit.radians**2)
    stack_force.addGlobalParameter("kphi2", 4.0/unit.radians**2)
    stack_force.addGlobalParameter("pi", np.pi)
    
    # ** Define U0(h, s, Tm), these variables will be assigned in the loop below ** #
    # ** Define r0, h, s, and Tm** #
    # * units: r0 (Å), h (kcal/mol), Tm (K)
    # * direction: 5' -- 3' (obtained from 2013 J. Physic. Chem (Table 2), pay attention to the arrow direction when finding these values)
    # * r0 values were obtained from the C code (stack, starting at line 338)
    # * A -- A r0 = 4.1806530 Å,  A -- C, r0 = 3.8260185 Å,  A -- G, r0 = 4.4255305 Å,  A -- U, r0 = 3.8260185 Å,
    # *         h = 4.980                  h = 4.970                  h = 5.735                  h = 4.970
    # *         s = -0.309                 s = -0.700                 s = 5.280                  s = -0.289  
    # *        Tm = 299.15 K               Tm = 299.15 K              Tm = 341.15 K              Tm = 299.15 K
    # * C -- A r0 = 4.7010580 Å,  C -- C, r0 = 4.2500910 Å,  C -- G, r0 = 4.9790760 Å,  C -- U, r0 = 4.2273615 Å,
    # *         h = 4.940                  h = 4.700                  h = 5.450                  h = 4.700
    # *         s = -0.309                 s = -1.567                 s = 0.300                  s = -1.567   
    # *        Tm = 299.15                Tm = 286.15                Tm = 315.15                Tm = 286.15
    # * G -- A r0 = 4.0128560 Å,  G -- C, r0 = 3.6784360 Å,  G -- G, r0 = 4.2427250 Å,  G -- U, r0 = 3.6616930 Å,
    # *         h = 5.732                  h = 5.783                  h = 6.198                  h = 5.700
    # *         s = 5.240                  s = 4.000                  s = 7.346                  s = 2.200
    # *        Tm = 341.15                Tm = 343.15                Tm = 366.15                Tm = 338.15
    # * U -- A r0 = 4.7010580 Å,  U -- C  r0 = 4.2679180 Å,  U -- G, r0 = 4.9977560 Å,  U -- U, r0 = 4.2453650 Å.      
    # *         h = 4.963                 h = 4.690                  h = 5.700                  h = 4.100
    # *         s = -0.319                s = -1.567                 s = 2.200                  s = -3.563
    # *        Tm = 299.15                Tm = 286.15                Tm = 338.15                Tm = 252.15
   
    for base in range (3, len(atoms) - 3, 3):
        prev_base = atoms[base - 3].name    # * Get the name of the previous atom in relation to the current atom of interest
        cur_base = atoms[base].name     # * Get the name of the current atom
        if prev_base == 'A':
            if cur_base == 'A':
                r0 = 4.1806530 * unit.angstrom
                h = 4.980 * unit.kilocalories_per_mole
               # h = (-5.194 + ss_D / 0.7093838769) * unit.kilocalories_per_mole
                s = -0.309 
                Tm = 299.15
            elif cur_base == 'C':
                r0 = 3.8260185 * unit.angstrom
                h = 4.970 * unit.kilocalories_per_mole
               #  h = (-5.146 + ss_D / 0.7185955600) * unit.kilocalories_per_mole 
                s = -0.700
                Tm = 299.15
            elif cur_base == 'G':
                r0 = 4.4255305 * unit.angstrom
                h = 5.735 * unit.kilocalories_per_mole 
               #  h = (-5.735 + ss_D / 0.6968019829) * unit.kilocalories_per_mole
                s = 5.280
                Tm = 341.15
            elif cur_base == 'U':
                r0 = 3.8260185 * unit.angstrom
                h = 4.970 * unit.kilocalories_per_mole
               #  h = (-5.146 + ss_D / 0.7185955600) * unit.kilocalories_per_mole
                s = -0.289
                Tm = 299.15
        elif prev_base == 'C':
            if cur_base == 'A':
                r0 = 4.7010580 * unit.angstrom
                h = 4.940 * unit.kilocalories_per_mole
                # h = (-5.163 + ss_D / 0.6847830171) * unit.kilocalories_per_mole
                s = -0.309
                Tm = 299.15
            elif cur_base == 'C':
                r0 = 4.2500910 * unit.angstrom
                h = 4.700 * unit.kilocalories_per_mole
                # h = (-4.873 + ss_D / 0.6991615586) * unit.kilocalories_per_mole
                s = -1.567
                Tm = 286.15
            elif cur_base == 'G':
                r0 = 4.9790760 * unit.angstrom
                h = 5.450 * unit.kilocalories_per_mole
                # h = (-5.482 + ss_D / 0.6816268897) * unit.kilocalories_per_mole
                s = 0.300
                Tm = 315.15
            elif cur_base == 'U':
                r0 = 4.2273615 * unit.angstrom 
                h = 4.700 * unit.kilocalories_per_mole
                # h = (-4.873 + ss_D / 0.6832570771) * unit.kilocalories_per_mole
                s = -1.567 
                Tm = 286.15
        elif prev_base == 'G':
            if cur_base == 'A':
                r0 = 4.0128560 * unit.angstrom
                h = 5.732 * unit.kilocalories_per_mole
                # h = (-5.732 + ss_D / 0.6903176657) * unit.kilocalories_per_mole
                s = 5.240
                Tm = 341.15
            elif cur_base == 'C':
                r0 = 3.6784360 * unit.angstrom
                h = 5.783 * unit.kilocalories_per_mole
                # h = (-5.927 + ss_D / 0.7042060343) * unit.kilocalories_per_mole
                s = 4.000
                Tm = 343.15
            elif cur_base == 'G':
                r0 = 4.2427250 * unit.angstrom
                h = 6.198 * unit.kilocalories_per_mole
                # h = (-6.416 + ss_D / 0.6971421514) * unit.kilocalories_per_mole
                s = 7.346
                Tm = 366.15
            elif cur_base == 'U':
                r0 = 3.6616930 * unit.angstrom
                h = 5.700 * unit.kilocalories_per_mole
                # h = (-5.163 + ss_D / 0.6847830171) * unit.kilocalories_per_mole
                s = 2.200
                Tm = 338.15
        elif prev_base == 'U':
            if cur_base == 'A':
                r0 = 4.7010580 * unit.angstrom
                h = 4.963 * unit.kilocalories_per_mole
                # h = (-5.163 + ss_D / 0.6847830171) * unit.kilocalories_per_mole
                s = -0.319
                Tm = 299.15
            elif cur_base == 'C':
                r0 = 4.2679180 * unit.angstrom
                h = 4.690 * unit.kilocalories_per_mole
                # h = (-4.880 + ss_D / 0.6758595771) * unit.kilocalories_per_mole
                s = -1.567
                Tm = 286.15
            elif cur_base == 'G':
                r0 = 4.9977560 * unit.angstrom
                h = 5.700 * unit.kilocalories_per_mole
                # h = (-5.886 + ss_D / 0.7025528229) * unit.kilocalories_per_mole
                s = 2.200
                Tm = 338.15
            elif cur_base == 'U':
                r0 = 4.2453650 * unit.angstrom
                h = 4.100 * unit.kilocalories_per_mole
                #  h = (-4.267 + ss_D / 0.6686014771) * unit.kilocalories_per_mole
                s = -3.563
                Tm = 252.15
        temperature_difference = simu.temp - Tm * unit.kelvin
        U0 = -h + (0.001985875 * unit.kilocalories_per_mole/unit.kelvin * temperature_difference) * s  # * kB = 0.001985875 kcal/(mol * K)
        phi10 = -2.58684 * unit.radians # * equilibrium dihedarl angle
        phi20 = 3.07135 * unit.radians # * equilibrium dihedral angle
        
        p1 = base - 1  # * phosphate
        #print("p1 = ", atoms[p1])
        p2 = base - 2  # * sugar
        #print("p2 = ", atoms[p2])
        p3 = base + 2  # * phosphate
        #print("p3 = ", atoms[p3])
        p4 = base + 1  # * sugar
        #print("p4 = ", atoms[p4])
        p6 = base - 3  # * base
        #print("p6 = ", atoms[p6])
        p7 = base   
        #print("p7 = ", atoms[p7])
        if (p3 < len(atoms) - 2):
            p5 = base + 5 # * phosphate
            #print("p5 (P) = ", atoms[p5])
        elif (p3 == len(atoms) - 2):
            p5 = len(atoms) - 1      # * 3' phosphate
            #print("p5 (P3) = ", atoms[p5])
        group_add = [p1, p2, p3, p4, p5, p6, p7]
        #print("group_add = ", group_add)
        stack_force.addBond(group_add, [U0, r0, phi10, phi20])
    stack_force.setUsesPeriodicBoundaryConditions(True)
    stack_force.setForceGroup(forcegroup)
    system.addForce(stack_force)

# ** Hydrogen bond ** #


