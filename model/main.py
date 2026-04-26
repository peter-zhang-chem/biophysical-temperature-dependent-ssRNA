# written by Peter Zhang and Hung Nguyen
# TODO This code simulate RNA with explicit ions

from openmm import app
import openmm as omm
from openmm import unit
from openmm.app import *
from openmm import *
import numpy as np
import time
import sys
import argparse
import build
import force

KELVIN_TO_KT = unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB / unit.kilocalorie_per_mole

parser = argparse.ArgumentParser(description='Coarse-grained SOP_IDP simulation using OpenMM')
parser.add_argument('-f','--sequence', type=str, help='input sequence')
parser.add_argument('-p','--pdb', type=str, help='input pdb')
parser.add_argument('-pc','--pdb_coordinates', type=str, help='input pdb coordinates') # * used for benchmarking against C code
parser.add_argument('-c','--cutoff', type=float, default='40.',
                    help='Cutoff distance for electrostatics (A) [40.0]')
parser.add_argument('-T','--temperature', type=float, default='20.',
                    help='Temperature (oC) [20.0]')
parser.add_argument('-t','--traj', type=str, default='md.dcd',
                    help='trajectory output')
parser.add_argument('-e','--energy', type=str, default='energy.out',
                    help='energy decomposition')
parser.add_argument('-o','--output', type=str, default='md.out',
                    help='status and energy output')
parser.add_argument('-r','--res_file', type=str, default="checkpnt.chk",
                    help='checkpoint file for restart')
parser.add_argument('-R','--resume', action="store_true",
                    help='flag to resume a simulation from checkpoint')
parser.add_argument('-fr','--from_res_file', type=str, default="checkpnt.chk",
                    help='checkpoint file for keep running the simulation')
parser.add_argument('-d', '--append_dcd', type=str, help='append to this existing dcd file')
parser.add_argument('-x','--frequency', type=int, default='10000',
                    help='output frequency')
parser.add_argument('-s','--step', type=int, default='10000',
                    help='Number of step [10000]')
parser.add_argument('-K','--monovalent_concentration', type=float, default='150.',
                    help='Monovalent concentration (mM) [150.0]')
parser.add_argument('-M','--divalent_concentration', type=float, default='1.0',
                    help='Divalent concentration (mM) [1.0]')
parser.add_argument('-v', '--box_size', type=float, default='80.',
                    help='Box length (A) [80.0]')
parser.add_argument('-n','--pdb_name', type=str, default='input.pdb',
                    help='PDB output file name')
parser.add_argument('-ts', '--time_step', type=int, default='10',
                    help='time step(fs/step)')
parser.add_argument('-info', '--simu_info', type=str, default="system.out",
                    help="simulation input saved in a text file")
args = parser.parse_args()

class simu:    ### structure to group all simulation parameter
    temp = 0.
    Kconc = 0.
    Mconc = 0.
    Nstep = 0
    epsilon = 0.
    cutoff = 40. * unit.angstrom
    b = 4.38178046 * unit.angstrom / unit.elementary_charge
    b_unitless = b / (unit.angstrom / unit.elementary_charge)
    ion_q = 2
    box_size = 0
    box_size_unitless = 0
    time_step = 10
    resume = False

simu.resume = args.resume
simu.temp = (args.temperature + 273.15) * unit.kelvin
#simu.temp = args.temperature * unit.kelvin
simu.Nstep = args.step
simu.Kconc = args.monovalent_concentration
simu.Mconc = args.divalent_concentration
simu.box_size = args.box_size * unit.angstrom
simu.box_size_unitless = simu.box_size/unit.angstrom
simu.cutoff = args.cutoff * unit.angstrom
simu.time_step = args.time_step

T_unitless = simu.temp * KELVIN_TO_KT
simu.epsilon = 296.0736276 - 619.2813716 * T_unitless + 531.2826741 * T_unitless**2 - 180.0369914 * T_unitless**3
simu.l_Bjerrum = 332.0637*unit.angstroms / simu.epsilon
simu.kappa = unit.sqrt (4* np.pi * simu.l_Bjerrum * 2*simu.Kconc*6.022e-7 / (T_unitless * unit.angstrom**3))
simu.Q = simu.b * T_unitless * unit.elementary_charge**2 / simu.l_Bjerrum
simu.Q_unitless = simu.Q/unit.elementary_charge
N_Mg = int(args.divalent_concentration * 6.022e-7 * pow(args.box_size, 3))

with open(args.simu_info, 'w') as f:
    f.write(f"Box size = {simu.box_size} Å\n")
    f.write(f"Electrostatic cut off distance = {simu.cutoff} Å\n")
    f.write(f"Number of steps = {simu.Nstep}\n")
    f.write(f"Monovalent concentration = {simu.Kconc} mM\n")
    f.write(f"Divalent concentration = {simu.Mconc} mM\n")
    f.write(f"Time step = {simu.time_step} fs\n")
    f.write(f"Number of Mg2+ added = {N_Mg}\n")
    f.write(f"Phosphate charge = {-simu.Q}\n")
    f.write(f"kappa = {simu.kappa}\n")
    f.write(f"Bjerrum length = {simu.l_Bjerrum / T_unitless}\n")
    f.write(f"epsilon = {simu.epsilon}\n")
    f.write(f"T_unitless = {T_unitless}\n")
    f.write(f"Q_unitless = {simu.Q_unitless}\n")

forcefield = app.ForceField('ForceField.xml')
topology = None
positions = None

if args.sequence != None:
    print("Building from sequence %s ..." % args.sequence)
    topology, positions = build.build_by_seq(args.sequence, N_Mg, simu.box_size_unitless, forcefield)
elif args.pdb != None:
    print("Building from pdb file ... %s" % args.pdb)
    topology, positions = build.build_by_pdb(args.pdb, N_Mg, simu.box_size_unitless, forcefield)
elif args.pdb_coordinates != None:
    topology, positions = build.build_by_coordinates(args.pdb_coordinates, N_Mg, forcefield)
else:
    print("Need sequence !!!")
    sys.exit()

#~~~~~~~~~~~~~~ for constraint RNA simulations~~~~~~~~~~~~~~~~#
#system = forcefield.createSystem(topology, constraints=AllBonds)
#for atom in topology.atoms():
#    if atom.residue.name != 'Mg':
#        system.setParticleMass(atom.index, 0 * unit.amu)
#        print(f"{atom} mass was set to 0")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
system = forcefield.createSystem(topology)
system.setDefaultPeriodicBoxVectors(omm.Vec3(simu.box_size, 0, 0), omm.Vec3(0, simu.box_size, 0), omm.Vec3(0, 0, simu.box_size))

#~~~~~~~~~~~~~~~ add force ~~~~~~~~~~~~~~~~~~~~~~~~#
force.add_bond_force (topology, system, 0)
force.add_angle_force (topology, system, 1)
force.add_DH_force   (topology, system, 0, simu, 2)
force.add_WCA_force (topology, system, 3)
force.add_Stack_force(topology, system, 0, simu, 4)
force.add_MP_Potential(topology, system, 0, simu, 5)
totalforcegroup = 6
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class EnergyReporter(object):
    def __init__ (self, file, reportInterval):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval

    def __del__ (self):
        self._out.close()

    def describeNextReport(self, simulation):
        step = self._reportInterval - simulation.currentStep%self._reportInterval
        return (step, False, False, False, True)
        #return (step, position, velocity, force, energy)

    def report(self, simulation, state):
        energy = []
        self._out.write(str(simulation.currentStep))
        for i in range(totalforcegroup):
            state = simulation.context.getState(getEnergy=True, getForces=True, groups=2**i)
            energy = state.getPotentialEnergy() / unit.kilocalorie_per_mole
            self._out.write("   " + str(energy))
        self._out.write("\n")

integrator = omm.LangevinMiddleIntegrator(simu.temp, 0.01/unit.picosecond, simu.time_step*unit.femtoseconds)
platform = omm.Platform.getPlatformByName('CPU')
#platform = omm.Platform.getPlatformByName('CUDA')
#properties = {'CudaPrecision': 'mixed'}
#properties["DeviceIndex"] = "0,1"
simulation = app.Simulation(topology, system, integrator)
#simulation = app.Simulation(topology, system, integrator, platform, properties)

#~~~~~~~~~~~~~~~~~~~~~~~To run simulation from checkpoint~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# remove the simulation.context.setPositions(positions) line
# remove DCD reporter
# remove energy minimization line
# make a copy of the original checkpoint file then read that new file to avoid losing data
# uncomment the lines below 
#with open(args.from_res_file, 'rb') as f:
#   simulation.context.loadCheckpoint(f.read())
#   simulation.reporters.append(app.DCDReporter(args.append_dcd, args.frequency, append=True))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if (simu.resume == False):
    simulation.context.setPositions(positions)
    print("Initial energy   %f   kcal/mol" % (simulation.context.getState(getEnergy=True).getPotentialEnergy() / unit.kilocalorie_per_mole))
    state = simulation.context.getState(getPositions=True)
    app.PDBFile.writeFile(topology, state.getPositions(), open(args.pdb_name, "w"), keepIds=True)
    for i in range(totalforcegroup):
                state = simulation.context.getState(getEnergy=True, getForces=True, groups=2**i)
                energy = state.getPotentialEnergy() / unit.kilocalorie_per_mole
                if i == 0:
                    print("Bond Force = ", energy)
                elif i == 1:
                    print("Angle Force = ", energy)
                elif i == 2:
                    print("DH Force = ", energy)
                elif i == 3:
                    print("WCA force = ", energy)
                elif i == 4:
                    print("Stack force = ", energy)
                elif i == 5:
                    print("MP potential = ", energy)
    print("OpenMM version = ", platform.getOpenMMVersion())
    if args.pdb_coordinates == None:
        print('Minimizing ...')
        simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(simu.temp)
    simulation.reporters.append(app.DCDReporter(args.traj, args.frequency))

elif (simu.resume == True):
    with open(args.from_res_file, 'rb') as f:
        simulation.context.loadCheckpoint(f.read())
        simulation.reporters.append(app.DCDReporter(args.traj, args.frequency, append=False))

simulation.reporters.append(app.StateDataReporter(args.output, args.frequency, step=True, potentialEnergy=True, temperature=True, remainingTime=True, totalSteps=simu.Nstep, separator='  '))
simulation.reporters.append(EnergyReporter(args.energy, args.frequency))
simulation.reporters.append(app.CheckpointReporter(args.res_file, args.frequency)) # update the Checkpoint reporter file the same frequency as the trajectory file

#~~~~~~~~~~~~~~For openMM folks to trouble shoot GPU issue~~~~~~~~~~~~~~#
#simulation.saveState('output.xml') # to save state of the simulation

#with open('integrator.xml', 'w') as output:
#    output.write(XmlSerializer.serialize(integrator))

#with open('system.xml', 'w') as output:
#    output.write(XmlSerializer.serialize(system))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

print('Running simulation ...')
t0 = time.time()
simulation.currentStep = 0
simulation.step(simu.Nstep)
prodtime = time.time() - t0
secs_in_day = 60.*60.*24
simulation_speed_steps = secs_in_day*simu.Nstep/(prodtime) 
print("Simulation speed: % .2e steps/day" % (simulation_speed_steps))
step_size_fs_unitless = integrator.getStepSize()/unit.femtoseconds
step_size_ns_unitless = step_size_fs_unitless/1000000
simulation_speed_ns = simulation_speed_steps * step_size_ns_unitless
print("Simulation speed: % .2e ns/day" % (simulation_speed_ns))
