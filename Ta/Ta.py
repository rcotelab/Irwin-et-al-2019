"""
#############################################
##  Transducin Refinement
##
#############################################

"""
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.topology

import os
import sys


#---------------------------
# Define Input Files
#---------------------------
datadirectory = ""
topology_file = datadirectory+"topologyTa.txt"


#--------------------------
# Set MC Sampling Parameters
#--------------------------
num_frames = 20000
if '--test' in sys.argv: num_frames=100
num_mc_steps = 10

#--------------------------
# Create movers
#--------------------------


rb_max_trans = 20.00
rb_max_rot = 1.0


bead_max_trans = 20.00


rigid_bodies = [["Ta"],["ntah"]]
super_rigid_bodies = [["Ta"],["ntah"],["Ta","ntah"]]
chain_of_super_rigid_bodies = [["Ta","ntah"]]
#
################################################
#


m = IMP.Model()


topology = IMP.pmi.topology.TopologyReader(topology_file)
domains = topology.component_list


bm = IMP.pmi.macros.BuildModel(m,
                    component_topologies=domains,
                    list_of_rigid_bodies=rigid_bodies,
                    list_of_super_rigid_bodies=super_rigid_bodies,
                    chain_of_super_rigid_bodies=chain_of_super_rigid_bodies)
representation = bm.get_representation()


for nc,component in enumerate(domains):
    name = component.name
    sel = IMP.atom.Selection(representation.prot,molecule=name)
    ps = sel.get_selected_particles()
    clr = IMP.display.get_rgb_color(float(nc)/len(domains))
    for p in ps:
        if not IMP.display.Colored.get_is_setup(p):
            IMP.display.Colored.setup_particle(p,clr)
        else:
            IMP.display.Colored(p).set_color(clr)



representation.shuffle_configuration(50)




representation.set_rigid_bodies_max_rot(rb_max_rot)
representation.set_floppy_bodies_max_trans(bead_max_trans)
representation.set_rigid_bodies_max_trans(rb_max_trans)

outputobjects = [] 
sampleobjects = [] 


outputobjects.append(representation)
sampleobjects.append(representation)


ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         representation, resolution=20)
ev.add_to_model()
outputobjects.append(ev)



columnmap={}
columnmap["Protein1"]="pep1"
columnmap["Protein2"]="pep2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None

xl1 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representation,
                                   datadirectory+'TarefineBS3.csv',
                                   length=26.0,
                                   slope=0.01,
                                   columnmapping=columnmap,
                                   resolution=1.0,
                                   label="BS3",
                                   csvfile=True )

xl1.add_to_model()             

sampleobjects.append(xl1) #
outputobjects.append(xl1)


# Crosslinks - dataset 2
#
columnmap={}
columnmap["Protein1"]="prot1"
columnmap["Protein2"]="prot2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None

xl2 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representation,
                                   datadirectory+'TarefineSDA.csv',
                                   length=19.0,
                                   slope=0.01,
                                   columnmapping=columnmap,
                                   resolution=1.0,
                                   label="SDA",
                                  csvfile=True )

xl2.add_to_model()             

sampleobjects.append(xl2) 
outputobjects.append(xl2)

#Crosslink dataset 3
columnmap={}
columnmap["Protein1"]="prot1"
columnmap["Protein2"]="prot2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None

xl3 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representation,
                                   datadirectory+'TarefineEDC.csv',
                                   length=13.0,
                                   slope=0.01,
                                   columnmapping=columnmap,
                                   resolution=1.0,
                                   label="EDC",
                                   csvfile=True)

xl3.add_to_model()             

sampleobjects.append(xl3) 
outputobjects.append(xl3)


representation.optimize_floppy_bodies(10)
####

# 
#--------------------------
# Monte-Carlo Sampling
#--------------------------

mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    representation,
                                    monte_carlo_sample_objects=sampleobjects,
                                    output_objects=outputobjects,
                                    crosslink_restraints=[xl1,xl2,xl3],    
                                    monte_carlo_temperature=1.0,
                                    simulated_annealing=True,
                                    simulated_annealing_minimum_temperature=1.0,
                                    simulated_annealing_maximum_temperature=2.5,
                                    simulated_annealing_minimum_temperature_nframes=200,
                                    simulated_annealing_maximum_temperature_nframes=20,
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=2.5,
                                    number_of_best_scoring_models=100,
                                    monte_carlo_steps=num_mc_steps,
                                    number_of_frames=num_frames,
                                    global_output_directory="output")
force_restraint=True


mc1.execute_macro()
