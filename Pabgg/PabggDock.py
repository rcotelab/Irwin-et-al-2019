"""
#############################################
## P gamma PDE Docking
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
topology_file = datadirectory+"PabggTopology.txt"


#--------------------------
# Set MC Sampling Parameters
#--------------------------
num_frames = 20000
if '--test' in sys.argv: num_frames=100
num_mc_steps = 10

#--------------------------
# Create movers
#--------------------------


rb_max_trans = 4.00
rb_max_rot = 0.04


bead_max_trans = 4.00


rigid_bodies = [["PgA1"],
                ["PgA2"],["PgA3"],["PgA4"],["PgA5"],["PgA6"],["PgB1"],["PgB2"],["PgB3"],["PgB4"],["PgB5"],["PgB6"]]
super_rigid_bodies = [["PgA1"],
                ["PgA2"],["PgA3"],["PgA4"],["PgA5"],["PgA6"],["PgB1"],["PgB2"],["PgB3"],["PgB4"],["PgB5"],["PgB6"]]
chain_of_super_rigid_bodies = [["PgA1"],
                ["PgA2"],["PgA3"],["PgA4"],["PgA5"],["PgA6"],["PgB1"],["PgB2"],["PgB3"],["PgB4"],["PgB5"],["PgB6"]]
#
################################################
#
def add_excluded_volume(self,prot,kappa):
        m=prot.get_model()
        rs = IMP.RestraintSet(m,'excluded_volume')
        atoms=IMP.atom.get_by_type(prot, IMP.atom.ATOM_TYPE)
        for atom in atoms:
            restype=IMP.atom.Residue(IMP.atom.Atom(atom).get_parent()).get_residue_type()
            vol=IMP.atom.get_volume_from_residue_type(restype)
            radius=IMP.algebra.get_ball_radius_from_volume_3d(vol)
            IMP.core.XYZR(atom).set_radius(radius)
        lsa=IMP.container.ListSingletonContainer(self.m)
        lsa.add(atoms)
        evr=IMP.core.ExcludedVolumeRestraint(lsa,kappa)
        rs.add_restraint(evr)
        self.all_restraints.add_restraint(rs)
        self.rest_set[('exvo',prot)]=rs
def add_external_barrier(self,rad,prot):
        m=prot.get_model()
        rs = IMP.RestraintSet(m,'barrier')
        c3= IMP.algebra.Vector3D(0,0,0)
        ub3= IMP.core.HarmonicUpperBound(rad, 10.0)
        ss3= IMP.core.DistanceToSingletonScore(ub3, c3)
        lsc= IMP.container.ListSingletonContainer(self.m)
        IMP.atom.get_by_type

        lsc.add(IMP.atom.get_leaves(prot))
        r3= IMP.container.SingletonsRestraint(ss3, lsc)
        rs.add_restraint(r3)
        self.all_restraints.add_restraint(rs)
        self.rest_set[('barr',prot)]=rs		

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


# Crosslinks - dataset 1

columnmap={}
columnmap["Protein1"]="prot1"
columnmap["Protein2"]="prot2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None

xl1 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representation,
                                   datadirectory+'PABggBS3.csv',
                                   length=30.0,
                                   slope=0.6,
                                   columnmapping=columnmap,
                                   resolution=1.0,
                                   label="BS3",
                                   csvfile=True)

xl1.add_to_model()             

sampleobjects.append(xl1) 
outputobjects.append(xl1)


# Crosslinks - dataset 2

columnmap={}
columnmap["Protein1"]="prot1"
columnmap["Protein2"]="prot2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None

xl2 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representation,
                                   datadirectory+'PABggMBS.csv',
                                   length=20.0,
                                   slope=0.7,
                                   columnmapping=columnmap,
                                   resolution=1.0,
                                   label="Peptide",
                                   csvfile=True)
xl2.add_to_model()

sampleobjects.append(xl2)
outputobjects.append(xl2)






#Crosslinks - dataset 5
columnmap={}
columnmap["Protein1"]="prot1"
columnmap["Protein2"]="prot2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None

xl5 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representation,
                                   datadirectory+'PABggEDC.csv',
                                   length=15.0,
                                   slope=0.5,
                                   columnmapping=columnmap,
                                   resolution=1.0,
                                   label="Peptide",
                                   csvfile=True)
xl5.add_to_model()

sampleobjects.append(xl5)
outputobjects.append(xl5)


#Crosslinks - dataset 6
columnmap={}
columnmap["Protein1"]="prot1"
columnmap["Protein2"]="prot2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None



#Crosslinks - dataset 8
columnmap={}
columnmap["Protein1"]="prot1"
columnmap["Protein2"]="prot2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None

xl8 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representation,
                                   datadirectory+'PABggpeg9.csv',
                                   length=50.0,
                                   slope=0.01,
                                   columnmapping=columnmap,
                                   resolution=1.0,
                                   label="Peptide",
                                   csvfile=True)
xl8.add_to_model()

sampleobjects.append(xl8)
outputobjects.append(xl8)

#Crosslinks - dataset 7
columnmap={}
columnmap["Protein1"]="prot1"
columnmap["Protein2"]="prot2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None




representation.optimize_floppy_bodies(10)
####


mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                    representation,
                                    monte_carlo_sample_objects=sampleobjects,
                                    output_objects=outputobjects,
                                    crosslink_restraints=[xl1,xl2,xl5,xl8],    
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


mc1.execute_macro()
