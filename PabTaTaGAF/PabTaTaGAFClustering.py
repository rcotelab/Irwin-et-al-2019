import matplotlib as mpl
mpl.use('Agg')

import IMP
import IMP.pmi
import IMP.pmi.macros
import sys,os


num_clusters = 1
num_top_models = 100
merge_directories = [""]
prefiltervalue = 2900.0
out_dir = "kmeans_%i_%i/" %(num_top_models,num_clusters)
if '--test' in sys.argv: prefiltervalue=8000.0



model=IMP.Model()


mc=IMP.pmi.macros.AnalysisReplicaExchange0(model,
                                           merge_directories=merge_directories)


feature_list=["ISDCrossLinkMS_Distance_intrarb",
              "ISDCrossLinkMS_Distance_interrb",
              "ISDCrossLinkMS_Data_Score",
              "GaussianEMRestraint_None",
              "SimplifiedModel_Linker_Score_None",
              "ISDCrossLinkMS_Psi",
              "ISDCrossLinkMS_Sigma"]


density_names = {"TaA":["TaA"],
               "Pa":["Pa"],"Pb":["Pb"],"TaB":["TaB"] }


rmsd_names = {"Ta":"Ta",
              "Pa":"Pa","Pb":"Pb","TaB":"TaB"}


align_names ={"Ta":"Ta",
              "Pa":"Pa","Pb":"Pb","TaB":"TaB"} 

mc.clustering(prefiltervalue=prefiltervalue,                   
              number_of_best_scoring_models=num_top_models,    
              alignment_components=None,                       
              rmsd_calculation_components=rmsd_names,          
              distance_matrix_file="distance.rawmatrix.pkl",   
              outputdir=out_dir,                              
              feature_keys=feature_list,                       
              load_distance_matrix_file=False,                 
              display_plot=True,                               
              exit_after_display=False,                        
              get_every=1,                                     
              number_of_clusters=num_clusters,                 
              voxel_size=3.0,                                  
              density_custom_ranges=density_names)             