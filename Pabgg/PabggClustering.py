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

#################################
# should not have to change below
##################################

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


density_names = {"PgA1":["PgA1"],"PgA2":["PgA2"],"PgA3":["PgA3"],"PgA4":["PgA4"],"PgA5":["PgA5"],"PgB1":["PgB1"],"PgB2":["PgB2"],"PgB3":["PgB3"],"PgB4":["PgB4"],"PgB5":["PgB5"],
               "PgA6":["PgA6"],"PgB6":["PgB6"] }


rmsd_names = {"PgA1":"PgA1","PgA2":"PgA2","PgA3":"PgA3","PgA4":"PgA4","PgA5":"PgA5","PgB1":"PgB1","PgB2":"PgB2","PgB3":"PgB3","PgB4":"PgB4","PgB5":"PgB5",
              "PgA6":"PgA6","PgB6":"PgB6"}


align_names ={"PgA1":"PgA1","PgA2":"PgA2","PgA3":"PgA3","PgA4":"PgA4","PgA5":"PgA5","PgB1":"PgB1","PgB2":"PgB2","PgB3":"PgB3","PgB4":"PgB4","PgB5":"PgB5",
              "PgA6":"PgA6","PgB6":"PgB6"} 

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