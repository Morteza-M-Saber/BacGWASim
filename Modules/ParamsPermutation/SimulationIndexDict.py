
## Parameters required to be defined for each model 
#SubstitutionModel(name:string, parameters:list,frequencies:list,neutralDNA:bolean)
Global_models={'substModel_codon':{'name':['substModel_codon']},\
               'substModel_TN93':{'name':['substModel_TN93'],'parameters':['substModel_TN93_alpha1','substModel_TN93_alpha2','substModel_TN93_beta'],'frequencies':['substModel_TN93_frequencies'],'neutral_DNA':['substModel_TN93_neutral_DNA']},
               'substModel_HKY':{'name':['substModel_HKY'],'parameters':['substModel_HKY_alpha','substModel_HKY_beta'],'frequencies':['substModel_HKY_frequencies'],'neutral_DNA':['substModel_HKY_neutral_DNA']},
               'substModel_F84':{'name':['substModel_F84'],'parameters':['substModel_F84_kappa','substModel_F84_beta'],'frequencies':['substModel_F84_frequencies'],'neutral_DNA':['substModel_F84_neutral_DNA']},
               'substModel_GTR':{'name':['substModel_GTR'],'parameters':['substModel_GTR_alpha','substModel_GTR_beta','substModel_GTR_gamma','substModel_GTR_delta','substModel_GTR_epsilon','substModel_GTR_zeta'],'frequencies':['substModel_GTR_frequencies'],'neutral_DNA':['substModel_GTR_neutral_DNA']},         

               #IndelModel(Gainrate:nonnegative,model:string,parametersLlist,maxlen:posint,Lossrate:nonnegative)
               'indel_model_ZIPF':{'Gainrate':['indel_gainRate_ZIPF'], 'model':['indel_model_ZIPF'],'parameters':['indel_gainParameters_ZIPF'],'maxLen':['indel_gainMaxLen_ZIPF'],'Lossrate':['indel_lossRate_ZIPF']},
               'indel_model_QG':{'Gainrate':['indel_gainRate_QG'],'model':['indel_model_QG'],'parameters':['indel_gainParameters_QG'],'maxLen':['indel_gainMaxLen_QG'],'Lossrate':['indel_lossRate_QG']},
               'indel_model_GEOM':{'Gainrate':['indel_gainRate_GEOM'],'model':['indel_model_GEOM'],'parameters':['indel_gainParameters_GEOM'],'maxLen':['indel_gainMaxLen_GEOM'],'Lossrate':['indel_lossRate_GEOM']},
               'indel_model_Custom':{'Gainrate':['indel_gainRate_Custom'],'model':['indel_model_Custom'],'parameters':['indel_gainParameters_Custom'],'maxLen':['indel_gainMaxLen_Custom']},
               
               #RateVarModel(model:string,areas:posint,motifFreq:nonnegative,alpha:nonnegative)
               'rateVar_model_Gamma':{'model':['rateVar_model_Gamma'],'areas':['rateVar_model_Gamma_areas'],'motifFreq':['rateVar_model_Gamma_motifFreq'],'alpha':['rateVar_model_Gamma_alpha']},
               'rateVar_model_poisson':{'model':['rateVar_model_poisson'],'areas':['rateVar_model_poisson_areas'],'motifFreq':['rateVar_model_poisson_motifFreq']}}

Global_models_noncoding={'substModel_GTR':{'name':['substModel_GTR_nc'],'parameters':['substModel_GTR_alpha_nc','substModel_GTR_beta_nc','substModel_GTR_gamma_nc','substModel_GTR_delta_nc','substModel_GTR_epsilon_nc','substModel_GTR_zeta_nc'],'frequencies':['substModel_GTR_frequencies_nc']},
                         'substModel_K2P':{'name':['substModel_K2P_nc'],'parameters':['substModel_K2P_transition_nc','substModel_K2P_transversion_nc'],'frequencies':['substModel_K2P_frequencies_nc']},
                         'substModel_K3P':{'name':['substModel_K3P_nc'],'parameters':['substModel_K3P_alpha_nc','substModel_K3P_beta','substModel_K3P_gamma_nc'],'frequencies':['substModel_K3P_frequencies_nc']},
                         'substModel_HKY':{'name':['substModel_HKY_nc'],'parameters':['substModel_HKY_alpha_nc','substModel_HKY_beta_nc'],'frequencies':['substModel_HKY_frequencies_nc']},
                         'substModel_F84':{'name':['substModel_F84_nc'],'parameters':['substModel_F84_kappa_nc'],'frequencies':['substModel_F84_frequencies_nc']},
                         'substModel_TN':{'name':['substModel_TN93_nc'],'parameters':['substModel_TN93_alpha1_nc','substModel_TN93_alpha2_nc','substModel_TN93_beta_nc'],'frequencies':['substModel_TN93_frequencies_nc']},
                         
                         'indel_model_GEOM':{'model':['indel_model_GEOM_nc'],'Gainrate':['indel_gainRate_GEOM_nc'],'parameters':['indel_gainParameters_GEOM_nc'],'maxLen':['indel_gainMaxLen_GEOM_nc']},
                         'indel_model_Power_law':{'model':['indel_model_PL_nc'],'Gainrate':['indel_gainRate_PL_nc'],'parameters':['indel_gainParameters_PL_nc'],'maxLen':['indel_gainMaxLen_PL_nc']},
                         'indel_model_US':{'model':['indel_model_US_nc'],'Gainrate':['indel_gainRate_US_nc'],'parameters':['indel_gainParameters_US_nc']},
                         
                         'rateVar_model_Gamma':{'model':['rateVar_model_Gamma_nc'],'parameters':['rateVar_model_Gamma_alpha_nc','rateVar_model_Gamma_invariant_site_ratio_nc']}}
                         
                         
                         

user_input=config['user_input']
##Reading Models and parameters defined by the user     
with open (user_input,'r') as file:
    line=file.readline()
    while line:
        if line.split(':')[0] == 'Protein_coding_genes':
            Protein_coding_genes_models=line.split(':')[1].strip().split(',')
        elif line.split(':')[0] == 'Intergenic_regions':
            Intergenic_regions_models=line.split(':')[1].strip().split(',')
        elif line.split(':')[0] == 'RNA_genes':
            RNA_genes_models=line.split(':')[1].strip().split(',')        
        elif line.split(':')[0] == 'Repeat_regions':
            Repeat_regions_models=line.split(':')[1].strip().split(',')
        line=file.readline()

##Parameter ranges defined in 'Garage.py', 
#important: The Garage.py should be repalced by yaml file importing parameters automatically 
import os
#dirname, filename = os.path.split(os.path.abspath(sys.argv[0]))
#garage= os.path.join(dirname, "garage.py")

ModuleDir=config['ModuleDir']
garage=os.path.join(ModuleDir,'ParamsPermutation/garage.py')
            
from garage import *           
         
#Inserting required dictionary into a dictionary 
permute_dict={}
for param_vars in Protein_coding_genes_models:
    for protein_model_sub in Global_models[param_vars]:
        for sub_protein_model_sub in Global_models[param_vars][protein_model_sub]:
                permute_dict[sub_protein_model_sub+'_prot']= range_dict[sub_protein_model_sub][0]

for param_vars in Intergenic_regions_models:
    for intergenic_model_sub in Global_models_noncoding[param_vars]:
        for sub_intergenic_model_sub in Global_models_noncoding[param_vars][intergenic_model_sub]:
                permute_dict[sub_intergenic_model_sub+'_intergenic']= range_dict[sub_intergenic_model_sub][0]
                
for param_vars in RNA_genes_models:
    for RNA_genes_sub in Global_models_noncoding[param_vars]:
        for sub_RNA_genes_model_sub in Global_models_noncoding[param_vars][RNA_genes_sub]:
                permute_dict[sub_RNA_genes_model_sub+'_RNA_genes']= range_dict[sub_RNA_genes_model_sub][1]
                
for param_vars in Repeat_regions_models:
    for Repeat_regions_sub in Global_models_noncoding[param_vars]:
        for sub_Repeat_regions_model_sub in Global_models_noncoding[param_vars][Repeat_regions_sub]:
                permute_dict[sub_Repeat_regions_model_sub+'_repeat_regions']= range_dict[sub_Repeat_regions_model_sub][2]

permute_dict['amongGeneDistr']= range_dict['amongGeneDistr']
permute_dict['aGAlpha']= range_dict['aGAlpha']

permute_dict['geneDuplRate']= range_dict['geneDuplRate']
permute_dict['numberDupl']= range_dict['numberDupl']
permute_dict['transDupl']= range_dict['transDupl']
permute_dict['fissionDupl']= range_dict['fissionDupl']
permute_dict['fusionDupl']= range_dict['fusionDupl']
permute_dict['P_pseudogene']= range_dict['P_pseudogene']
permute_dict['P_neofunc']= range_dict['P_neofunc']
permute_dict['P_subfunc']= range_dict['P_subfunc']

permute_dict['geneLossRate']= range_dict['geneLossRate']
permute_dict['numberLoss']= range_dict['numberLoss']

permute_dict['lgtRate']= range_dict['lgtRate']
permute_dict['orthRep']= range_dict['orthRep']
permute_dict['lgtGRate']= range_dict['lgtGRate']
permute_dict['lgtGSize']= range_dict['lgtGSize']

permute_dict['invers']= range_dict['invers']
permute_dict['invSize']= range_dict['invSize']
permute_dict['transloc']= range_dict['transloc']
permute_dict['transSize']= range_dict['transSize']
permute_dict['invtrans']= range_dict['invtrans']

permute_dict['fissionRate']= range_dict['fissionRate']
permute_dict['fusionRate']= range_dict['fusionRate']
permute_dict['numberFusion']= range_dict['numberFusion']


from itertools import product
permute_list_iteration=[dict(zip(permute_dict, v)) for v in product(*permute_dict.values())]



simulation_index_dict={}
for item in permute_list_iteration:
        simulation_index_dict[permute_list_iteration.index(item)]=item
