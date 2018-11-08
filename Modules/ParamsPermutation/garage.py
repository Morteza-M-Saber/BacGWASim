import os, sys, numpy


import os
dirname = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
dirname='/'.join(dirname.split('/')[:-1])

with open(os.path.join(dirname,'ConfigFile.yaml'),'r') as file:
	line=file.readline()
	while line:
		if line.split(':')[0]== 'Params_input':
			Params_input=line.split(':')[1].strip()
		line=file.readline()



range_dict={}
alfsim_input=open (Params_input,'r').readlines()
for line in alfsim_input:
        line=line.split('#')[0]
        # Codon model
        if line.split(':')[0] == 'substModel_codon':
            range_dict['substModel_codon']= [line.split(':')[1].strip().split('-')]
        #TN93
        if line.split(':')[0] == 'substModel_TN93':
           range_dict['substModel_TN93']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'substModel_TN93_alpha1':
            substModel_TN93_alpha1_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_TN93_alpha1']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_TN93_alpha1_pre]
        if line.split(':')[0] == 'substModel_TN93_alpha2':
            substModel_TN93_alpha2_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_TN93_alpha2']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_TN93_alpha2_pre]
        if line.split(':')[0] == 'substModel_TN93_beta':
            substModel_TN93_beta_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_TN93_beta']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_TN93_beta_pre]
        if line.split(':')[0] == 'substModel_TN93_frequencies':
            substModel_TN93_frequencies_pre= [item.replace('{','').replace("}",'') for item in line.split(':')[1].strip().split(',')]
            range_dict['substModel_TN93_frequencies']=[[item] for item in substModel_TN93_frequencies_pre]
        if line.split(':')[0] == 'substModel_TN93_neutral_DNA':
            substModel_TN93_neutral_DNA_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_TN93_neutral_DNA']= [item.split('-') for item in substModel_TN93_neutral_DNA_pre]
        #HKY 
        if line.split(':')[0] == 'substModel_HKY':
           range_dict['substModel_HKY']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'substModel_HKY_alpha':
            substModel_HKY_alpha_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_HKY_alpha']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_HKY_alpha_pre]
        if line.split(':')[0] == 'substModel_HKY_beta':
            substModel_HKY_beta_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_HKY_beta']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_HKY_beta_pre]
        if line.split(':')[0] == 'substModel_HKY_frequencies':
            substModel_HKY_frequencies_pre= [item.replace('{','').replace("}",'') for item in line.split(':')[1].strip().split(',')]
            range_dict['substModel_HKY_frequencies']=[[item] for item in substModel_HKY_frequencies_pre]
        if line.split(':')[0] == 'substModel_HKY_neutral_DNA':
            substModel_HKY_neutral_DNA_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_HKY_neutral_DNA']= [item.split('-') for item in substModel_HKY_neutral_DNA_pre]
        #F84
        if line.split(':')[0] == 'substModel_F84':
           range_dict['substModel_F84']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'substModel_F84_kappa':
            substModel_F84_kappa_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_F84']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_F84_kappa_pre]
        if line.split(':')[0] == 'substModel_F84_beta':
            substModel_F84_beta_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_F84_beta']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_F84_beta_pre]
        if line.split(':')[0] == 'substModel_F84_frequencies':
            substModel_F84_frequencies_pre= [item.replace('{','').replace("}",'') for item in line.split(':')[1].strip().split(',')]
            range_dict['substModel_F84_frequencies']=[[item] for item in substModel_F84_frequencies_pre]
        if line.split(':')[0] == 'substModel_F84_neutral_DNA':
            substModel_F84_neutral_DNA_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_F84_neutral_DNA']= [item.split('-') for item in substModel_F84_neutral_DNA_pre]
        #GTR
        if line.split(':')[0] == 'substModel_GTR':
           range_dict['substModel_GTR']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'substModel_GTR_alpha':
            substModel_GTR_alpha_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_GTR_alpha']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_GTR_alpha_pre]
        if line.split(':')[0] == 'substModel_GTR_beta':
            substModel_GTR_beta_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_GTR_beta']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_GTR_beta_pre]
        if line.split(':')[0] == 'substModel_GTR_gamma':
            substModel_GTR_gamma_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_GTR_gamma']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_GTR_gamma_pre]
        if line.split(':')[0] == 'substModel_GTR_delta':
            substModel_GTR_delta_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_GTR_delta']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_GTR_delta_pre]
        if line.split(':')[0] == 'substModel_GTR_epsilon':
            substModel_GTR_epsilon_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_GTR_epsilon']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_GTR_epsilon_pre]
        if line.split(':')[0] == 'substModel_GTR_zeta':
            substModel_GTR_zeta_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_GTR_zeta']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_GTR_zeta_pre]
        if line.split(':')[0] == 'substModel_GTR_frequencies':
            substModel_GTR_frequencies_pre= [item.replace('{','').replace("}",'') for item in line.split(':')[1].strip().split(',')]
            range_dict['substModel_GTR_frequencies']=[[item] for item in substModel_GTR_frequencies_pre]
        if line.split(':')[0] == 'substModel_GTR_neutral_DNA':
            substModel_GTR_neutral_DNA_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_GTR_neutral_DNA']= [item.split('-') for item in substModel_GTR_neutral_DNA_pre]
        #INDEL models
        #ZIPF
        if line.split(':')[0] == 'indel_model_ZIPF':
           range_dict['indel_model_ZIPF']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'indel_gainRate_ZIPF':
            indel_gainRate_ZIPF_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainRate_ZIPF']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainRate_ZIPF_pre]
        if line.split(':')[0] == 'indel_gainParameters_ZIPF':
            indel_gainParameters_ZIPF_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainParameters_ZIPF']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainParameters_ZIPF_pre]
        if line.split(':')[0] == 'indel_gainMaxLen_ZIPF':
            indel_gainMaxLen_ZIPF_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainMaxLen_ZIPF']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainMaxLen_ZIPF_pre]
        if line.split(':')[0] == 'indel_lossRate_ZIPF':
            indel_lossRate_ZIPF_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_lossRate_ZIPF']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_lossRate_ZIPF_pre]
        #QG
        if line.split(':')[0] == 'indel_model_QG':
           range_dict['indel_model_QG']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'indel_gainRate_QG':
            indel_gainRate_QG_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainRate_QG']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainRate_QG_pre]
        if line.split(':')[0] == 'indel_gainParameters_QG':
            indel_gainParameters_QG_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainParameters_QG']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainParameters_QG_pre]
        if line.split(':')[0] == 'indel_gainMaxLen_QG':
            indel_gainMaxLen_QG_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainMaxLen_QG']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainMaxLen_QG_pre]
        if line.split(':')[0] == 'indel_lossRate_QG':
            indel_lossRate_QG_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_lossRate_QG']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_lossRate_QG_pre]
        # GEOM
        if line.split(':')[0] == 'indel_model_GEOM':
           range_dict['indel_model_GEOM']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'indel_gainRate_GEOM':
            indel_gainRate_GEOM_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainRate_GEOM']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainRate_GEOM_pre]
        if line.split(':')[0] == 'indel_gainParameters_GEOM':
            indel_gainParameters_GEOM_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainParameters_GEOM']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainParameters_GEOM_pre]
        if line.split(':')[0] == 'indel_gainMaxLen_GEOM':
            indel_gainMaxLen_GEOM_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainMaxLen_GEOM']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainMaxLen_GEOM_pre]
        if line.split(':')[0] == 'indel_lossRate_GEOM':
            indel_lossRate_GEOM_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_lossRate_GEOM']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_lossRate_GEOM_pre]
            
        # Custom
        if line.split(':')[0] == 'indel_model_Custom':
           range_dict['indel_model_Custom']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'indel_gainRate_Custom':
            range_dict['indel_gainRate_Custom']=[[float(item)] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'indel_gainParameters_Custom':
            indel_gainParameters_Custom_pre= line.split(':')[1].strip()
            range_dict['indel_gainParameters_Custom']=[['['+indel_gainParameters_Custom_pre+']']]
        if line.split(':')[0] == 'indel_gainMaxLen_Custom':
            indel_gainMaxLen_Custom_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainMaxLen_Custom']=[[int(item)] for item in indel_gainMaxLen_Custom_pre]

            
        #Site rate variability model
        #Gamma
        if line.split(':')[0] == 'rateVar_model_Gamma':
           range_dict['rateVar_model_Gamma']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'rateVar_model_Gamma_alpha':
            rateVar_model_Gamma_alpha_pre= line.split(':')[1].strip().split(',')
            range_dict['rateVar_model_Gamma_alpha']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in rateVar_model_Gamma_alpha_pre]
        if line.split(':')[0] == 'rateVar_model_Gamma_areas':
            rateVar_model_Gamma_areas_pre= line.split(':')[1].strip().split(',')
            range_dict['rateVar_model_Gamma_areas']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in rateVar_model_Gamma_areas_pre]
        if line.split(':')[0] == 'rateVar_model_Gamma_motifFreq':
            rateVar_model_Gamma_motifFreq_pre= line.split(':')[1].strip().split(',')
            range_dict['rateVar_model_Gamma_motifFreq']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in rateVar_model_Gamma_motifFreq_pre]
        
       #Poisson
        if line.split(':')[0] == 'rateVar_model_poisson':
           range_dict['rateVar_model_poisson']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'rateVar_model_poisson_areas':
            rateVar_model_poisson_areas_pre= line.split(':')[1].strip().split(',')
            range_dict['rateVar_model_poisson_areas']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in rateVar_model_poisson_areas_pre]
        if line.split(':')[0] == 'rateVar_model_poisson_motifFreq':
            rateVar_model_poisson_motifFreq_pre= line.split(':')[1].strip().split(',')
            range_dict['rateVar_model_poisson_motifFreq']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in rateVar_model_poisson_motifFreq_pre]
            
        #Sequence rate variability model
        if line.split(':')[0] == 'amongGeneDistr':
            range_dict['amongGeneDistr']= line.split(':')[1].strip().split('-')
        if line.split(':')[0] == 'aGAlpha':
            range_dict['aGAlpha']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        
        # Gene duplication
        if line.split(':')[0] == 'geneDuplRate':
            range_dict['geneDuplRate']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'transDupl':
            range_dict['transDupl']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'numberDupl':
            range_dict['numberDupl']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'fissionDupl':
            range_dict['fissionDupl']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'fusionDupl':
            range_dict['fusionDupl']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'P_pseudogene':
            range_dict['P_pseudogene']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'P_neofunc':
            range_dict['P_neofunc']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'P_subfunc':
            range_dict['P_subfunc']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]                      
        # Gene deletion
        if line.split(':')[0] == 'geneLossRate':
            range_dict['geneLossRate']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'numberLoss':
            range_dict['numberLoss']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        #Lateral Gene Transfer
        if line.split(':')[0] == 'lgtRate':
            range_dict['lgtRate']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'orthRep':
            range_dict['orthRep']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'lgtGRate':
            range_dict['lgtGRate']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'lgtGSize':
            range_dict['lgtGSize']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        #Genome Rearrangement
        if line.split(':')[0] == 'invers':
            range_dict['invers']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'invSize':
            range_dict['invSize']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'transloc':
            range_dict['transloc']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'transSize':
            range_dict['transSize']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'invtrans':
            range_dict['invtrans']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        #Gene fission and fusion without previous duplication
        if line.split(':')[0] == 'fissionRate':
            range_dict['fissionRate']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'fusionRate':
            range_dict['fusionRate']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        if line.split(':')[0] == 'numberFusion':
            range_dict['numberFusion']=numpy.linspace(float(line.split(':')[1].strip().split('-')[0]),float(line.split(':')[1].strip().split('-')[1]),float(line.split(':')[1].strip().split('-')[2])).tolist() if \
                              len(line.split(':')[1].strip().split('-'))==3 else [float(line.split(':')[1].strip())]
        #DAWG parameters
        #GTR
        if line.split(':')[0] == 'substModel_GTR_nc':
           range_dict['substModel_GTR_nc']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'substModel_GTR_alpha_nc':
            substModel_TN93_alpha1_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_GTR_alpha_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_TN93_alpha1_pre]
        if line.split(':')[0] == 'substModel_GTR_beta_nc':
            substModel_GTR_beta_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_GTR_beta_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_GTR_beta_pre]
        if line.split(':')[0] == 'substModel_GTR_gamma_nc':
            substModel_GTR_gamma_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_GTR_gamma_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_GTR_gamma_pre]
        if line.split(':')[0] == 'substModel_GTR_delta_nc':
            substModel_GTR_delta_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_GTR_delta_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_GTR_delta_pre]
        if line.split(':')[0] == 'substModel_GTR_epsilon_nc':
            substModel_GTR_epsilon_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_GTR_epsilon_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_GTR_epsilon_pre]
        if line.split(':')[0] == 'substModel_GTR_zeta_nc':
            substModel_GTR_zeta_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_GTR_zeta_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_GTR_zeta_pre]
        if line.split(':')[0] == 'substModel_GTR_frequencies_nc':
            substModel_GTR_frequencies_pre= [item.replace('{','').replace("}",'') for item in line.split(':')[1].strip().split(',')]
            range_dict['substModel_GTR_frequencies_nc']=[[item] for item in substModel_GTR_frequencies_pre]
        
        #TN93
        if line.split(':')[0] == 'substModel_TN93_nc':
           range_dict['substModel_TN93_nc']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'substModel_TN93_alpha1_nc':
            substModel_TN93_alpha1_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_TN93_alpha1_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_TN93_alpha1_pre]
        if line.split(':')[0] == 'substModel_TN93_alpha2_nc':
            substModel_TN93_alpha2_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_TN93_alpha2_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_TN93_alpha2_pre]
        if line.split(':')[0] == 'substModel_TN93_beta_nc':
            substModel_TN93_beta_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_TN93_beta_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_TN93_beta_pre]
        if line.split(':')[0] == 'substModel_TN93_frequencies_nc':
            substModel_TN93_frequencies_pre= [item.replace('{','').replace("}",'') for item in line.split(':')[1].strip().split(',')]
            range_dict['substModel_TN93_frequencies_nc']=[[item] for item in substModel_TN93_frequencies_pre]
            
        #HKY 
        if line.split(':')[0] == 'substModel_HKY_nc':
           range_dict['substModel_HKY_nc']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'substModel_HKY_alpha_nc':
            substModel_HKY_alpha_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_HKY_alpha_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_HKY_alpha_pre]
        if line.split(':')[0] == 'substModel_HKY_beta_nc':
            substModel_HKY_beta_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_HKY_beta_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_HKY_beta_pre]
        if line.split(':')[0] == 'substModel_HKY_frequencies_nc':
            substModel_HKY_frequencies_pre= [item.replace('{','').replace("}",'') for item in line.split(':')[1].strip().split(',')]
            range_dict['substModel_HKY_frequencies_nc']=[[item] for item in substModel_HKY_frequencies_pre]
         #F84
        if line.split(':')[0] == 'substModel_F84_nc':
           range_dict['substModel_F84_nc']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'substModel_F84_kappa_nc':
            substModel_F84_kappa_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_F84_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_F84_kappa_pre]
        if line.split(':')[0] == 'substModel_F84_frequencies_nc':
            substModel_F84_frequencies_pre= [item.replace('{','').replace("}",'') for item in line.split(':')[1].strip().split(',')]
            range_dict['substModel_F84_frequencies_nc']=[[item] for item in substModel_F84_frequencies_pre]
         #K2P
        if line.split(':')[0] == 'substModel_K2P_nc':
           range_dict['substModel_K2P_nc']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'substModel_K2P_transition_nc':
            substModel_F84_kappa_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_K2P_transition_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_F84_kappa_pre]
        if line.split(':')[0] == 'substModel_K2P_transversion_nc':
            substModel_F84_kappa_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_K2P_transversion_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_F84_kappa_pre]
        if line.split(':')[0] == 'substModel_K2P_frequencies_nc':
            substModel_F84_frequencies_pre= [item.replace('{','').replace("}",'') for item in line.split(':')[1].strip().split(',')]
            range_dict['substModel_K2P_frequencies_nc']=[[item] for item in substModel_F84_frequencies_pre]
            
         #K3P
        if line.split(':')[0] == 'substModel_K3P_nc':
           range_dict['substModel_K3P_nc']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'substModel_K3P_alpha_nc':
            substModel_F84_kappa_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_K3P_alpha_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_F84_kappa_pre]
        if line.split(':')[0] == 'substModel_K3P_beta_nc':
            substModel_F84_kappa_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_K3P_beta_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_F84_kappa_pre]
        if line.split(':')[0] == 'substModel_K3P_gamma_nc':
            substModel_F84_kappa_pre= line.split(':')[1].strip().split(',')
            range_dict['substModel_K3P_gamma_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in substModel_F84_kappa_pre]
        if line.split(':')[0] == 'substModel_K3P_frequencies_nc':
            substModel_F84_frequencies_pre= [item.replace('{','').replace("}",'') for item in line.split(':')[1].strip().split(',')]
            range_dict['substModel_K3P_frequencies_nc']=[[item] for item in substModel_F84_frequencies_pre]
        #Noncoding DNA indel model   
        # GEOM
        if line.split(':')[0] == 'indel_model_GEOM_nc':
           range_dict['indel_model_GEOM_nc']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'indel_gainRate_GEOM_nc':
            indel_gainRate_GEOM_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainRate_GEOM_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainRate_GEOM_pre]
        if line.split(':')[0] == 'indel_gainParameters_GEOM_nc':
            indel_gainParameters_GEOM_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainParameters_GEOM_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainParameters_GEOM_pre]
        if line.split(':')[0] == 'indel_gainMaxLen_GEOM_nc':
            indel_gainMaxLen_GEOM_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainMaxLen_GEOM_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainMaxLen_GEOM_pre]
        # Power_law
        if line.split(':')[0] == 'indel_model_PL_nc':
           range_dict['indel_model_PL_nc']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'indel_gainRate_PL_nc':
            indel_gainRate_PL_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainRate_PL_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainRate_PL_pre]
        if line.split(':')[0] == 'indel_gainParameters_PL_nc':
            indel_gainParameters_PL_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainParameters_PL_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainParameters_PL_pre]
        if line.split(':')[0] == 'indel_gainMaxLen_PL_nc':
            indel_gainMaxLen_PL_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainMaxLen_PL_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainMaxLen_PL_pre]
        # US
        if line.split(':')[0] == 'indel_model_US_nc':
           range_dict['indel_model_US_nc']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'indel_gainRate_US_nc':
            indel_gainRate_US_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainRate_US_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in indel_gainRate_US_pre]
        if line.split(':')[0] == 'indel_gainParameters_US_nc':
            indel_gainParameters_US_pre= line.split(':')[1].strip().split(',')
            range_dict['indel_gainParameters_US_nc']=[[','.join(item.split(';'))] for item in indel_gainParameters_US_pre]

            
        #Site rate variability model
        #Gamma
        if line.split(':')[0] == 'rateVar_model_Gamma_nc':
           range_dict['rateVar_model_Gamma_nc']= [[item] for item in line.split(':')[1].strip().split(',')]
        if line.split(':')[0] == 'rateVar_model_Gamma_alpha_nc':
            rateVar_model_Gamma_alpha_pre= line.split(':')[1].strip().split(',')
            range_dict['rateVar_model_Gamma_alpha_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),int(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in rateVar_model_Gamma_alpha_pre]
        if line.split(':')[0] == 'rateVar_model_Gamma_invariant_site_ratio_nc':
            rateVar_model_Gamma_motifFreq_pre= line.split(':')[1].strip().split(',')
            range_dict['rateVar_model_Gamma_invariant_site_ratio_nc']=[numpy.linspace(float(item.split('-')[0]),float(item.split('-')[1]),float(item.split('-')[2])).tolist() if \
                              len(item.split('-'))==3 else [float(item)] for item in rateVar_model_Gamma_motifFreq_pre]
            

            
            
        
        
        