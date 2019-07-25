import os
import pandas as pd
import physcalc as pc
import physfiles as pf
import physstats as ps
import physplot as pp
import pickle
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np


# Class data for project
dir_in = "C:\\Users\\haley\\Dropbox\\Projects\\FMRP\\Spreadsheets"
#dir_in = "D:\\Dropbox\\Projects\\FMRP\\Spreadsheets"
dir_group = "\\Descriptive_Stats\\By_Group"
dir_measure = "\\Descriptive_Stats\\By_Measurement"
dir_tukey = '\\Posthoc_Tukey\\' 
dir_anova = '\\One-Way_ANOVA'
dir_ttest = '\\T-Test'
dir_isteps = '\\Isteps'
dir_vsteps = '\\Vsteps'
dir_sEPSCs = "\\sEPSC"
dir_mEPSCs = '\\mEPSC'
dir_sEPSPs = '\\sEPSP'
dir_mEPSPs = '\\mEPSP'
    
# Experiment specfic data
exp_type = 'mEPSC' 
exp_name = 'wt'
exp_groups = list(['ipsilateral', 'contralateral'])

between_factor = 'lateral'
measures_common = ['exp','file','genotype','slice','shaved_side','hemisphere','lateral']
measures_minis = ['amp','q', 'half_width', 'decay', 'freq']
measures_vsteps = ['vm-50','vm-40','vm-30','vm-20','vm-10','vm0','vm+10','vm+20','vm+30','vm+40','vm+50','vm+60']
measures_isteps = ['holding-150','holding-125','holding-100','holding-75','holding-50','holding-25','holding','holding+25','holding+50','holding+75']
confidence_interval = 0.95

# Initiate an experiment class and a mini class
print ('Initializing ........')
exp = pc.mini(exp_type, exp_name, exp_groups, between_factor, confidence_interval, dir_in, dir_group, dir_measure,
                dir_tukey, dir_anova, dir_ttest, dir_isteps, dir_vsteps, dir_sEPSCs, dir_mEPSCs, 
                dir_sEPSPs, dir_mEPSPs, measures_minis, measures_isteps, measures_vsteps, measures_common)       


print ('Importing data ........')
# Import each CSV in the directory into dataframes (*****alpha-numeric order!******)
list_dfs = pf.get_csv (dir_in)

# Do basic calculations for cc, gj, rin, ect
print('Calculating Descriptive Stats.......')
list_desc = list()
for df in list_dfs:
    list_groupdfs = pf.get_groups (exp, df)
    if df.iloc[0]['exp'] == 'mEPSC':
        exp.dir_out = exp.dir_in + exp.dir_mEPSCs 
        
        # Get t-tests for minis
        df_t_mEPSC = ps.ttest_minis (exp, list_groupdfs) 
        file_out = exp.name + '_mEPSC_t-test.csv'
        pf.save_csv (df_t_mEPSC, file_out, exp.dir_out, True) 
        
        # Get descriptive stats
        for df_group in list_groupdfs:
            df_desc = ps.get_desc_mini (exp, df_group)  
            list_desc.append(df_desc)   
        print('pass: ', df.iloc[0]['exp'])  

    elif df.iloc[0]['exp'] == 'mEPSP':
        exp.dir_out = exp.dir_in + exp.dir_mEPSPs 
        
        # Get t-tests for minis
        df_t_mEPSP = ps.ttest_minis (exp, list_groupdfs)
        file_out = exp.name + '_mEPSP_t-test.csv'
        pf.save_csv (df_t_mEPSP, file_out, exp.dir_out, True)
        
        # Get descriptive stats
        for df_group in list_groupdfs:
            df_desc = ps.get_desc_mini (exp, df_group)
            list_desc.append(df_desc)
        print('pass: ', df.iloc[0]['exp'])

    elif df.iloc[0]['exp'] == 'sEPSC':
        exp.dir_out = exp.dir_in + exp.dir_sEPSCs 
        
        # Get t-tests for minis 
        df_t_sEPSC = ps.ttest_minis (exp, list_groupdfs)
        file_out = exp.name + '_sEPSC_t-test.csv'
        pf.save_csv (df_t_mEPSC, file_out, exp.dir_out, True)
        
        # Get descriptive stats
        for df_group in list_groupdfs:
            df_desc = ps.get_desc_mini (exp, df_group)
            list_desc.append(df_desc)
        print('pass: ', df.iloc[0]['exp'])

    elif df.iloc[0]['exp'] == 'sEPSP':
        exp.dir_out = exp.dir_in + exp.dir_sEPSPs
        
        # Get t-tests for minis
        df_t_sEPSP = ps.ttest_minis (exp, list_groupdfs)
        file_out = exp.name + '_sEPSP_t-test.csv'
        pf.save_csv (df_t_sEPSP, file_out, exp.dir_out, True)
        
        # Get descriptive stats
        for df_group in list_groupdfs:
            df_desc = ps.get_desc_mini (exp, df_group)
            list_desc.append(df_desc)
        print('pass: ', df.iloc[0]['exp'])

    elif df.iloc[0]['exp'] == 'isteps':
        exp.dir_out = exp.dir_in + '//isteps'

        # Make a spreadsheet with indexed data for 2-way RM ANOVA in OriginLab Pro 2016
        df_isteps = ps.get_indexed_seps(exp, df, exp.measures_isteps) # Get indexed data for anova
        file_out = exp.name + '_isteps_indexed.csv'
        pf.save_csv (df_isteps, file_out, exp.dir_out, True) 

        # Get t-test for Vm
        df_t_isteps = ps.ttest_isteps (exp, list_groupdfs) 
        file_out = exp.name + '_isteps_t-test.csv'
        pf.save_csv (df_t_isteps, file_out, exp.dir_out, True) 
        
        # Get descriptive stats on steps
        for df_group in list_groupdfs:
            df_desc = ps.get_desc_isteps (exp, df_group)
            list_desc.append(df_desc)
            print('pass: ', df.iloc[0]['exp'])

    elif df.iloc[0]['exp'] == 'vsteps':
        exp.dir_out = exp.dir_in + '//vsteps'
        
        # Make a spreadsheet with indexed data for 2-way RM ANOVA in OriginLab Pro 2016
        df_vsteps = ps.get_indexed_seps(exp, df, exp.measures_vsteps) 
        file_out = exp.name + '_vsteps_indexed.csv'
        pf.save_csv (df_vsteps, file_out, exp.dir_out, True) 

        # Get t-test for holding current
        df_t_vsteps = ps.ttest_vsteps (exp, list_groupdfs)
        file_out = exp.name + '_vsteps_t-test.csv'  
        pf.save_csv (df_t_vsteps, file_out, exp.dir_out, True) 

        ## Calculate Rinput and save back to the original file
        #df['r_input'] = np.nan
        #df['r_input'] = -20 / df['vm-20'] * 1000
        #file_out = '190702_vsteps.csv'
        #pf.save_csv (df, file_out, exp.dir_in, True) 

        # Get t-test for Input Resistance
        #df_t_rin = ps.ttest_rinput (exp, list_groupdfs)
        #file_out = exp.name + '_rinput_t-test.csv'  
        #pf.save_csv (df_t_rin, file_out, exp.dir_out, True)

        # Get descriptive stats for steps
        for df_group in list_groupdfs:
            df_desc = ps.get_desc_vsteps (exp, df_group)
            list_desc.append(df_desc)
            print('pass: ', df.iloc[0]['exp'])





print ('Saving experiment parameters to file.')
pf.save_class (exp)


print ('Analysis complete.')