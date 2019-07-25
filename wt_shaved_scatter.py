import matplotlib.pyplot as plt
import pandas as pd
import os
import pickle
import physplot as pp
import physcalc as pc
import physfiles as pf
import numpy as np
import math


# import experiment class
dir_in = "C:\\Users\\haley\\Dropbox\\Projects\\FMRP\\Spreadsheets"
os.chdir(dir_in)
class_exp = 'wt_mEPSC.p'
exp = pf.get_class(class_exp)
exp.group_labels = ['Ipsilateral', 'Contralateral']
exp.dir_in = dir_in

# Default Settings
exp.dir_plots = dir_in + '\\Plots'
exp.x_label_text = exp.group_labels
exp.x_label_size = 12
exp.y_label_size = 14
exp.x_ticklabel_angle = 45
exp.x_title_size = 16
exp.y_title_size = 16
exp.x_title_text = ''
exp.scatter_marker_size = 3
exp.marker_colors = ['black', 'dimgray']
exp.markers = ['s','o']
exp.avg_marker_size = 25
exp.avg_marker = '_'
exp.avg_marker_color = 'black'
exp.avg_err_size = 1
exp.avg_err_color = 'black'
exp.err_cap_size = 5
exp.err_line_width = 1
exp.posthoc_stars_size = 12
exp.plot_dpi = 80     
exp.save_dpi = 200

# Import raw data
list_scatter = pf.get_csv (exp.dir_in)

for df_scatter in list_scatter:   
    # Import desc stats and t-tests
    if df_scatter.iloc[0]['exp'] == 'vsteps':
        df_group1 = pd.read_csv(exp.dir_in + exp.dir_vsteps + '\\' + 'wt_vsteps_ipsilateral_desc.csv')
        df_group2 = pd.read_csv(exp.dir_in + exp.dir_vsteps + '\\' + 'wt_vsteps_contralateral_desc.csv')
        
        pp.iv_plot (exp, df_group1, df_group2)

    elif df_scatter.iloc[0]['exp'] == 'mEPSC':
        df_group1 = pd.read_csv(exp.dir_in + exp.dir_mEPSCs + '\\' + 'wt_mEPSC_ipsilateral_desc.csv')
        df_group2 = pd.read_csv(exp.dir_in + exp.dir_mEPSCs + '\\' + 'wt_mEPSC_contralateral_desc.csv')
        df_ttest = pd.read_csv(exp.dir_in + exp.dir_mEPSCs + '\\' + 'wt_mEPSC_t-test.csv')
        df_group1 = df_group1.set_index(['stat'])
        df_group2 = df_group2.set_index(['stat'])
        df_ttest = df_ttest.set_index(['measure'])

        # Plot each measure for this experiment type
        for measure in exp.measures_minis:
            if measure == 'amp':
                exp.y_label_text = 'mEPSC Amplitude (pA)'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
            elif measure == 'q':
                exp.y_label_text = 'mEPSC Charge (pA*ms)'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
            elif measure == 'half_width':
                exp.y_label_text = 'mEPSC Half-Width (ms)'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
            elif measure == 'decay':
                exp.y_label_text = 'mEPSC Decay (s' + r'$^-$' + r'$^1$' + ')'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
            elif measure == 'freq':
                exp.y_label_text = 'mEPSC Frequency (Hz)'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)

    
    
   
    elif df_scatter.iloc[0]['exp'] == 'sEPSC':
        df_group1 = pd.read_csv(exp.dir_in + exp.dir_sEPSCs + '\\' + 'wt_sEPSC_ipsilateral_desc.csv')
        df_group2 = pd.read_csv(exp.dir_in + exp.dir_sEPSCs + '\\' + 'wt_sEPSC_contralateral_desc.csv')
        df_ttest = pd.read_csv(exp.dir_in + exp.dir_sEPSCs + '\\' + 'wt_sEPSC_t-test.csv')
        df_group1 = df_group1.set_index(['stat'])
        df_group2 = df_group2.set_index(['stat'])
        df_ttest = df_ttest.set_index(['measure'])
        
        for measure in exp.measures_minis:
            if measure == 'amp':
                exp.y_label_text = 'sEPSC Amplitude (pA)'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
            elif measure == 'q':
                exp.y_label_text = 'sEPSC Charge (pA*ms)'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
            elif measure == 'half_width':
                exp.y_label_text = 'sEPSC Half-Width (ms)'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
            elif measure == 'decay':
                exp.y_label_text = 'sEPSC Decay (s' + r'$^-$' + r'$^1$' + ')'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
            elif measure == 'freq':
                exp.y_label_text = 'sEPSC Frequency (Hz)'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
    
                

    elif df_scatter.iloc[0]['exp'] == 'sEPSP':
        df_group1 = pd.read_csv(exp.dir_in + exp.dir_sEPSPs + '\\' + 'wt_sEPSP_ipsilateral_desc.csv')
        df_group2 = pd.read_csv(exp.dir_in + exp.dir_sEPSPs + '\\' + 'wt_sEPSP_contralateral_desc.csv')
        df_ttest = pd.read_csv(exp.dir_in + exp.dir_sEPSPs + '\\' + 'wt_sEPSP_t-test.csv')
        df_group1 = df_group1.set_index(['stat'])
        df_group2 = df_group2.set_index(['stat'])
        df_ttest = df_ttest.set_index(['measure'])
        
        for measure in exp.measures_minis:
            # If amp, then these things
            if measure == 'amp':
                exp.y_label_text = 'sEPSP Amplitude (pA)'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
            elif measure == 'q':
                exp.y_label_text = 'sEPSP Charge (pA*ms)'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
            elif measure == 'half_width':
                exp.y_label_text = 'sEPSP Half-Width (ms)'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
            elif measure == 'decay':
                exp.y_label_text = 'sEPSP Decay (s' + r'$^-$' + r'$^1$' + ')'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
            elif measure == 'freq':
                exp.y_label_text = 'sEPSP Frequency (Hz)'
                pp.mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure)
    