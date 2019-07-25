import matplotlib.pyplot as plt
import pandas as pd
import os
import pickle
import physplot as pp
import physcalc as pc
import physfiles as pf
import phystrace as pt
import numpy as np
import math
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

# import experiment class
dir_in = "C:\\Users\\haley\\Dropbox\\Projects\\FMRP\\Spreadsheets"
os.chdir(dir_in)
class_exp = 'wt_mEPSC.p'
exp = pf.get_class(class_exp)
exp.group_labels = ['Ipsilateral', 'Contralateral']
exp.dir_in = dir_in
exp.dir_traces = '\\Traces'
# Default Settings
exp.dir_plots = dir_in + '\\Plots'
exp.x_label_text = exp.group_labels
exp.x_label_size = 12
exp.y_label_size = 14
exp.y_tick_length = 10
exp.y_tick_units = 'pA'
exp.x_tick_length_20 = 2000
exp.x_tick_units = 's'
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

# Get spreadsheets for mEPSCs, mEPSPs, sEPSCs
df_mC = pd.read_csv(exp.dir_in + exp.dir_traces + '\\' + 'wt_mEPSCs_traces.csv')
df_sC = pd.read_csv(exp.dir_in + exp.dir_traces + '\\' +'wt_sEPSCs_traces.csv')
df_sP = pd.read_csv(exp.dir_in + exp.dir_traces + '\\' +'wt_sEPSPs_traces.csv')

# Normalize all traces to 0
df_mC = pt.sub_baseline (df_mC)
df_mP = pt.sub_baseline (df_sC)
df_sP = pt.sub_baseline (df_sP) 

pt.mEPSC_traces(exp, df_mC)
pt.sEPSC_traces(exp,df_sC)
pt.sEPSP_traces(exp, df_sP)
