import pandas as pd 
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import numpy as np 
import math
import physcalc as pc
import physfiles as pf

class trace:
    def __init__(self):
        self.plot_dpi = 80                  
        self.save_dpi = 300 


def set_padding(trace):
    trace.x_min = trace.x_min - trace.x_min * trace.x_pad
    trace.x_max = trace.x_max + trace.x_max + trace.x_pad
    trace.y_stim_min = trace.y_stim_min - trace.y_stim_min * trace.y_pad
    trace.y_stim_max = trace.y_stim_max + trace.y_stim_max * trace.y_pad

def set_range (trace, df):
    # Get original vm
    trace.x_range = trace.x_max - trace.x_min
    trace.vm_stim = df.iloc[0:100]['stim'].mean()
    trace.vm_rec = df.iloc[0:100]['rec'].mean() 
    trace.vm_stim_nmda = df.iloc[0:100]['stim_nmda'].mean()
    trace.vm_rec_nmda = df.iloc[0:100]['rec_nmda'].mean() 

    # Set the largest scale that will accommodate all traces
    if df['stim'].min() <= df['stim_nmda'].min():
        y_stim_min = df['stim']
    else:
        y_stim_min = df['stim_nmda']
    if df['stim'].max() >= df['stim_nmda'].max():
        y_stim_max = df['stim']
    else:
        y_stim_max = df['stim_nmda']

    # Set the scales for stim and record
    trace.y_stim_min = pc.round_down(min(y_stim_min), 5)
    trace.y_stim_max = pc.round_up(max(y_stim_max), 5)
    trace.y_stim_range = abs(trace.y_stim_max - trace.y_stim_min)

    trace.y_rec_min = pc.round_down(df['rec'].min(), 5)
    trace.y_rec_max = pc.round_up(df['rec'].max(), 5)
    trace.y_rec_range = abs(trace.y_rec_max - trace.y_rec_min)

    trace.y_stim_nmda_min = trace.y_stim_min
    trace.y_stim_nmda_max = trace.y_stim_max
    trace.y_stim_nmda_range = trace.y_stim_nmda_max - trace.y_stim_nmda_min

    # Adjust Vm of 1st rec to match 1st stim
    diff_rec = trace.vm_rec - trace.vm_stim
    sf = trace.y_stim_range/trace.y_rec_range
    df['adj_rec'] = df['rec'] - diff_rec/sf
    trace.y_rec_min = pc.round_down(df['adj_rec'].min(), 5)
    trace.y_rec_max = pc.round_up(df['adj_rec'].max(), 5)
    trace.y_rec_range = abs(trace.y_rec_max - trace.y_rec_min)
    trace.vm_rec = df.iloc[0]['adj_rec']
    
    # Adjust Vm of 2nd stim to be level with 1st stim
    diff_stim = trace.vm_stim_nmda - trace.vm_stim
    df['adj_stim_nmda'] = df['stim_nmda'] - diff_stim
    trace.y_stim_nmda_range = max(df['adj_stim_nmda']) - min(df['adj_stim_nmda'])
    trace.vm_stim_nmda = df.iloc[0]['adj_stim_nmda']

    # Adjust Vm and scale of 2nd rec to line up with Vm of stim
    diff_rec = trace.vm_rec_nmda - trace.vm_rec
    df['adj_rec_nmda'] = df['rec_nmda'] - (-8.916928517799995)
    trace.y_rec_nmda_min = trace.y_rec_min
    trace.y_rec_nmda_max = trace.y_rec_max
    trace.y_rec_nmda_range = abs(trace.y_rec_max - trace.y_rec_min)
    trace.vm_rec_nmda = df.iloc[0]['adj_rec_nmda']
    print('adj: ', trace.vm_rec_nmda)
    
    # Set horizontal line
    trace.hline_stim_label = str(df.iloc[0:100]['stim'].mean()) + ' mV'
    trace.hline_stim_nmda_label = str(df.iloc[0:100]['rec'].mean()) + 'mV'
    trace.hline_rec_label = str(df.iloc[0:100]['stim_nmda'].mean()) + ' mV'
    trace.hline_rec_nmda_label = str(df.iloc[0:100]['stim_nmda'].mean()) + 'mV'
    return df

def set_axes (trace, ax1, ax2, ax3, ax4):
    ax1.set_ylim(bottom = trace.y_stim_min, top = trace.y_stim_max)
    ax1.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax1.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax2.set_ylim(bottom =trace.y_rec_min, top = trace.y_rec_max)
    ax2.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax2.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax3.set_ylim(bottom = trace.y_stim_nmda_min, top = trace.y_stim_nmda_max)
    ax3.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax3.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax4.set_ylim(bottom = trace.y_rec_nmda_min, top = trace.y_rec_nmda_max)
    ax4.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax4.yaxis.set_major_locator(plt.MaxNLocator(5))
   
def hide_bounds (list_ax):
    for ax in list_ax:
        ax.spines['top'].set_visible (False)
        ax.spines['right'].set_visible (False)
        ax.spines['left'].set_visible (False)
        ax.spines['bottom'].set_visible (False)

def hide_labels (list_ax):
    for ax in list_ax:
        # Remove ticks and labels
        ax.yaxis.set_ticks([])
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_ticks([])
        ax.xaxis.set_ticklabels([])

def get_legend (list_ax):
    i = 0
    for ax in list_ax:
        if i == 0 or i == 2:
            ax.legend(loc = 'lower left', bbox_to_anchor=(0, -0.1), frameon = False) 
        if i == 1 or i == 3:
            ax.legend(loc = 'lower left', bbox_to_anchor=(0, -0.2), frameon = False)
        i += 1
    
def get_scalebar (trace, list_ax):
    trace.x_tick_length = int(trace.x_range / 2)
    trace.y_stim_tick_length = trace.y_stim_range / 5
    #trace.y_rec_tick_length = trace.y_rec_range / 0.005     

    x_label = str(trace.x_tick_length) + ' ' + trace.x_tick_units
    # needs to be individualized for each ax ------------------------------------------------------------------------------------------------
    y_label1 = ' ' + str(y_tick_length1) + ' ' + y_tick_units1
    #y_label2 = ' ' + str(y_tick_length2) + ' ' + y_tick_units2

    x_scalebar = AnchoredSizeBar(list_ax[0].transData, trace.x_tick_length, '', 'center right', pad=0, color='black', frameon=False, size_vertical= .2)
    #y_scalebar1 = AnchoredSizeBar(ax1.transData, 1.4, '', 'lower right', pad=0, color='black', frameon=False, size_vertical= y_tick_length1)
    list_ax[0].add_artist(x_scalebar)
    #list_ax[0].add_artist(y_scalebar1)(self, parameter_list)

    # Add scale and units to error bars
    list_ax[0].text(trace.x_max + 75, trace.y_stim_max - trace.y_stim_range/2 - 1, x_label, size = 6)
    #ax1.text(x_max, y_min1 + 8, y_label1, size = 12, color = line_color1)
    #ax2.text(x_max, y_min2 +.5, y_label2, size = 12, color = line_color2)

def sub_baseline (df):
    
    for i in df.columns.values:
        if i != 'time':
            diff = 0 - df.iloc[0][i]
            df[i] += diff   
    return df

def mEPSC_traces(exp, df):
    x_data = df['time']

    exp.x_min = 0
    exp.x_max_2 = 900
    exp.x_max_20 = 20000
    exp.y_min = -30
    exp.y_max = 10
    exp.y_range = exp.y_max-exp.y_min

    # Plot the data
    fig, ax = plt.subplots(2,2, figsize=(4, 2), dpi = exp.save_dpi)
    plt.subplots_adjust(hspace = 0.5)

    ax[0,0].plot(x_data, df.ipsilateral_20s, ls = 'solid', lw = .5, color = 'black', label = 'Ipsilateral') 
    ax[0,0].axis ([exp.x_min, exp.x_max_20, exp.y_min, exp.y_max]) 
    ax[0,0].spines['top'].set_visible (False)
    ax[0,0].spines['right'].set_visible (False)
    ax[0,0].spines['left'].set_visible (False)
    ax[0,0].spines['bottom'].set_visible (False)


    # Set tick lengths
    ax[0,0].xaxis.set_ticks(np.arange(exp.x_min, exp.x_max_20, exp.x_tick_length_20))
    ax[0,0].yaxis.set_ticks(np.arange(exp.y_min, exp.y_max, exp.y_tick_length))

    # Create scale bar
    x_label = str(exp.x_tick_length_20/1000) + ' ' + exp.x_tick_units
    y_label = ' ' + str(exp.y_tick_length) + ' ' + exp.y_tick_units
    x_scalebar = AnchoredSizeBar(ax[0,0].transData, exp.x_tick_length_20, '', 'lower right', pad=-2, color='black', frameon=False, size_vertical= 0.2)
    y_scalebar = AnchoredSizeBar(ax[0,0].transData, 0, '', 'lower right', pad=-2, color='black', frameon=False, size_vertical= exp.y_tick_length)
    ax[0,0].add_artist(x_scalebar)
    ax[0,0].add_artist(y_scalebar)

    # SB = 2s of 20s trace
    # SB = 0.1333 s of a 20s trace

    ax[0,0].text(exp.x_min, exp.y_max+2, 'Ipsilateral', size = 8)
    #ax[0,0].text(exp.x_max_20 + 200, exp.y_min + 5, y_label, size = 6, color = 'black')

    ax[0,0].yaxis.set_ticks([]) 
    ax[0,0].yaxis.set_ticklabels([])
    ax[0,0].xaxis.set_ticks([])
    ax[0,0].xaxis.set_ticklabels([])

    ax[1,0].plot(x_data, df.contralateral_20s, ls = 'solid', lw = .5, color = 'dimgray', label = 'Contralateral') 
    ax[1,0].axis ([exp.x_min, exp.x_max_20, exp.y_min, exp.y_max]) 
    ax[1,0].spines['top'].set_visible (False)
    ax[1,0].spines['right'].set_visible (False)
    ax[1,0].spines['left'].set_visible (False)
    ax[1,0].spines['bottom'].set_visible (False)
    ax[1,0].yaxis.set_ticks([])
    ax[1,0].yaxis.set_ticklabels([])
    ax[1,0].xaxis.set_ticks([])
    ax[1,0].xaxis.set_ticklabels([])
    ax[1,0].text(exp.x_min, exp.y_max+ 2, 'Contralateral', size = 8)

    ax[0,1].plot(x_data, df.ipsilateral_2s, ls = 'solid', lw = .5, color = 'black', label = 'Ipsilateral') 
    ax[0,1].axis ([exp.x_min, exp.x_max_2, exp.y_min, exp.y_max]) 
    ax[0,1].spines['top'].set_visible (False)
    ax[0,1].spines['right'].set_visible (False)
    ax[0,1].spines['left'].set_visible (False)
    ax[0,1].spines['bottom'].set_visible (False)
    ax[0,1].yaxis.set_ticks([])
    ax[0,1].yaxis.set_ticklabels([])
    ax[0,1].xaxis.set_ticks([])
    ax[0,1].xaxis.set_ticklabels([])

    ax[1,1].plot(x_data, df.contralateral_2s, ls = 'solid', lw = .5, color = 'dimgray', label = 'Contralateral') 
    ax[1,1].axis ([exp.x_min, exp.x_max_2, exp.y_min, exp.y_max]) 
    ax[1,1].spines['top'].set_visible (False)
    ax[1,1].spines['right'].set_visible (False)
    ax[1,1].spines['left'].set_visible (False)
    ax[1,1].spines['bottom'].set_visible (False)
    ax[1,1].yaxis.set_ticks([])
    ax[1,1].yaxis.set_ticklabels([])
    ax[1,1].xaxis.set_ticks([])
    ax[1,1].xaxis.set_ticklabels([])

    #plt.show()
    exp.file_fig = exp.name + '_mEPSC_traces.png' 
    pf.save_fig (exp, fig, exp.file_fig, exp.dir_in + exp.dir_traces)


def sEPSC_traces(exp, df):
    x_data = df['time']

    exp.x_min = 0
    exp.x_max_2 = 900
    exp.x_max_20 = 20000
    exp.y_min = -30
    exp.y_max = 10
    exp.y_range = exp.y_max-exp.y_min
    exp.y_tick_length = 10

    # Plot the data
    fig, ax = plt.subplots(2,2, figsize=(4, 2), dpi = exp.save_dpi)
    plt.subplots_adjust(hspace = 0.5)

    ax[0,0].plot(x_data, df.ipsilateral_20s, ls = 'solid', lw = .5, color = 'black', label = 'Ipsilateral') 
    ax[0,0].axis ([exp.x_min, exp.x_max_20, exp.y_min, exp.y_max]) 
    ax[0,0].spines['top'].set_visible (False)
    ax[0,0].spines['right'].set_visible (False)
    ax[0,0].spines['left'].set_visible (False)
    ax[0,0].spines['bottom'].set_visible (False)


    # Set tick lengths
    ax[0,0].xaxis.set_ticks(np.arange(exp.x_min, exp.x_max_20, exp.x_tick_length_20))
    ax[0,0].yaxis.set_ticks(np.arange(exp.y_min, exp.y_max, exp.y_tick_length))

    # Create scale bar
    x_label = str(exp.x_tick_length_20/1000) + ' ' + exp.x_tick_units
    y_label = ' ' + str(exp.y_tick_length) + ' ' + exp.y_tick_units
    x_scalebar = AnchoredSizeBar(ax[0,0].transData, exp.x_tick_length_20, '', 'lower right', pad=-2, color='black', frameon=False, size_vertical= 0.2)
    y_scalebar = AnchoredSizeBar(ax[0,0].transData, 0, '', 'lower right', pad=-2, color='black', frameon=False, size_vertical= exp.y_tick_length)
    ax[0,0].add_artist(x_scalebar)
    ax[0,0].add_artist(y_scalebar)

    # SB = 2s of 20s trace
    # SB = 0.1333 s of a 20s trace

    ax[0,0].text(exp.x_min, exp.y_max+2, 'Ipsilateral', size = 8)
    #ax[0,0].text(exp.x_max_20 + 200, exp.y_min + 5, y_label, size = 6, color = 'black')

    ax[0,0].yaxis.set_ticks([]) 
    ax[0,0].yaxis.set_ticklabels([])
    ax[0,0].xaxis.set_ticks([])
    ax[0,0].xaxis.set_ticklabels([])

    ax[1,0].plot(x_data, df.contralateral_20s, ls = 'solid', lw = .5, color = 'dimgray', label = 'Contralateral') 
    ax[1,0].axis ([exp.x_min, exp.x_max_20, exp.y_min, exp.y_max]) 
    ax[1,0].spines['top'].set_visible (False)
    ax[1,0].spines['right'].set_visible (False)
    ax[1,0].spines['left'].set_visible (False)
    ax[1,0].spines['bottom'].set_visible (False)
    ax[1,0].yaxis.set_ticks([])
    ax[1,0].yaxis.set_ticklabels([])
    ax[1,0].xaxis.set_ticks([])
    ax[1,0].xaxis.set_ticklabels([])
    ax[1,0].text(exp.x_min, exp.y_max+ 2, 'Contralateral', size = 8)

    ax[0,1].plot(x_data, df.ipsilateral_2s, ls = 'solid', lw = .5, color = 'black', label = 'Ipsilateral') 
    ax[0,1].axis ([exp.x_min, exp.x_max_2, exp.y_min, exp.y_max]) 
    ax[0,1].spines['top'].set_visible (False)
    ax[0,1].spines['right'].set_visible (False)
    ax[0,1].spines['left'].set_visible (False)
    ax[0,1].spines['bottom'].set_visible (False)
    ax[0,1].yaxis.set_ticks([])
    ax[0,1].yaxis.set_ticklabels([])
    ax[0,1].xaxis.set_ticks([])
    ax[0,1].xaxis.set_ticklabels([])

    ax[1,1].plot(x_data, df.contralateral_2s, ls = 'solid', lw = .5, color = 'dimgray', label = 'Contralateral') 
    ax[1,1].axis ([exp.x_min, exp.x_max_2, exp.y_min, exp.y_max]) 
    ax[1,1].spines['top'].set_visible (False)
    ax[1,1].spines['right'].set_visible (False)
    ax[1,1].spines['left'].set_visible (False)
    ax[1,1].spines['bottom'].set_visible (False)
    ax[1,1].yaxis.set_ticks([])
    ax[1,1].yaxis.set_ticklabels([])
    ax[1,1].xaxis.set_ticks([])
    ax[1,1].xaxis.set_ticklabels([])

    #plt.show()
    exp.file_fig = exp.name + '_sEPSC_traces.png' 
    pf.save_fig (exp, fig, exp.file_fig, exp.dir_in + exp.dir_traces)

def sEPSP_traces(exp, df):
    x_data = df['time']

    exp.x_min = 10
    exp.x_max_2 = 900
    exp.x_max_20 = 20000
    exp.y_min = -10
    exp.y_max = 5
    exp.y_range = exp.y_max-exp.y_min

    # Plot the data
    fig, ax = plt.subplots(2,2, figsize=(4, 2), dpi = exp.save_dpi)
    plt.subplots_adjust(hspace = 0.5)

    ax[0,0].plot(x_data, df.ipsilateral_20s, ls = 'solid', lw = .5, color = 'black', label = 'Ipsilateral') 
    ax[0,0].axis ([exp.x_min, exp.x_max_20, exp.y_min, exp.y_max]) 
    ax[0,0].spines['top'].set_visible (False)
    ax[0,0].spines['right'].set_visible (False)
    ax[0,0].spines['left'].set_visible (False)
    ax[0,0].spines['bottom'].set_visible (False)


    # Set tick lengths
    ax[0,0].xaxis.set_ticks(np.arange(exp.x_min, exp.x_max_20, exp.x_tick_length_20))
    ax[0,0].yaxis.set_ticks(np.arange(exp.y_min, exp.y_max, exp.y_tick_length))

    # Create scale bar
    x_label = str(exp.x_tick_length_20/1000) + ' ' + exp.x_tick_units
    y_label = ' ' + str(exp.y_tick_length) + ' ' + exp.y_tick_units
    x_scalebar = AnchoredSizeBar(ax[0,0].transData, exp.x_tick_length_20, '', 'lower right', pad=-2, color='black', frameon=False, size_vertical= 0.2)
    y_scalebar = AnchoredSizeBar(ax[0,0].transData, 0, '', 'lower right', pad=-2, color='black', frameon=False, size_vertical= exp.y_tick_length)
    ax[0,0].add_artist(x_scalebar)
    ax[0,0].add_artist(y_scalebar)

    # SB = 2s of 20s trace
    # SB = 0.1333 s of a 20s trace

    ax[0,0].text(exp.x_min, exp.y_max+2, 'Ipsilateral', size = 8)
    #ax[0,0].text(exp.x_max_20 + 200, exp.y_min + 5, y_label, size = 6, color = 'black')

    ax[0,0].yaxis.set_ticks([]) 
    ax[0,0].yaxis.set_ticklabels([])
    ax[0,0].xaxis.set_ticks([])
    ax[0,0].xaxis.set_ticklabels([])

    ax[1,0].plot(x_data, df.contralateral_20s, ls = 'solid', lw = .5, color = 'dimgray', label = 'Contralateral') 
    ax[1,0].axis ([exp.x_min, exp.x_max_20, exp.y_min, exp.y_max]) 
    ax[1,0].spines['top'].set_visible (False)
    ax[1,0].spines['right'].set_visible (False)
    ax[1,0].spines['left'].set_visible (False)
    ax[1,0].spines['bottom'].set_visible (False)
    ax[1,0].yaxis.set_ticks([])
    ax[1,0].yaxis.set_ticklabels([])
    ax[1,0].xaxis.set_ticks([])
    ax[1,0].xaxis.set_ticklabels([])
    ax[1,0].text(exp.x_min, exp.y_max+ 2, 'Contralateral', size = 8)

    ax[0,1].plot(x_data, df.ipsilateral_2s, ls = 'solid', lw = .5, color = 'black', label = 'Ipsilateral') 
    ax[0,1].axis ([exp.x_min, exp.x_max_2, exp.y_min, exp.y_max]) 
    ax[0,1].spines['top'].set_visible (False)
    ax[0,1].spines['right'].set_visible (False)
    ax[0,1].spines['left'].set_visible (False)
    ax[0,1].spines['bottom'].set_visible (False)
    ax[0,1].yaxis.set_ticks([])
    ax[0,1].yaxis.set_ticklabels([])
    ax[0,1].xaxis.set_ticks([])
    ax[0,1].xaxis.set_ticklabels([])

    ax[1,1].plot(x_data, df.contralateral_2s, ls = 'solid', lw = .5, color = 'dimgray', label = 'Contralateral') 
    ax[1,1].axis ([exp.x_min, exp.x_max_2, exp.y_min, exp.y_max]) 
    ax[1,1].spines['top'].set_visible (False)
    ax[1,1].spines['right'].set_visible (False)
    ax[1,1].spines['left'].set_visible (False)
    ax[1,1].spines['bottom'].set_visible (False)
    ax[1,1].yaxis.set_ticks([])
    ax[1,1].yaxis.set_ticklabels([])
    ax[1,1].xaxis.set_ticks([])
    ax[1,1].xaxis.set_ticklabels([])

    #plt.show()
    exp.file_fig = exp.name + '_sEPSP_traces.png' 
    pf.save_fig (exp, fig, exp.file_fig, exp.dir_in + exp.dir_traces)
