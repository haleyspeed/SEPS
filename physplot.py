import math
import numpy as np
import os
from scipy.stats import sem
from scipy import stats
import matplotlib.pyplot as plt
import physfiles as pf

class scatter:
    def __init__(self):
        self.x_min = float()
        self.x_max = float()
        self.x_label_text = list() 
        self.x_label_size = 12
        self.y_label_size = 14
        self.x_ticklabel_angle = 45
        self.x_title_size = 16
        self.y_title_size = 16
        self.x_title_text = ''
        self.y_title_text = str()
        self.y_label_text = str()
        self.y_min = float()
        self.y_max = float()
        self.hline_y = 0
        self.hline_text = 12
        self.hline_loc = float()
        self.hline_color = 'silver'
        self.hline_size = 1
        self.hline_style = ':'
        self.scatter_marker_size = 2
        self.scatter_marker_colors = ['dimgray', 'teal', 'blueviolet', 'black', ]
        self.scatter_plot_markers = ['s','o','v','_']
        self.avg_marker_size = 50
        self.avg_marker_color = 'black'
        self.avg_err_size = 5
        self.avg_err_color = 'black'
        self.err_cap_size = 5
        self.err_line_width = 2
        self.anova_text = str()
        self.anova_stars = ''
        self.anova_text_size = 8
        self.sig_between = ''
        self.factor_between = 'Genotype'
        self.facotr_within = 'Treatment'
        self.sig_within = ''
        self.posthoc_stars1 = ['','','']
        self.posthoc_stars1_pos = list()
        self.posthoc_stars1_size = 12
        self.posthoc_stars2 = ['','','']
        self.posthoc_stars2_pos = list()
        self.posthoc_stars2_size = 12  
        self.posthoc_stars3 = ['','','']
        self.posthoc_stars3_pos = list()
        self.posthoc_stars3_size = 12     
        self.plot_dpi = 80     
        self.save_dpi = 300 

def set_scatter_x (scatter, x_data, x_label_text):
    scatter.x_max = scatter.x_min + len(x_data)
    scatter.x_label_text = x_label_text
    scatter.x_min = -.5
    scatter.x_max = len(x_label_text) - .5
    return scatter

def set_ylim (y_range):
    i = 1
    while i < 2000:
        if y_range < 5 * i:
            return 5
        i += 5

def set_scatter_anova (scatter, df, measure): #, anova_stars, posthoc_stars1, posthoc_stars1_pos): 
    p = df.loc[measure, 'p_value']     
    if p < 0.00001:
        scatter.anova_text =  'Genotype  ' +  '****'
    elif p < 0.001:
        scatter.anova_text =  'Genotype  ' +  '***'
    elif p < 0.01: 
        scatter.anova_text =  'Genotype  ' +  '**'
    elif p < 0.05:
        scatter.anova_text =  'Genotype  ' +  '*'
    elif p >= 0.05:
        scatter.anova_text =  'Genotype  ' +  'n.s.'      


def get_format_coupling (scatter, measure):
    if measure == 'cc':
        scatter.y_min = -0.5
        scatter.y_max = 5
        scatter.hline_y = 0.1
        scatter.hline_text = '0.1%'
        scatter.y_label_text = 'Coupling Coefficient (%)'
    elif measure == 'cc_nmda':
        scatter.y_min = -0.2
        scatter.y_max = 5
        scatter.hline_y = 0.2
        scatter.hline_text = '0.1%'
        scatter.y_label_text = 'Coupling Coefficient (%) \n with NMDA'
    elif measure == 'cc_change':
        scatter.y_min = -3
        scatter.y_max = 3
        scatter.hline_y = 0
        scatter.hline_text = ''
        scatter.y_label_text = '\u0394 Coupling Coefficient \n with NMDA'
    elif measure == 'cc_pchange':
        scatter.y_min = -200
        scatter.y_max = 500
        scatter.hline_y = 0
        scatter.hline_text = ''
        scatter.y_label_text = '% \u0394 Coupling Coefficient (%) \n with NMDA'
    elif measure == 'gj':
        scatter.y_min = -100
        scatter.y_max = 300
        scatter.hline_y = 0
        scatter.hline_text = ''
        scatter.y_label_text = 'Junctional Conductance (pS) '
    elif measure == 'gj_nmda':
        scatter.y_min = -100
        scatter.y_max = 300
        scatter.hline_y = 0
        scatter.hline_text = ''
        scatter.y_label_text = 'Junctional Conductance (pS)\n with NMDA '
    elif measure == 'gj_change':
        scatter.y_min = -100
        scatter.y_max = 200
        scatter.hline_y = 0
        scatter.hline_text = ''
        scatter.y_label_text = '\u0394 Junctional Conductance (pS) \n with NMDA'
    elif measure == 'gj_pchange':
        scatter.y_min = -200
        scatter.y_max = 1000
        scatter.hline_y = 0
        scatter.hline_text = ''
        scatter.y_label_text = '% \u0394 Junctional Conductance \n with NMDA'
    elif measure == 'rin1':
        scatter.y_min = -100
        scatter.y_max = 500
        scatter.hline_y = 0
        scatter.hline_text = ''
        scatter.y_label_text = 'Input Resistance (M\u03A9)'
    elif measure == 'rin1_nmda':
        scatter.y_min = -100
        scatter.y_max = 500
        scatter.hline_y = 0
        scatter.hline_text = ''
        scatter.y_label_text = 'Input Resistance (M\u03A9) \n with NMDA'
    elif measure == 'rin_change':
        scatter.y_min = -500
        scatter.y_max = 500
        scatter.hline_y = 0
        scatter.hline_text = ''
        scatter.y_label_text = '\u0394 Input Resistance (M\u03A9)  \n with NMDA'
    elif measure == 'rin_pchange':
        scatter.y_min = -200
        scatter.y_max = 200
        scatter.hline_y = 0
        scatter.hline_text = ''
        scatter.y_label_text = '% \u0394 Input Resistance (M\u03A9) \n with NMDA'
    elif measure == 'vm1':
        scatter.y_min = -100
        scatter.y_max = 0
        scatter.hline_y = -40
        scatter.hline_text = ''
        scatter.y_label_text = 'Resting Membrane Potential (mV)'
    elif measure == 'vm_change':
        scatter.y_min = -50
        scatter.y_max = 50
        scatter.hline_y = 0
        scatter.hline_text = ''
        scatter.y_label_text = '\u0394 Resting Membrane Potential (mV) \n with NMDA'


def set_stars_posthoc (exp, df, compare_to): # needs to stand alone, no assoc with scatter
    # how many groups? 2,3,4,6,8
    # if < 4 compare to 1 (wt) [*]
    # if < 6 compare to 2 (wt, control) [#]
    # if < 8 compare to 3 (wt, control, hemisphere) [@]
    if len (df['group1']) == 2: # for later when I have t-test data make sure stats file has same structure as dunnett/tukey
        stars1, stars1_pos = set_stars(df.iloc[0,4])
        if compare_to == exp.exp_groups[0]:
            exp.posthoc_stars1_pos += 1
    if len (df['group1']) == 3:
        exp.posthoc_stars1, exp.posthoc_stars1_pos = set_stars(df.iloc[0,4])
        stars2, stars2_pos = set_stars(df.iloc[1,4])
        if compare_to == exp.exp_groups[0]:
            exp.posthoc_stars1_pos += 1
            exp.posthoc_stars2_pos += 2
        elif compare_to == exp.exp_groups[1]:
            exp.posthoc_stars1_pos += 0
            exp.posthoc_stars2_pos += 2
        elif compare_to == exp.exp_groups[2]:
            exp.posthoc_stars1_pos += 0
            exp.posthoc_stars2_pos += 1

    if len (df['group1']) == 4 and len(compare_to) == 1:
        exp.posthoc_stars1, exp.posthoc_stars1_pos = set_stars(df.iloc[0,4])
        exp.posthoc_stars2, exp.posthoc_stars2_pos = set_stars(df.iloc[1,4])
        exp.posthoc_stars3, exp.posthoc_stars3_pos = set_stars(df.iloc[2,4])
        if compare_to == exp.exp_groups[0]:
            exp.posthoc_stars1_pos += 1
            exp.posthoc_stars2_pos += 2
            exp.posthoc_stars3_pos += 3
        elif compare_to == exp.exp_groups[1]:
            exp.posthoc_stars1_pos += 0
            exp.posthoc_stars2_pos += 2
            exp.posthoc_stars2_pos += 3
        elif compare_to == exp.exp_groups[2]:
            exp.posthoc_stars1_pos += 0
            exp.posthoc_stars2_pos += 1
            exp.posthoc_stars3_pos += 3
        elif compare_to == exp.exp_groups[3]:
            exp.posthoc_stars1_pos += 0
            exp.posthoc_stars2_pos += 1
            exp.posthoc_stars3_pos += 2
    # Add function here to actually plot the stars

def get_t_stars (exp, df_ttest, measure):    
        if df_ttest.loc[measure]['p_value'] < 0.0001:
            exp.posthoc_stars = '****'
            exp.posthoc_stars_pos = '1.5' 
        elif df_ttest.loc[measure]['p_value'] < 0.001:
            exp.posthoc_stars = '***'
            exp.posthoc_stars_pos = '1.5'
        elif df_ttest.loc[measure]['p_value'] < 0.01:
            exp.posthoc_stars = '**'
            exp.posthoc_stars_pos = '1.5'
        elif df_ttest.loc[measure]['p_value'] < 0.05:
            exp.posthoc_stars = '*'
            exp.posthoc_stars_pos = '1.5'
        elif df_ttest.loc[measure]['p_value'] > 0.05:
            exp.posthoc_stars = ''
            exp.posthoc_stars_pos = '1.5'
        

def save_fig (exp, sct, fig):
    try:                         
        os.stat(exp.dir_plots)
    except:
        os.mkdir(exp.dir_plots)
    
    # Change the current directory to the output directory
    os.chdir(exp.dir_plots)

    # Write the figure to file
    fig.savefig(sct.file_fig, dpi=sct.save_dpi, facecolor='w', edgecolor='w',  
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        metadata=None)

    # Reset the working directory to your input directory
    os.chdir(exp.dir_plots)

def make_plot (exp, sct, fig, ax, measure, df_desc, list_groupdf):
    
    ax.spines['top'].set_visible(False) 
    ax.spines['right'].set_visible(False)
   
    x_data = list()
    y_mean = list()
    y_err = list()
    i = 0
    for df in list_groupdf:
        # Plot the scatter data
        ax.plot(df[exp.factor1], df[measure], linestyle = 'none', marker = sct.scatter_plot_markers [i], color = sct.scatter_marker_colors[i], markersize = sct.scatter_marker_size) # plot first scatter dataset     
        x_data.append(df.iloc[1][exp.factor1])
        y_mean.append(df_desc.iloc[1][i+1])
        y_err.append(df_desc.iloc[3][i+1])
        i += 1
    
    get_format_coupling(sct, measure)    
    ax.plot(x_data, y_mean, linestyle = 'none', marker = '_', color = 'black', markersize = sct.avg_marker_size)    # plot first mean
    ax.errorbar(x_data, y_mean, yerr = y_err, linestyle = 'none', elinewidth = 2, capsize = 5, color = 'black' )
    ax.axis ([sct.x_min, sct.x_max, sct.y_min, sct.y_max]) 
    # Draw the horizontal line (line y-value is determined in get_measure())
    ax.axhline(sct.hline_y, 0, 1 , color = sct.hline_color, lw = sct.hline_size, linestyle = sct.hline_style) # Add a horizontal line
    ax.set_ylabel(sct.y_label_text, fontsize = sct.y_label_size)    
    ax.set_xticklabels (sct.x_label_text, fontsize = sct.x_label_size, rotation = sct.x_ticklabel_angle)
    
    # ANOVA significance
    ax.text(-0.4, sct.y_max, sct.anova_text, size = sct.anova_text_size, color = 'black')


    sct.file_fig = exp.name + '_' + exp.exp_type + '_' + measure + '_scatter.png' 
    save_fig (exp, sct, fig)
    return ax

def mini_scatter (exp, df_group1, df_group2, df_scatter, df_ttest, measure):
    x_min = -.5
    x_data = exp.groups
    x_max = x_min + len(x_data) 
    # for x in x_data sort x, if x1-x < x *0.1 then x = x + x*.01   
    y_mean = list([abs(df_group1.loc['mean', measure]), abs(df_group2.loc['mean', measure])])
    y_err = list([abs(df_group1.loc['se', measure]), abs(df_group2.loc['se', measure])])      
    y_data =  abs(df_scatter[measure])
    y_min = 0
    y_max = max(y_data)
    y_range = y_max-y_min
    round_to = set_ylim (y_range)
    y_int = math.ceil(y_range/round_to)
    y_max = (y_int) * round_to      
    df_y1 =  df_scatter[df_scatter.lateral != 'contralateral']
    df_y2 =  df_scatter[df_scatter.lateral != 'ipsilateral']
    get_t_stars(exp, df_ttest, measure)

    # Plot the data
    fig, ax = plt.subplots(figsize=(3, 4), dpi = exp.save_dpi)
    ax.set_prop_cycle(color = exp.marker_colors, marker = exp.markers)
    print(measure, df_scatter.iloc[0]['exp'])
    ax.plot(df_y1['lateral'], abs(df_y1[measure]) , linestyle = 'none', markersize = exp.scatter_marker_size)
    ax.plot(df_y2['lateral'], abs(df_y2[measure]), linestyle = 'none', markersize = exp.scatter_marker_size) 
    ax.plot(x_data, y_mean, linestyle = 'none', marker = exp.avg_marker, color = 'black', markersize = exp.avg_marker_size)    # plot first mean
    ax.errorbar(x_data, y_mean, yerr = y_err, linestyle = 'none', elinewidth = exp.err_line_width, capsize = 5, color = 'black' )
                
    # Common formatting variables
    ax.axis ([x_min, x_max, y_min, y_max]) 
    ax.spines['top'].set_visible(False) 
    ax.spines['right'].set_visible(False)
    ax.set_ylabel(exp.y_label_text, fontsize = exp.y_label_size)    
    ax.set_xticklabels (exp.x_label_text, fontsize = exp.x_label_size, rotation = exp.x_ticklabel_angle)

    # Add significance and text labels
    get_t_stars (exp, df_ttest, measure) 
    
    plt.gcf().text(0.2, -0.025, '(' + str(int(df_group1.iloc[0][measure])) + ')', fontsize=10, color = 'black', rotation = exp.x_ticklabel_angle)
    plt.gcf().text(0.55, -0.05, '(' + str(int(df_group2.iloc[0][measure])) + ')', fontsize=10, color = 'black', rotation = exp.x_ticklabel_angle)
    ax.set_xticklabels (exp.group_labels, rotation = exp.x_ticklabel_angle)
    plt.setp (ax.xaxis.get_majorticklabels(), ha='right')
    exp.file_fig = exp.name + '_' + df_scatter.iloc[0]['exp'] + '_' + measure + '_scatter.png' 
    pf.save_fig (exp, fig, exp.file_fig, exp.dir_plots)

def iv_plot (exp, df_group1, df_group2):
    # Get rid of unnecessary columns
    cols = df_group1['stat']
    df_group1 = df_group1[exp.measures_vsteps]
    df_group2 = df_group2[exp.measures_vsteps]
    
    # Transpose columns
    df_trans1 = df_group1.transpose()
    df_trans2 = df_group2.transpose()
    df_trans1.columns = cols
    df_trans2.columns = cols
    df_trans1.reset_index()
    
    # Format x column to numeric
    x_data = list()
    x_data1 = df_trans1.index.values
    for x in x_data1:
        x = str(x).replace('vm', '')
        x_data.append(int(x))

    # Assign y_data to individual variables for linear regression
    y1_data = df_trans1['mean']
    y2_data = df_trans2['mean']
    y1_err = df_trans1['se']
    y2_err = df_trans2['se']
    x_min = -80
    x_max = 80
    y_min = -600
    y_max = 800
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x_data,y1_data)
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(x_data,y2_data)   
    label1 = 'Ipsilateral'
    label2 = 'Contralateral'
    list_data = [y1_data, y2_data]
    list_err = [y1_err, y2_err]
    list_slope = [slope1, slope2]
    list_intercept = [intercept1, intercept2]
    list_labels = [label1, label2]

    # Initiate the figure
    fig, ax = plt.subplots(figsize=(4, 4), dpi = exp.save_dpi) 
    
    
    m = 0
    while m < len(list_data): 
        # Set origin
        ax = plt.gca()  
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.spines['bottom'].set_position(('data',0))
        ax.yaxis.set_ticks_position('left')
        ax.spines['left'].set_position(('data',0))
    
        # Plot the data
        ax.scatter(x_data, list_data[m], color = exp.marker_colors[m], s = 10, label = list_labels[m])  # plot first mean
        ax.errorbar(x_data, list_data[m], yerr = list_err[m], linestyle = 'none', elinewidth = 1, capsize = 5, color = exp.marker_colors[m]) # plot error bars
        line = [i * list_slope[m] + list_intercept[m] for i in x_data]
        ax.plot (x_data, line, color = exp.marker_colors[m])
        ax.axis ([x_min, x_max, y_min, y_max]) 
        
        m += 1
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels[0:2], frameon = False, loc='upper left', bbox_to_anchor=(-.15, 1.2), shadow=False, fontsize = 10)
    #plt.gcf().text(0.0, 0.4, 'mV', fontsize=10, color = 'black')
    plt.gcf().text(0.35, .11, '-600 pA', fontsize=10, color = 'black')
    plt.gcf().text(0.35, .875, '800 pA', fontsize=10, color = 'black')
    plt.gcf().text(0.05, 0.4, '-80 mV', fontsize=10, color = 'black')
    plt.gcf().text(0.85, 0.4, '80 mV', fontsize=10, color = 'black')
    
    ax.xaxis.set_major_locator(plt.MultipleLocator(20))
    x_lab = ax.get_xticklabels()
    ax.set_xticklabels(x_lab)
    ax.yaxis.set_major_locator(plt.MultipleLocator(200))
    y_lab = ax.get_yticklabels()
    y_lab[:] = ' '
    ax.set_yticklabels(y_lab)
    exp.file_fig = exp.name + '_iv_plot.png' 
    pf.save_fig (exp, fig, exp.file_fig, exp.dir_plots)
