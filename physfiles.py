
import os
import pandas as pd
import pickle
from pathlib import Path
import matplotlib.pyplot as plt

# Save the dataframe to a csv
def save_csv (df_out, file_out, dir_out, ignore):
    try:
        os.stat(dir_out)
    except:
        os.makedirs(dir_out)
    # Write data to file
    os.chdir(dir_out)
    df_out.to_csv(file_out, index = ignore)

    # Read a csv file into a dataframe with all data
def get_all (exp): 
    os.chdir(exp.dir_in)
    df_all = pd.read_csv (exp.file_in)
    return df_all

###################################################################
def get_desc (exp, dir_analysis, measure): 
    os.chdir(dir_analysis + exp.dir_measure)
    
    filename = exp.name + '_coupling_' + measure +'_desc.csv' 
    #print ('Loading decriptive stats from : dir_analysis + exp.dir_measure + '//' + filename)
    df_measure = pd.read_csv (filename)
    return df_measure

def set_dirs (list_dirs):
    for dir_out in list_dirs:
        try:
            os.stat(dir_out)
        except:
            os.makedirs(dir_out)
    
def get_groups (exp, df_all): 
    list_all = list ()
    for group in exp.groups:
        df_group = df_all[df_all[exp.factor1] == group]
        list_all.append(df_group)
    return list_all

def set_subdir (exp, subset_name):
    if subset_name == 'attempted':
        dir_in = exp.dir_in + exp.dir_attempted
    elif subset_name == 'connected':
        dir_in = exp.dir_in + exp.dir_connected
    elif subset_name == 'created':
        dir_in = exp.dir_in + exp.dir_created
    elif subset_name == 'strcreated':
        dir_in = exp.dir_in + exp.dir_strcreated
    elif subset_name == 'strengthened':
        dir_in = exp.dir_in + exp.dir_strengthened
    elif subset_name == 'strong':
        dir_in = exp.dir_in + exp.dir_strong
    elif subset_name == 'weak':
        dir_in = exp.dir_in + exp.dir_weak
    elif subset_name == 'weakened':
        dir_in = exp.dir_in + exp.dir_weakened
    return dir_in

def save_class (exp):
    os.chdir(exp.dir_in)
    filename = exp.name + '_' + exp.exp_type + '.p'
    handle = open(filename, 'wb')
    pickle.dump(exp, handle)

def get_anova (exp):
    os.chdir(exp.dir_anova) 
    filename = exp.file_anova
    df_anova = pd.read_csv(filename)
    return df_anova
 
def get_class (filename):
   new_class = pickle.load( open( filename, "rb" ))
   return new_class

def get_x_range(df):
    min = 0.5
    max = 0.5 + len(df.groupby('strain').count())
    x_range = [min, max]
    return x_range 

# Get all files in a directory into a list and convert that list to a str
def get_csv (dir_in):
    dir_current = os.getcwd()
    os.chdir(dir_in)
    list_dfs = list()
    for filename in os.listdir(dir_in):
        if filename.endswith('.csv'):
            new_df = pd.read_csv(filename)
            list_dfs.append(new_df)
    os.chdir(dir_current)
    return list_dfs
def save_fig (exp,fig, file_name, dir_out):
    try:
        os.stat(dir_out)
    except:
        os.makedirs(dir_out)
    os.chdir(dir_out)
    fig.savefig (file_name, dpi= exp.save_dpi, facecolor='w', edgecolor='w',  
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        metadata=None)
