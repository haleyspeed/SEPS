import pandas as pd
import os
import numpy as np
import math
import scipy as sp
from scipy import stats
import physfiles as pfiles

class exp:
    def __init__(self, exp_type, exp_name, groups, dir_in, dir_group, dir_measure,
                dir_tukey, dir_anova, dir_ttest):
        
        self.exp_type = exp_type
        self.name = exp_name
        self.groups = groups
        self.dir_in = dir_in
        self.dir_group = dir_group
        self.dir_measure = dir_measure
        self.dir_tukey = dir_tukey
        self.dir_anova = dir_anova
        self.dir_ttest = dir_ttest

class mini():
    def __init__(self, exp_type, exp_name, groups, factor1, conf, dir_in, dir_group, dir_measure,
                dir_tukey, dir_anova, dir_ttest, isteps, vsteps, sEPSCs, mEPSCs, sEPSPs, mEPSPs, 
                measures_minis, measures_isteps, measures_vsteps, measures_common):
        self.exp_type = exp_type
        self.name = exp_name
        self.groups = groups
        self.factor1 = factor1
        self.conf = conf
        self.dir_in = dir_in
        self.dir_group = dir_group
        self.dir_measure = dir_measure
        self.dir_tukey = dir_tukey
        self.dir_anova = dir_anova
        self.dir_ttest = dir_ttest
        self.dir_isteps = isteps
        self.dir_vsteps = vsteps
        self.dir_sEPSCs = sEPSCs
        self.dir_mEPSCs = mEPSCs
        self.dir_sEPSPs = sEPSPs
        self.dir_mEPSPs = mEPSPs
        self.measures_minis = measures_minis
        self.measures_isteps = measures_isteps
        self.measures_vsteps = measures_vsteps
        self.measures_common = measures_common

# Use to assign the values to an exp, whether it was created locally or from 
# pre-calculated files
def set_exp (exp, exp_type, exp_name, exp_groups, dir_in, file_in, 
                dir_group, dir_measure, dir_tukey, dir_anova, dir_attempted,
                dir_connected, dir_created, dir_nmda5, dir_nmda50, dir_strong,
                dir_strcreated, dir_strengthened, dir_weak, dir_weakened, exp_file, measures, 
                between_factor, confidence_interval, filters, 
                connection, bin_width, strong_connection): 
    exp.exp_type = exp_type

    exp.name = exp_name
    exp.groups = exp_groups
    exp.dir_in = dir_in
    exp.file_in = file_in
    exp.dir_group = dir_group
    exp.dir_measure = dir_measure
    exp.dir_tukey = dir_tukey
    exp.dir_anova = dir_anova
    exp.dir_attempted = dir_attempted
    exp.dir_connected = dir_connected
    exp.dir_created = dir_created
    exp.dir_nmda5 = dir_nmda5
    exp.dir_nmda50 = dir_nmda50
    exp.dir_strong = dir_strong
    exp.dir_strcreated = dir_strcreated
    exp.dir_strengthened = dir_strengthened
    exp.dir_weak = dir_weak
    exp.dir_weakened = dir_weakened 
    exp.exp_file = exp_file
    exp.measures = measures
    exp.factor1 = between_factor
    exp.conf = confidence_interval
    exp.filters = filters
    exp.connection = connection
    exp.bin = bin_width
    exp.strong = strong_connection
    return exp


# Calculates Rin from CC data (not from VSteps to be consistent with xfer)
# Validated in Origin
def get_rin(data_Rin):
    data_Rin['rin1'] = np.nan
    data_Rin['rin2'] = np.nan
    data_Rin['rin1_nmda'] = np.nan
    data_Rin['rin2_nmda'] = np.nan

    for index, row in data_Rin.iterrows():	
        if row['istim1'] != 0 and math.isnan(row['istim1']) == False:
            data_Rin.loc[index:index:, 'rin1'] = row['vstim1'] / row['istim1'] * 1000
        if row['istim2'] != 0 and math.isnan(row['istim2']) == False:
            data_Rin.loc[index:index:, 'rin2'] = row['vstim2'] / row['istim2'] * 1000
        if row['istim1_nmda'] != 0 and math.isnan(row['vstim1_nmda']) == False:
            data_Rin.loc[index:index:, 'rin1_nmda'] = row['vstim1_nmda'] / row['istim1_nmda'] * 1000    
        if row['istim2_nmda'] != 0 and math.isnan(row['istim2_nmda']) == False:
            data_Rin.loc[index:index:, 'rin2_nmda'] = row['vstim2_nmda'] / row['istim2_nmda'] * 1000
    return data_Rin

# Combine 1>2 and 2>1
# Validated in Origin
def get_connections (data_connect):
    connections = pd.DataFrame({'id':[],'strain':[],'distance':[], 'conc_nmda': [], 'vm1':[],'vm2':[],
        'vrec1':[],'vstim1':[],'istim1':[],'vm1_nmda':[],'vm2_nmda':[],'vrec2_nmda':[],'vstim1_nmda':[],'istim1_nmda':[],
        'rin1':[],'rin2':[],'rin1_nmda':[],'rin2_nmda':[]}) 
    for index, row in data_connect.iterrows():
        connections = connections.append({'id':row['id'],'strain':row['strain'],'distance':row['distance'], 'conc_nmda': row['conc_nmda'],
            'vrec1':row['vrec1'],'vstim1':row['vstim1'],'istim1':row['istim1'],
            'vrec1_nmda':row['vrec1_nmda'],'vstim1_nmda':row['vstim1_nmda'],'istim1_nmda':row['istim1_nmda'],
            'rin1':row['rin1'],'rin2':row['rin2'],'rin1_nmda':row['rin1_nmda'],'rin2_nmda':row['rin2_nmda'],
            'vm1':row['vm1'],'vm1_nmda':row['vm1_nmda']}, ignore_index = True) 
        connections = connections.append({'id':row['id'],'strain':row['strain'],'distance':row['distance'],'conc_nmda': row['conc_nmda'],
            'vrec1':row['vrec2'],'vstim1':row['vstim2'],'istim1':row['istim2'],
            'vrec1_nmda':row['vrec2_nmda'],'vstim1_nmda':row['vstim2_nmda'],'istim1_nmda':row['istim2_nmda'],
            'rin1':row['rin2'],'rin2':row['rin1'],'rin1_nmda':row['rin2_nmda'],'rin2_nmda':row['rin1_nmda'], 'vm1':row['vm2'],'vm1_nmda':row['vm2_nmda']}, ignore_index = True)
        connections = connections[['id','strain','distance','conc_nmda', 'vrec1','vstim1','istim1',
            'vrec1_nmda','vstim1_nmda','istim1_nmda','rin1','rin2','rin1_nmda','rin2_nmda', 'vm1','vm2','vm1_nmda', 'vm2_nmda']] 
    return connections

# Get CC
# Validated in Origin
def get_cc(data_cc):
    data_cc['cc'] = np.nan
    data_cc['cc_nmda'] = np.nan
    
    # Checking for NaN in vrec rather than vstim because vrec appears first in the calculation
    # Otherwise you would get a type error on the next line because vrec is called first
    for index, row in data_cc.iterrows():
        if row['vstim1'] != 0 and math.isnan(row['vstim1']) == False: 
            data_cc.loc[index:index:, 'cc'] = row['vrec1'] / row['vstim1'] * 100
        if row['vstim1_nmda'] != 0 and math.isnan(row['vrec1_nmda']) == False:
            data_cc.loc[index:index:, 'cc_nmda'] = row['vrec1_nmda'] / row['vstim1_nmda'] * 100
    return data_cc

# Get xfer
# Validated in Origin
def get_xfer(data_xfer):
    data_xfer['xfer'] = np.nan
    data_xfer['xfer_nmda'] = np.nan
    for index, row in data_xfer.iterrows():	
        if row['istim1'] != 0 and math.isnan(row['vrec1']) == False:
            data_xfer.loc[index:index:, 'xfer'] = row['vrec1'] / row['istim1'] * 1000
        if row['istim1_nmda'] != 0 and math.isnan(row['vrec1_nmda']) == False:
            data_xfer.loc[index:index:, 'xfer_nmda'] = row['vrec1_nmda'] / row['istim1_nmda'] * 1000
    return data_xfer

# Get gj
# Validated in Origin
def get_gj(data_gj):
    data_gj['gj'] = np.nan
    data_gj['gj_nmda'] = np.nan
    for index, row in data_gj.iterrows():	
        if row['xfer'] != 0:
            data_gj.loc[index:index:, 'gj'] = 1/(((row['rin1'] * row['rin2'])-(row['xfer']*row['xfer']))/row['xfer']) * 1000000
        if row['xfer_nmda'] != 0:
            data_gj.loc[index:index:, 'gj_nmda'] = 1/(((row['rin1_nmda'] * row['rin2_nmda'])-(row['xfer_nmda']*row['xfer_nmda']))/row['xfer_nmda']) * 1000000
    return data_gj

# Get_change in gj, cc, and Rin
# Validated in Origin
def get_changes(data_change):
    data_change['cc_change'] = np.nan
    data_change['cc_pchange'] = np.nan
    data_change['rin_change'] = np.nan
    data_change['rin_pchange'] = np.nan
    data_change['gj_change'] = np.nan
    data_change['gj_pchange'] = np.nan
    data_change['vm_change'] = np.nan
    for index, row in data_change.iterrows():
        if row['cc'] != 0:
            data_change.loc[index:index:, 'cc_change'] = row['cc_nmda']-row['cc']
            data_change.loc[index:index:,'cc_pchange'] = ((row['cc_nmda']-row['cc'])/row['cc']) * 100
        if row['rin1'] != 0:
            data_change.loc[index:index:,'rin_change'] = row['rin1_nmda']-row['rin1']
            data_change.loc[index:index:,'rin_pchange'] = ((row['rin1_nmda']-row['rin1'])/row['rin1']) * 100
        if row['gj'] != 0:
            data_change.loc[index:index:,'gj_change'] = row['gj_nmda']-row['gj']
            data_change.loc[index:index:,'gj_pchange'] = ((row['gj_nmda']-row['gj'])/row['gj']) * 100
        data_change.loc[index:index:,'vm_change'] = row['vm1_nmda']-row['vm1']
    return data_change

def get_bin (data, exp):
    splitter = exp.bin
    data_bin = data.copy()
    data_bin['bins'] = np.nan
    for index, row in data_bin.iterrows():
        if row['distance'] >= 0 and row['distance'] < splitter:
            data_bin.loc[index:index:, 'bins'] = splitter/2
        elif row['distance'] >= splitter and row['distance'] < 2*splitter:
            data_bin.loc[index:index:, 'bins'] = 2*splitter-(splitter/2)
        elif row['distance'] >= 2*splitter and row['distance'] < 3*splitter:
            data_bin.loc[index:index:, 'bins'] = 3*splitter - (2*splitter/2)
        elif row['distance'] >= 3*splitter and row['distance'] < 4*splitter:
            data_bin.loc[index:index:, 'bins'] = 4*splitter - (3*splitter/2)
        elif row['distance'] >= 4*splitter and row['distance'] < 5*splitter:
            data_bin.loc[index:index:, 'bins'] = 5*splitter - (4*splitter/2)
        elif row['distance'] >= 5*splitter and row['distance'] < 6*splitter:
            data_bin.loc[index:index:, 'bins'] = 6*splitter - (5*splitter/2)
    return data_bin


def get_calc (data, exp):
    data = get_bin (data, exp)
    data = get_cc (data)
    data = get_xfer (data)
    data = get_gj (data)
    data = get_changes (data)
    return data

def round_down(data, round_to):
    rounded = data - (data % round_to)
    return rounded

def round_up(data, round_to):
    rounded = data + (round_to - data % round_to)
    return rounded

