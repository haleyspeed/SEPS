import pandas as pd
import os
import numpy as np
import physfiles as pfiles



def get_filtered (data, exp):
    # Filter data based on Vm and Rinput
    # filter [0] = vm1
    # filter [1] = rin1
    data = data.drop(data[data.vm1 > exp.filters[0]].index)
    data = data.drop(data[np.isnan(data.vstim1) == True].index)
    file_out = exp.name + '_coupling_attempted.csv'
    pfiles.save_csv (data, file_out, exp.dir_in, True)
    return data, file_out

def get_clean(data, exp):
    # Remove negative CC values and set to 0
    # vrec1	vrec2	vrec1_nmda	vrec2_nmda
    data['vrec1'] = np.where(data['vrec1'] > 0, 0, data['vrec1'])
    data['vrec2'] = np.where(data['vrec2'] > 0, 0, data['vrec2'])
    data['vrec1_nmda'] = np.where(data['vrec1_nmda'] > 0, 0, data['vrec1_nmda'])
    data['vrec2_nmda'] = np.where(data['vrec2_nmda'] > 0, 0, data['vrec2_nmda'])
    file_out = exp.name + '_coupling_clean.csv'
    pfiles.save_csv (data, file_out, exp.dir_in, True)
    return data

def sort_connected (data, exp):
    data_connected = data.copy()
    for index, row in data_connected.iterrows():
        if row['cc'] < exp.connection:
            data_connected.drop(index, inplace = True)
    file_out = exp.name + '_coupling_connected.csv'
    pfiles.save_csv (data_connected, file_out, exp.dir_in, True)
    return data_connected,file_out

def sort_strong (data, exp):
    data_strong = data.copy()
    for index, row in data_strong.iterrows():
        if row['cc'] >= exp.strong and row['cc_nmda'] >= 0:
            pass
        else:
            data_strong.drop(index, inplace = True)
    file_out = exp.name + '_coupling_strong.csv'
    pfiles.save_csv (data_strong, file_out, exp.dir_in, True)
    return data_strong, file_out

def sort_weak (data, exp):
    data_weak = data.copy()
    for index, row in data_weak.iterrows():
        if row['cc'] < exp.strong and row['cc'] >=0.1 and row['cc_nmda'] >= 0:
            pass
        else:
            data_weak.drop(index, inplace = True)
    file_out = exp.name + '_coupling_weak.csv'
    pfiles.save_csv (data_weak, file_out, exp.dir_in, True)
    return data_weak, file_out

def sort_strengthened (data, exp):
    data_strengthened = data.copy()
    for index, row in data_strengthened.iterrows():
        if row['cc'] >= 0.1 and row['cc_nmda'] >= 0.1 and row['cc_pchange'] > 0:
            pass
        else:
            data_strengthened.drop(index, inplace = True)
    file_out = exp.name + '_coupling_strengthened.csv'
    pfiles.save_csv (data_strengthened, file_out, exp.dir_in, True)
    return data_strengthened, file_out

def sort_weakened (data, exp):
    data_weakened = data.copy()
    for index, row in data_weakened.iterrows():
        if row['cc'] >= 0.1 and row['cc_pchange'] < 0:
            pass
        else:
            data_weakened.drop(index, inplace = True)
    file_out = exp.name + '_coupling_weakened.csv'
    pfiles.save_csv (data_weakened, file_out, exp.dir_in, True)
    return data_weakened, file_out

def sort_NMDA (data, exp):
    data_NMDA = data.copy()
    for index, row in data_NMDA.iterrows():
        if row['cc'] < 0.1 and row['cc_nmda'] > 0.1 and row['cc_pchange'] > 0:
            pass
        else:
            data_NMDA.drop(index, inplace = True)
    file_out = exp.name + '_coupling_created.csv'
    pfiles.save_csv (data_NMDA, file_out, exp.dir_in, True)
    return data_NMDA, file_out

def sort_strength_created (data, exp):
    data_strength_created = data.copy()
    for index, row in data_strength_created.iterrows():
        if row['cc_nmda'] >= 0.1 and row['cc_pchange'] > 0:
            pass
        else:
            data_strength_created.drop(index, inplace = True)
    file_out = exp.name + '_coupling_allstrengthened.csv'
    pfiles.save_csv (data_strength_created, file_out, exp.dir_in, True)
    return data_strength_created, file_out

def sort_nmda5 (data, exp):
    data_nmda5= data.copy()
    for index, row in data_nmda5.iterrows():
        if row['conc_nmda'] != 5:
            data_nmda5.drop(index, inplace = True)
    file_out = exp.name + '_coupling_nmda5.csv'
    pfiles.save_csv (data_nmda5, file_out, exp.dir_in, True)
    return data_nmda5, file_out

def sort_nmda50 (data, exp):
    data_nmda50 = data.copy()
    for index, row in data_nmda50.iterrows():
        if row['conc_nmda'] != 50:
            data_nmda50.drop(index, inplace = True)
    file_out = exp.name + '_coupling_nmda50.csv'
    pfiles.save_csv (data_nmda50, file_out, exp.dir_in, True)
    return data_nmda50, file_out
