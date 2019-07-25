import pandas as pd
import numpy as np
import scipy as sp
from scipy import stats
import math
import os
import physfiles as pfiles
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison
from statsmodels.stats.libqsturng import psturng

def get_qualitative (df_all, df_connected, df_created, 
            df_strcreated, df_strong, df_weak, df_weakened, exp):
    
    total = df_all.groupby('strain').count()
    total['category'] = 'total'
    connected = df_connected.groupby('strain').count()
    connected['category'] = 'connected'
    created = df_created.groupby('strain').count()
    created['category'] = 'created'
    strcreated = df_strcreated.groupby('strain').count()
    strcreated['category'] = 'strcreated'
    strong = df_strong.groupby('strain').count()
    strong['category'] = 'strong'
    weak = df_weak.groupby('strain').count()
    weak['category'] = 'weak'
    weakened = df_weakened.groupby('strain').count()
    weakened['category'] = 'weakened'
    list_quals = list([total, connected, created, strcreated, strong,
                    weak, weakened])
    
    # Separate by genotype
    cols = ['strain', 'category', 'n', 'total', 'percent']
    df_cat = pd.DataFrame(columns = cols)
    for df in list_quals:
        df.reset_index(level=df.index.names, inplace=True)
        df = df.rename(columns={ 'index': 'strain'})
        for index, row in df.iterrows():
            for group in exp.groups:
                if row ['strain'] == group:
                    tot = total.loc[index,'vrec1']
                    percent = row['vrec1']/tot * 100
                    new_row = [{'strain': group, 'category': row['category'], 'n': row['vrec1'],
                                'total': tot, 'percent': percent}]
                    df_cat = df_cat.append(new_row) 
    return df_cat

def get_desc_cc (exp, df_desc, dir_out):
    factor = exp.factor1
    confidence = exp.conf
    if df_desc.empty:
        # Make new dataframe for stats
        df_return = pd.DataFrame(columns = ['stat','strain','id','distance','conc_nmda',
                'vrec1', 'vstim1', 'istim1','vrec1_nmda','vstim1_nmda','istim1_nmda','rin1','rin2',
                'rin1_nmda','rin2_nmda','vm1',	'vm2',	'vm1_nmda',	'vm2_nmda',	'bins',	'cc','cc_nmda',	
                'xfer',	'xfer_nmda','gj','gj_nmda','cc_change',	'cc_pchange',	'rin_change',	'rin_pchange',
                'gj_change', 'gj_pchange',	'vm_change'])
        df_return['stat'] = ['n', 'mean', 'sd', 'se', 'sum', 'min','max','quart_25', 'quart_75','median','ci_5', 'ci_95']
        df_return.index = ['n', 'mean', 'sd', 'se', 'sum', 'min','max','quart_25', 'quart_75','median','ci_5', 'ci_95']
        df_return.reset_index(level=df_return.index.names, inplace=True)
        df_return = df_return.rename(columns={ 'index': 'stat'})
        
    else: 
        group = df_desc.iloc[0][factor]
        # Calculations for all measures at once
        n = df_desc.groupby(factor).count()
        avg = df_desc.groupby(factor).mean()
        sd = df_desc.groupby(factor).std()
        se = sd/np.sqrt(n.astype('int'))
        added = df_desc.groupby(factor).sum()
        minimum = df_desc.groupby(factor).min()
        maximum = df_desc.groupby(factor).max()
        quartile25 = df_desc.groupby(factor).quantile(q = 0.25, 
                axis = 0, numeric_only = True, interpolation = 'linear')
        quartile75 = df_desc.groupby(factor).quantile(q = 0.75, 
                axis = 0, numeric_only = True, interpolation = 'linear')
        median = df_desc.groupby(factor).median()
        conf = se * sp.stats.t._ppf((1+confidence)/2., n.astype('int')-1)
        conf_5 = avg - conf
        conf_95 = avg + conf
        n = n.reset_index()
        avg = avg.reset_index()
        sd = sd.reset_index()
        added = added.reset_index()
        minimum = minimum.reset_index()
        maximum = maximum.reset_index()
        quartile25 = quartile25.reset_index()
        quartile75 = quartile75.reset_index()
        median = median.reset_index()

        # Make new dataframe for stats
        df_return = pd.DataFrame()
        df_return = pd.concat([n,avg,sd,se,added,minimum,maximum,quartile25,quartile75,median,conf_5,conf_95], sort = False, ignore_index=False)
        df_return.index = ['n', 'mean', 'sd', 'se', 'sum', 'min','max','quart_25', 'quart_75','median','ci_5', 'ci_95']
        df_return.reset_index(level=df_return.index.names, inplace=True)
        df_return = df_return.rename(columns={ 'index': 'stat'})

        file_out = exp.name + '_coupling_' + group + '_desc.csv'
        pfiles.save_csv (df_return, file_out, dir_out, False)
    return df_return

def get_desc_mini (exp, df_desc):
    factor = exp.factor1
    confidence = exp.conf
    exp_type = df_desc.iloc[0]['exp']
    df_desc = df_desc [['lateral', 'genotype', 'amp', 'q', 'half_width', 'decay', 'freq']]
    if df_desc.empty:
        # Make new dataframe for stats
        df_return = pd.DataFrame(columns = ['stat','exp','file','genotype','slice',
                'shaved', 'hemisphere', 'lateral','amp','q','half_width','decay','freq'])
        df_return['stat'] = ['n', 'mean', 'sd', 'se', 'sum', 'min','max','quart_25', 'quart_75','median','ci_5', 'ci_95']
        df_return.index = ['n', 'mean', 'sd', 'se', 'sum', 'min','max','quart_25', 'quart_75','median','ci_5', 'ci_95']
        df_return.reset_index(level=df_return.index.names, inplace=True)
        df_return = df_return.rename(columns={ 'index': 'stat'})
        
    else: 
        group = df_desc.iloc[0][factor]
        # Calculations for all measures at once
        n = df_desc.groupby(factor).count()
        avg = df_desc.groupby(factor).mean()
        sd = df_desc.groupby(factor).std()
        se = df_desc.groupby(factor).sem()
        added = df_desc.groupby(factor).sum()
        minimum = df_desc.groupby(factor).min()
        maximum = df_desc.groupby(factor).max()
        quartile25 = df_desc.groupby(factor).quantile(q = 0.25, 
                axis = 0, numeric_only = True, interpolation = 'linear')
        quartile75 = df_desc.groupby(factor).quantile(q = 0.75, 
                axis = 0, numeric_only = True, interpolation = 'linear')
        median = df_desc.groupby(factor).median()
        conf = se * sp.stats.t._ppf((1+confidence)/2., n.astype('int')-1)
        conf_5 = avg - conf
        conf_95 = avg + conf
        n = n.reset_index()
        avg = avg.reset_index()
        sd = sd.reset_index()
        added = added.reset_index()
        minimum = minimum.reset_index()
        maximum = maximum.reset_index()
        quartile25 = quartile25.reset_index()
        quartile75 = quartile75.reset_index()
        median = median.reset_index()

        # Make new dataframe for stats
        df_return = pd.DataFrame()
        df_return = pd.concat([n,avg,sd,se,added,minimum,maximum,quartile25,quartile75,median,conf_5,conf_95], sort = False, ignore_index=True)
        df_return.index = ['n', 'mean', 'sd', 'se', 'sum', 'min','max','quart_25', 'quart_75','median','ci_5', 'ci_95']
        df_return.reset_index(level=df_return.index.names, inplace=True)
        df_return = df_return.rename(columns={ 'index': 'stat'})

        file_out = exp.name + '_' + exp_type + '_' + group + '_desc.csv'
        pfiles.save_csv (df_return, file_out, exp.dir_out, False)
    return df_return

def get_desc_isteps (exp, df_desc):
    factor = exp.factor1
    confidence = exp.conf
    exp_type = df_desc.iloc[0]['exp']
    
    df_desc = df_desc [['genotype', 'lateral', 'vm', 'holding-150', 'holding-125', 'holding-100', 'holding-75', 'holding-50', 
                        'holding-25', 'holding', 'holding+25', 'holding+50', 'holding+75']]
    if df_desc.empty:
        # Make new dataframe for stats
        df_return = pd.DataFrame(columns = ['stat','exp','genotype', 'lateral', 'vm', 'holding-150', 'holding-125', 'holding-100', 'holding-75', 'holding-50', 
                        'holding-25', 'holding', 'holding+25', 'holding+50', 'holding+75'])
        df_return['stat'] = ['n', 'mean', 'sd', 'se', 'sum', 'min','max','quart_25', 'quart_75','median','ci_5', 'ci_95']
        df_return.index = ['n', 'mean', 'sd', 'se', 'sum', 'min','max','quart_25', 'quart_75','median','ci_5', 'ci_95']
        df_return.reset_index(level=df_return.index.names, inplace=True)
        df_return = df_return.rename(columns={ 'index': 'stat'})
        
    else: 
        group = df_desc.iloc[0][factor]
        # Calculations for all measures at once
        n = df_desc.groupby(factor).count()
        avg = df_desc.groupby(factor).mean()
        sd = df_desc.groupby(factor).std()
        se = df_desc.groupby(factor).sem()
        added = df_desc.groupby(factor).sum()
        minimum = df_desc.groupby(factor).min()
        maximum = df_desc.groupby(factor).max()
        quartile25 = df_desc.groupby(factor).quantile(q = 0.25, 
                axis = 0, numeric_only = True, interpolation = 'linear')
        quartile75 = df_desc.groupby(factor).quantile(q = 0.75, 
                axis = 0, numeric_only = True, interpolation = 'linear')
        median = df_desc.groupby(factor).median()
        conf = se * sp.stats.t._ppf((1+confidence)/2., n.astype('int')-1)
        conf_5 = avg - conf
        conf_95 = avg + conf
        n = n.reset_index()
        avg = avg.reset_index()
        sd = sd.reset_index()
        added = added.reset_index()
        minimum = minimum.reset_index()
        maximum = maximum.reset_index()
        quartile25 = quartile25.reset_index()
        quartile75 = quartile75.reset_index()
        median = median.reset_index()

        # Make new dataframe for stats
        df_return = pd.DataFrame()
        df_return = pd.concat([n,avg,sd,se,added,minimum,maximum,quartile25,quartile75,median,conf_5,conf_95], sort = False, ignore_index=True)
        df_return.index = ['n', 'mean', 'sd', 'se', 'sum', 'min','max','quart_25', 'quart_75','median','ci_5', 'ci_95']
        df_return.reset_index(level=df_return.index.names, inplace=True)
        df_return = df_return.rename(columns={ 'index': 'stat'})

        file_out = exp.name + '_' + exp_type + '_' + group + '_desc.csv'
        pfiles.save_csv (df_return, file_out, exp.dir_out, False)
    return df_return

def get_desc_vsteps (exp, df_desc):
    factor = exp.factor1
    confidence = exp.conf
    exp_type = df_desc.iloc[0]['exp']
    
    df_desc = df_desc [['genotype','lateral', 'holding','vm-50','vm-40','vm-30','vm-20','vm-10','vm0','vm+10','vm+20','vm+30','vm+40','vm+50','vm+60']]
    if df_desc.empty:
        # Make new dataframe for stats
        df_return = pd.DataFrame(columns = ['stat','exp','genotype', 'lateral', 'vm', 'holding','vm-50','vm-40','vm-30','vm-20',
                        'vm-10','vm0','vm+10','vm+20','vm+30','vm+40','vm+50','vm+60'])
        df_return['stat'] = ['n', 'mean', 'sd', 'se', 'sum', 'min','max','quart_25', 'quart_75','median','ci_5', 'ci_95']
        df_return.index = ['n', 'mean', 'sd', 'se', 'sum', 'min','max','quart_25', 'quart_75','median','ci_5', 'ci_95']
        df_return.reset_index(level=df_return.index.names, inplace=True)
        df_return = df_return.rename(columns={ 'index': 'stat'})
        
    else: 
        group = df_desc.iloc[0][factor]
        # Calculations for all measures at once
        n = df_desc.groupby(factor).count()
        avg = df_desc.groupby(factor).mean()
        sd = df_desc.groupby(factor).std()
        se = df_desc.groupby(factor).sem()
        added = df_desc.groupby(factor).sum()
        minimum = df_desc.groupby(factor).min()
        maximum = df_desc.groupby(factor).max()
        quartile25 = df_desc.groupby(factor).quantile(q = 0.25, 
                axis = 0, numeric_only = True, interpolation = 'linear')
        quartile75 = df_desc.groupby(factor).quantile(q = 0.75, 
                axis = 0, numeric_only = True, interpolation = 'linear')
        median = df_desc.groupby(factor).median()
        conf = se * sp.stats.t._ppf((1+confidence)/2., n.astype('int')-1)
        conf_5 = avg - conf
        conf_95 = avg + conf
        n = n.reset_index()
        avg = avg.reset_index()
        sd = sd.reset_index()
        added = added.reset_index()
        minimum = minimum.reset_index()
        maximum = maximum.reset_index()
        quartile25 = quartile25.reset_index()
        quartile75 = quartile75.reset_index()
        median = median.reset_index()

        # Make new dataframe for stats
        df_return = pd.DataFrame()
        df_return = pd.concat([n,avg,sd,se,added,minimum,maximum,quartile25,quartile75,median,conf_5,conf_95], sort = False, ignore_index=True)
        df_return.index = ['n', 'mean', 'sd', 'se', 'sum', 'min','max','quart_25', 'quart_75','median','ci_5', 'ci_95']
        df_return.reset_index(level=df_return.index.names, inplace=True)
        df_return = df_return.rename(columns={ 'index': 'stat'})

        file_out = exp.name + '_' + exp_type + '_' + group + '_desc.csv'
        pfiles.save_csv (df_return, file_out, exp.dir_out, False)
    return df_return

def ttest_minis (exp, list_groupdfs):
    df_ttest = pd.DataFrame()
    for measure in exp.measures_minis:
        t, p = stats.ttest_ind(list_groupdfs[0][measure], list_groupdfs[1][measure])
        new_row = [{'exp': exp.name, 'measure': measure, 'group1': list_groupdfs[0].iloc[0]['lateral'], 'group2': list_groupdfs[1].iloc[0]['lateral'],
                    't_value': t, 'p_value': p}]
        df_ttest = df_ttest.append(new_row)
    return df_ttest

def ttest_isteps (exp, list_groupdfs):
    df_ttest = pd.DataFrame()
    t, p = stats.ttest_ind(list_groupdfs[0]['vm'], list_groupdfs[1]['vm'])
    new_row = [{'exp': exp.name, 'measure': 'vm', 'group1': list_groupdfs[0].iloc[0]['lateral'], 'group2': list_groupdfs[1].iloc[0]['lateral'],
                't_value': t, 'p_value': p}]
    df_ttest = df_ttest.append(new_row)
    return df_ttest

def ttest_vsteps (exp, list_groupdfs):
    df_ttest = pd.DataFrame()
    t, p = stats.ttest_ind(list_groupdfs[0]['holding'], list_groupdfs[1]['holding'])
    new_row = [{'exp': exp.name, 'measure': 'holding', 'group1': list_groupdfs[0].iloc[0]['lateral'], 'group2': list_groupdfs[1].iloc[0]['lateral'],
                't_value': t, 'p_value': p}]
    df_ttest = df_ttest.append(new_row)
    return df_ttest

def ttest_rinput (exp, list_groupdfs):
    df_ttest = pd.DataFrame()
    t, p = stats.ttest_ind(list_groupdfs[0]['r_input'], list_groupdfs[1]['r_input'])
    new_row = [{'exp': exp.name, 'measure': 'holding', 'group1': list_groupdfs[0].iloc[0]['lateral'], 'group2': list_groupdfs[1].iloc[0]['lateral'],
                't_value': t, 'p_value': p}]
    df_ttest = df_ttest.append(new_row)
    return df_ttest

def get_measures (list_desc, list_stats, exp, dir_out):
    # Get Desc by measure for all genotypes 
    # Routin for if one of the desc is empty
    i = 0
    list_measures = list()
    for measure in exp.measures:
        df_measure = pd.DataFrame()
        df_measure[measure] = list_stats # This will be the index column
        for desc in list_desc:                # one per group from the above desc stats section
            df_measure[exp.groups[i]] = desc[measure] 
            i = i + 1
        list_measures.append(df_measure)
        file_out = exp.name + '_coupling_' + measure + '_desc.csv'
        pfiles.save_csv (df_measure, file_out, dir_out, False)
        i = 0   
    return list_measures

def get_sorted (exp, dir_sort, df):
    dir_desc = exp.dir_in + dir_sort + exp.dir_group
    dir_meas = exp.dir_in + dir_sort + exp.dir_measure
    list_groupdfs = list() # creates list of group dfs for attempted
    list_groupdfs = pfiles.get_groups (exp, df) # Fills list of group dfs for attempted

    # Descriptive stats by group for attempted
    list_desc = list()
    list_meas = list()
    for groupdf in list_groupdfs:
        if groupdf.empty:
            pass
        else:
            new_desc = get_desc(exp, groupdf, dir_desc)
            list_desc.append(new_desc)
            list_meas = get_measures (list_desc, new_desc['stat'], exp, dir_meas)

def get_anova (exp, list_groups, measure):
    groups_anova = list()
    i = 0
    while i < len(list_groups):
        groups_anova.append(list_groups[i])
        i +=1

    list_anova = list()
    total_n = 0
    for df_group in groups_anova:
        new_measure = df_group[measure].dropna()
        list_anova.append(new_measure)
        total_n = total_n + len(new_measure)

    # One-way Anova with strain as the independent factor and 3 groups
    if len(groups_anova) == 3:
        f_value, p_anova = stats.f_oneway(list_anova[0], list_anova[1], list_anova[2])
    elif len(groups_anova) == 4:
        f_value, p_anova = stats.f_oneway(list_anova[0], list_anova[1], list_anova[2], list_anova[3])

    freedom1 = len(groups_anova)
    freedom2 = total_n - len(groups_anova)
    df_anova = pd.DataFrame({'measure': [measure], 'df_within': [freedom1], 'df_between': [freedom2], 
                            'f_value': [f_value], 'p_value': [p_anova]})
    anova_text =  'One-Way ANOVA: F(' + str(freedom1) + ',' + str(freedom2) + ') = ' + str(round(f_value,4)) + ', p = ' + str(round(p_anova, 4))
    df_anova['print'] = anova_text
    file_out = exp.name + '_coupling_' + measure + '_' + '_anova_' + '.csv'
    pfiles.save_csv (df_anova, file_out, exp.dir_anova, False)
    return df_anova
    

def get_tukey (exp, df_all, measure):
    # Tukey posthoc analysis
    # See https://jpktd.blogspot.com/2013/03/multiple-comparison-and-tukey-hsd-or_25.html
    # And https://code.google.com/archive/p/qsturng-py/
    # And https://stackoverflow.com/questions/48200699/how-can-i-get-p-values-of-each-group-comparison-when-applying-the-tukey-s-hones
    # q, res_table, std_pairs, etc can be found from print(dir(result)) which will list all possible calculations
    
    if len(df_all.groupby('strain').count()) >= 3:
        df_tukey = df_all[np.isfinite(df_all[measure])]
        mc = MultiComparison(df_tukey[measure], df_tukey['strain'])
        result = mc.tukeyhsd()
        p = psturng(np.abs(result.meandiffs/result.std_pairs), len(result.groupsunique), result.df_total) 
        df_pairs = pd.DataFrame({'group1': [result._results_table[1][0], result._results_table[2][0], result._results_table[3][0]],
                             'group2': [result._results_table[1][1], result._results_table[2][1], result._results_table[3][1]],
                             'p_value': [np.around(p[0], 4), np.around(p[1], 4), np.around(p[2], 4)]})
    else:
        df_pairs = pd.DataFrame({'group1': [], 'group2': [], 'p_value': []})
    
    file_out = exp.name + '_coupling_' + measure + '_' + '_tukey_' + '.csv'
    pfiles.save_csv (df_pairs, file_out, exp.dir_tukey, False)
    return df_pairs 

def get_hyp_coupling (exp, df_all, dir_in): 
    os.chdir(dir_in)
    list_groups = list() # creates list of group dfs for attempted
    list_groups = pfiles.get_groups (exp, df_all) # Fills list of group dfs for attempted
        # If strain = na then drop

    df_anova = pd.DataFrame()
    df_tukey = pd.DataFrame()
    for meas in exp.measures:
        df_new_anova = pd.DataFrame()
        df_new_tukey = pd.DataFrame()
        new_anova = get_anova(exp, list_groups, meas)
        new_anova['measure'] = meas
        new_anova.set_index ('measure', inplace = True)
        df_anova = df_anova.append(new_anova)
    
        new_tukey = get_tukey(exp, df_all, meas)
        new_tukey['measure'] = meas
        new_tukey.set_index ('measure', inplace = True)
        df_tukey = df_tukey.append(new_tukey)
    
    file_anova = exp.name + '_coupling_all_anovas.csv'
    file_tukey = exp.name + '_coupling_all_tukey.csv'
    pfiles.save_csv(df_anova, file_anova, dir_in + exp.dir_anova, True)
    pfiles.save_csv(df_tukey, file_tukey, dir_in + exp.dir_tukey, True)


def get_indexed_seps (exp, df, steps): # df contains all groups
    df_new = pd.DataFrame()
    i = 0 # keeps up with rows in df
    while i < len(df['exp']):
        for step in steps:
            if df.iloc[0]['exp'] == 'isteps':
                new_row = [{'exp':df.iloc[i]['exp'], 'file':df.iloc[i]['file'], 'genotype':df.iloc[i]['genotype'], 'slice':df.iloc[i]['slice'], 'shaved_side':df.iloc[i]['shaved_side'],
                        'hemisphere':df.iloc[i]['hemisphere'], 'lateral':df.iloc[i]['lateral'], 'vm':df.iloc[i]['vm'], 'step':step, 'value':df.iloc[i][step]}]
            elif df.iloc[0]['exp'] == 'vsteps':
                new_row = [{'exp':df.iloc[i]['exp'], 'file':df.iloc[i]['file'], 'genotype':df.iloc[i]['genotype'], 'slice':df.iloc[i]['slice'], 'shaved_side':df.iloc[i]['shaved_side'],
                        'hemisphere':df.iloc[i]['hemisphere'], 'lateral':df.iloc[i]['lateral'], 'holding':df.iloc[i]['holding'], 'step':step, 'value':df.iloc[i][step]}]
            df_new = df_new.append(new_row)
            #df_new = df_new.reindex()
        i += 1
    return df_new

 
        