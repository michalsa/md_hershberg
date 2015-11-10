from scipy import stats
import os, sqlite3, re
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pandas, matplotlib, sqlalchemy
from collections import OrderedDict
from itertools import cycle

from sklearn.cluster import AffinityPropagation
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs

OUTPUTDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/output'

THRESHOLDS = [0.1, 0.08, 0.05, 0.03, 0.02]

def plot_pangene_number_percentage(db_fname):
    con = sqlite3.connect(db_fname)
    matplotlib.style.use('ggplot')
    fig, axes = plt.subplots(2, len(THRESHOLDS))
    plt.subplots_adjust(hspace=0.2)
    with con:        
        fig.suptitle('Pangene Numbers by Origin\nand Percentages of Genes Annotated as Viral, by Fluidity Threshold', fontsize=22)
        c = ['r', 'g', 'b', 'y', 'c', 'm']
        for ix, threshold in enumerate(THRESHOLDS):
            threshold_str_file = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
            n = len(os.listdir(os.path.join(OUTPUTDIR, dir, 'fasta_out_pangenome', threshold_str_file, 'unified')))
            threshold_str_sql = '%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
            stats_table_name = 'pangenome_%s_stats' % threshold_str_sql
            stats_table_df = pandas.read_sql_query('SELECT pangene_id, n_genomes, viral_flag, origin FROM {} WHERE paralog_flag=0'.format(stats_table_name), 
                                                   con, index_col='pangene_id')
            stats_table_df['bin'] = pandas.cut(stats_table_df['n_genomes'], bins=[0, 1, round(0.25*n, 0), round(0.75*n, 0), n-1, n], 
                                               labels=['1_unique', '2_rare', '3_intermediate', '4_near-core', '5_core'])
            
            plot_data = stats_table_df.groupby(['bin', 'origin']).size()
            '''print plot_data.unstack().unstack()
            return'''
            #plot_data = stats_table_df.groupby(['n_genomes', 'viral_flag']).size()
            plot_data.unstack().plot(ax=axes[0, ix], kind='bar', stacked=True, ylim=(0,4000), colormap='spring', rot=325)#, rot=1)
            axes[0, ix].set_title('%s\n(n_genomes=%d, n_pangenes=%d)' % (threshold_str_file, n, plot_data.sum()), fontsize=12)
            '''patches, labels = axes[0,ix].get_legend_handles_labels()
            axes[0, ix].legend(patches,['non-viral', 'viral'], loc='upper left')'''
            pcnt_data = (stats_table_df[(stats_table_df['origin']=='viral')].groupby('bin').size())/(stats_table_df.groupby('bin').size())
            '''print pcnt_data
            return'''
            pcnt_data.plot(ax=axes[1, ix], kind='bar', stacked=True, ylim=(0,0.2), color=matplotlib.cm.spring(2), rot=325)#, rot=1)
            for k in range(2):
                #axes[k, ix].set_xticklabels(tuple([str(x+1) if (x==0 or (x+1)%5==0 or (x==n-1 and (x+1)%5>1)) else '' for x in range(n)]))
                axes[k, ix].set_xticklabels(['unique', 'rare', 'intermediate', 'near-core', 'core'])
                axes[k, ix].set_xlabel('')
        
        #axes[k, 2].set_xlabel('n_genomes')
        axes[1, 2].set_xlabel('bin')
        
        axes[0, 0].set_ylabel('n_pangenes')
        axes[1, 0].set_ylabel('frac_viral_pangenes')
        
        plt.show()
        return

def plot_enc_encp(db_fname):
    con = sqlite3.connect(db_fname)
    matplotlib.style.use('ggplot')
    fig, axes = plt.subplots(2, len(THRESHOLDS))
    plt.subplots_adjust(hspace=0.2)
    with con:
        p_val = 0.05
        for ix, threshold in enumerate(THRESHOLDS):
            threshold_str_file = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
            fig.suptitle('Average Viral and Non-Viral ENC and ENCprime Values by Fluidity Threshold', fontsize=22)
            n = len(os.listdir(os.path.join(OUTPUTDIR, dir, 'fasta_out_pangenome', threshold_str_file, 'unified')))
            threshold_str_sql = '%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
            stats_table_name = 'pangenome_%s_stats' % threshold_str_sql
            stats_table_df = pandas.read_sql_query('SELECT pangene_id, n_genomes, viral_flag, enc, encp, fop FROM {} WHERE paralog_flag=0'.format(stats_table_name), 
                                                   con, index_col='pangene_id')
            '''plot_data = stats_table_df.groupby(['n_genomes', 'viral_flag']).agg({'enc': np.average, 'encp': np.average})
            n_pangenes = stats_table_df.groupby(['n_genomes', 'viral_flag']).size().sum()'''
            stats_table_df['bin'] = pandas.cut(stats_table_df['n_genomes'], bins=[0, 1, round(0.25*n, 0), round(0.75*n, 0), n-1, n], 
                                               labels=['1_unique', '2_rare', '3_intermediate', '4_near-core', '5_core'])
            
            plot_data = stats_table_df.groupby(['bin', 'viral_flag']).agg({'enc': np.average, 'encp': np.average})
            n_pangenes = stats_table_df.groupby(['bin', 'viral_flag']).size().sum()
            plot_data['enc'].unstack().plot(ax=axes[0, ix], kind='bar', ylim=(20,61), colormap='summer', rot=325)#, rot=1)          
            
            axes[0, ix].set_title('%s\n(n_genomes=%d, n_pangenes=%d)' % (threshold_str_file, n, n_pangenes), fontsize=12)           
            plot_data['encp'].unstack().plot(ax=axes[1, ix], kind='bar', ylim=(20,61), colormap='summer', rot=325)#, rot=1)
            for k in range(2):
                #axes[k, ix].set_xticklabels(tuple([str(x+1) if (x==0 or (x+1)%5==0 or (x==n-1 and (x+1)%5>1)) else '' for x in range(n)]))
                axes[k, ix].set_xticklabels(['unique', 'rare', 'intermediate', 'near-core', 'core'])
                axes[k, ix].set_xlabel('')
                patches, labels = axes[k, ix].get_legend_handles_labels()
                axes[k, ix].legend(patches,['non-viral', 'viral'])
        
            for i, l in enumerate(['1_unique', '2_rare', '3_intermediate', '4_near-core', '5_core']):
                if True not in [str(m)=='nan' for m in plot_data['enc'].unstack().values[i-1]]:
                    if stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==l)]['enc'], 
                                          stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==l)]['enc'])[1]<p_val:
                        #print plot_data['enc'].unstack().values[i-1] 
                        axes[0, ix].annotate('*', (i-1, max(plot_data['enc'].unstack().values[i-1])*1.005))
                if True not in [str(m)=='nan' for m in plot_data['encp'].unstack().values[i-1]]:
                    if stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==l)]['encp'], 
                                          stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==l)]['encp'])[1]<p_val:
                        axes[1, ix].annotate('*', (i-1, max(plot_data['encp'].unstack().values[i-1])*1.005))
                        #print plot_data['encp'].unstack().values[i-1]
            '''
            for i in range(1, n+1):
                if True not in [str(m)=='nan' for m in plot_data['enc'].unstack().values[i-1]]:
                    if stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['n_genomes']==i)]['enc'], 
                                          stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['n_genomes']==i)]['enc'])[1]<p_val:
                        #print plot_data['enc'].unstack().values[i-1] 
                        axes[0, ix].annotate('*', (i-1, max(plot_data['enc'].unstack().values[i-1])*1.005))
                if True not in [str(m)=='nan' for m in plot_data['encp'].unstack().values[i-1]]:
                    if stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['n_genomes']==i)]['encp'], 
                                          stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['n_genomes']==i)]['encp'])[1]<p_val:
                        axes[1, ix].annotate('*', (i-1, max(plot_data['encp'].unstack().values[i-1])*1.005))
                        #print plot_data['encp'].unstack().values[i-1]
            '''
            for k in range(2):
                axes[k, 2].set_xlabel('bins')
                #axes[k, 2].set_xlabel('n_genomes')
            
            axes[0, 0].set_ylabel('enc_avg')
            axes[1, 0].set_ylabel('encp_avg')
            
            fig.text(0.05, 0.01, '* signifies p_value of less than 0.05 in Mann-Whitney U test', verticalalignment='bottom')
            
        plt.show()
        return

def plot_fop(db_fname):
    con = sqlite3.connect(db_fname)
    with con:
        matplotlib.style.use('ggplot')
        p_val = 0.05
        fig, axes = plt.subplots(5,1)
        plt.subplots_adjust(hspace=0.5)
        fig.suptitle('Average Viral and Non-Viral FOP Values by Fluidity Threshold', fontsize=22)
        for ix, threshold in enumerate(THRESHOLDS):
            threshold_str_file = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
            n = len(os.listdir(os.path.join(OUTPUTDIR, dir, 'fasta_out_pangenome', threshold_str_file, 'unified')))
            threshold_str_sql = '%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
            stats_table_name = 'pangenome_%s_stats' % threshold_str_sql
            stats_table_df = pandas.read_sql_query('SELECT pangene_id, n_genomes, viral_flag, enc, encp, fop FROM {} WHERE paralog_flag=0'.format(stats_table_name), 
                                                   con, index_col='pangene_id')
            '''plot_data = stats_table_df.groupby(['n_genomes', 'viral_flag']).agg({'fop': np.average})
            n_pangenes = stats_table_df.groupby(['n_genomes', 'viral_flag']).size().sum()'''
            
            stats_table_df['bin'] = pandas.cut(stats_table_df['n_genomes'], bins=[0, 1, round(0.25*n, 0), round(0.75*n, 0), n-1, n], 
                                               labels=['1_unique', '2_rare', '3_intermediate', '4_near-core', '5_core'])
            plot_data = stats_table_df.groupby(['bin', 'viral_flag']).agg({'fop': np.average})
            n_pangenes = stats_table_df.groupby(['bin', 'viral_flag']).size().sum()
            plot_data['fop'].unstack().plot(ax=axes[ix], kind='bar', ylim=(0,0.6), colormap='spring', rot=1)          
            axes[ix].set_title('%s\n(n_genomes=%d, n_pangenes=%d)' % (threshold_str_file, n, n_pangenes), fontsize=12)
            #axes[ix].set_xticklabels(tuple([str(x+1) if (x==0 or (x+1)%5==0 or (x==n-1 and (x+1)%5>1)) else '' for x in range(n)]))
            axes[ix].set_xticklabels(['unique', 'rare', 'intermediate', 'near-core', 'core'])
            axes[ix].set_xlabel('')
            patches, labels = axes[ix].get_legend_handles_labels()
            axes[ix].legend(patches,['non-viral', 'viral'], loc='lower right')
        
            '''for i in range(1, n+1):
                if True not in [str(m)=='nan' for m in plot_data['fop'].unstack().values[i-1]]:
                    if stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['n_genomes']==i)]['fop'], 
                                          stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['n_genomes']==i)]['fop'])[1]<p_val:
                        #print plot_data['enc'].unstack().values[i-1] 
                        axes[ix].annotate('*', (i-1, max(plot_data['fop'].unstack().values[i-1])*1.005))'''
            
            for i, l in enumerate(['1_unique', '2_rare', '3_intermediate', '4_near-core', '5_core']):
                if True not in [str(m)=='nan' for m in plot_data['fop'].unstack().values[i-1]]:
                    if stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==l)]['fop'], 
                                          stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==l)]['fop'])[1]<p_val:
                        #print plot_data['enc'].unstack().values[i-1] 
                        axes[ix].annotate('*', (i-1, max(plot_data['fop'].unstack().values[i-1])*1.005))
            
        #axes[4].set_xlabel('n_genomes')
        axes[4].set_xlabel('bins')
        axes[2].set_ylabel('fop_avg')
        
        fig.text(0.05, 0.05, '* signifies p_value of less than 0.05 in Mann-Whitney U test', verticalalignment='bottom')
        
        plt.show()
        return

def plot_viral_pcnt(db_fname):
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    matplotlib.style.use('ggplot')
    fig, axes = plt.subplots(2, 3)
    plt.subplots_adjust(hspace=0.2)
    fig.suptitle('Percentages of Pangenes Annotated as Viral,\nFlagged as Viral vs. Non-Viral by PhiSpy, by Fluidity Threshold', fontsize=22)
    c = ['r', 'g', 'b', 'y', 'c', 'm']
    for ix, threshold in enumerate(THRESHOLDS):
        threshold_str_file = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        n = len(os.listdir(os.path.join(OUTPUTDIR, dir, 'fasta_out_pangenome', threshold_str_file, 'unified')))
        threshold_str_sql = '%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        stats_table_name = 'pangenome_%s_stats' % threshold_str_sql
        stats_table_df = pandas.read_sql_query('SELECT pangene_id, n_genomes, viral_flag, origin FROM {} WHERE paralog_flag=0'.format(stats_table_name), 
                                               engine, index_col='pangene_id')
        
        stats_table_df['bin'] = pandas.cut(stats_table_df['n_genomes'], bins=[0, 1, round(0.25*n, 0), round(0.75*n, 0), n-1, n], 
                                           labels=['1_unique', '2_rare', '3_intermediate', '4_near-core', '5_core'])
        
        plot_data = stats_table_df.groupby(['bin', 'viral_flag', 'origin']).size()
        n_pangenes = plot_data.sum()
        
        print plot_data
        #return
        plot_df = (plot_data.unstack()['viral']/plot_data.unstack(level=[0, 1]).sum()).unstack()
        print plot_df
        
        plot_df.plot(ax=axes[ix//3, ix%3], kind='bar', ylim=(0,1.05), color=[c[ix], c[-1]], rot=325)
        
        axes[ix//3, ix%3].set_title('%s\n(n_genomes=%d, n_pangenes=%d)' % (threshold_str_file, n, n_pangenes), fontsize=12)
        
        patches, labels = axes[ix//3, ix%3].get_legend_handles_labels()
        axes[ix//3, ix%3].legend(patches,['non-viral', 'viral'], loc='lower right')
        
        axes[ix//3, ix%3].set_xticklabels(['unique', 'rare', 'intermediate', 'near-core', 'core'])
        axes[ix//3, ix%3].set_xlabel('')
    
    axes[1,1].set_xlabel('bins')
    axes[-1, -1].axis('off')
            
    plt.show()
    return

def plot_origin_piechart(db_fname):
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    matplotlib.style.use('ggplot')
    fig, axes = plt.subplots(2, 3)
    plt.subplots_adjust(hspace=0.2)
    fig.suptitle('Annotations of Pangenes Flagged as Non-Viral by Fluidity Threshold', fontsize=22)
    for ix, threshold in enumerate(THRESHOLDS):
        threshold_str_file = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        n = len(os.listdir(os.path.join(OUTPUTDIR, dir, 'fasta_out_pangenome', threshold_str_file, 'unified')))
        threshold_str_sql = '%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        stats_table_name = 'pangenome_%s_stats' % threshold_str_sql
        stats_table_df = pandas.read_sql_query('SELECT pangene_id, viral_flag, origin FROM {} WHERE paralog_flag=0'.format(stats_table_name), 
                                               engine, index_col='pangene_id')
        plot_data = stats_table_df.groupby(['viral_flag', 'origin']).size()
        n_pangenes = plot_data.sum()
        
        print plot_data[0]
        #return
        
        plot_data[0].plot(ax=axes[ix//3, ix%3], kind='pie', autopct='%.2f', fontsize=16, labels=['', '', ''])
        if ix//3==1 and ix%3==1:
            axes[ix//3, ix%3].legend(bbox_to_anchor=(0.85,0.30), fontsize=18,
                                     bbox_transform=plt.gcf().transFigure, labels=plot_data[1].index)
        
        axes[ix//3, ix%3].set_title('%s\n(n_genomes=%d, n_pangenes=%d)' % (threshold_str_file, n, n_pangenes), fontsize=18)
        
        axes[ix//3, ix%3].set_xlabel('')
    '''
    patches, labels = axes[1,1].get_legend_handles_labels()
    fig.legend(patches, ('bacterial', 'unknown', 'viral'), bbox_to_anchor=(0.85,0.30), fontsize=18,
           bbox_transform=plt.gcf().transFigure)          
    '''
    axes[-1, -1].axis('off')
    plt.show()
    return

def plot_annotation_piechart(db_fname):
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    matplotlib.style.use('ggplot')
    #fig, axes = plt.subplots(2,3)
    fig, axes = plt.subplots(2,5)
    #plt.subplots_adjust(hspace=0.4, wspace=0.4)
    #plt.subplots_adjust(hspace=0.4)
    p_val = 0.1
    fig.suptitle('10 Most Common Annotations for Core Pangenes Flagged as Non-Viral vs. Viral by Fluidity Threshold', fontsize=22)
    for ix, threshold in enumerate(THRESHOLDS):
        threshold_str_file = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        n = len(os.listdir(os.path.join(OUTPUTDIR, dir, 'fasta_out_pangenome', threshold_str_file, 'unified')))
        threshold_str_sql = '%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        stats_table_name = 'pangenome_%s_stats' % threshold_str_sql
        stats_table_df = pandas.read_sql_query('SELECT pangene_id, n_genomes, viral_flag, annotation FROM {} WHERE paralog_flag=0'.format(stats_table_name), 
                                               engine, index_col='pangene_id')
        
        cats_dict = OrderedDict([('ribosomal', ['ribosom', 'rrna']),
                                 ('phage', ['phage', 'virion', 'tail', 'capsid', 'head', '\Wq protein', 
                                            'anti[ -]?terminat', '\Wearly\W', '\Wlate\W', 'anti[ -]?repress', 'holin', 'packag',
                                            'endolysin', '\Wlysis\W', '\Wlytic\W', '\Wvlp\W', 'polarity suppress', 'dna inject',
                                            'lambda', '\Wcro\W', '\Wnin\W']),
                                 ('virulence', ['invasi', 'adhesi', 'secret', 'fim?b', 'pil', 'ha?emolys', 
                                                'necro', 'tox', 't\dss', 'effector inject', 'opacity', 'effect', 
                                                'competence', 'virulen', 'pathogen', 'infect', 'intimin', 'usher',
                                                'von willebrand factor', 'vwf', 'entericidin', 'mch', 'curli', 'pap',
                                                'rata', 'injection protein', 'pho']),
                                 ('antibiotic', ['antibiotic', '\w+cin', 'drug', 'penicillin']),
                                 ('trna', ['trna']),
                                 ('common metab.', ['carb', 'lip', 'amin', 'saccaride', 'fatty', 'nitr',
                                                    'sugar', 'glyco', '\w+ose', '\w+ol', 'pts', 'lactam', 
                                                    'galactos']),
                                 ('rare metab.', ['nickel', 'chloride', 'cobalt', 'copper', 'metal', 'zinc', 
                                                  'hem', 'ferr', 'iron', 'mercury', 'tellurite', 'molybdenum',
                                                  'silver']),
                                 ('cell cycle', ['replicat', 'translat', 'divi', 'transcript', 'polymerase', 'growth', 
                                                 'chromosom', 'primosom', 'elongat', 'dna \w+ process']),
                                 ('defense', ['restriction', 'modification', 'crispr', 'cas', 'immunity', 'resist',
                                              'radical', 'dedd', 'tery', 'tolera']),
                                 ('transport', ['port', 'transfer', 'flux', 'channel', 'pump', 'carr', 
                                                'porin', 'transloc', 'mobiliz', 'facilitat', 'relay', 'siderophore']),
                                 ('dna related', ['dna', 'nucleo', 'strand[ -]bind', 'helix[- ]?destab', 'histone', 'mismatch',
                                                  'anneal', 'rec']),
                                 ('expression', ['modulat', 'initiat', 'regulator', 'repress', 'leader', 'rho',
                                                 'sigma', 'inhibit', 'modif', 'dead box']),
                                 ('protein assembly', ['chaperon', 'assembly']),
                                 #('MGEs', ['transpos', 'insertion', 'plasmid', 'recombin', 'addiction', 'dna uptake', 
                                           #'entry exclu', 'rhs', 'is[\d| |$]']),
                                 ('structure', ['membrane', 'flagel', 'structur', 'capsul', 'shape', 'sept',
                                                'scaffold', 'microcompartment', 'peptidoglycan']),
                                 ('lifestyle', ['biofilm', 'filamentation', 'specificity', 'host', 'motil', 'nodul',
                                                'spor']),
                                 #('pseudogene', ['pseudo']),
                                 ('stress response', ['stress', 'sens', 'shock', 'starv']),
                                 ('communication', ['signal', 'transduc', 'pix']),
                                 ('enzymes', ['synthe', '\w+ase', 'metabolism', '[ag]tp', 'biogene', 'enzym',
                                              'interchange', 'doxin', 'degrad']),
                                 ('unknown', ['hypothetical', 'uncharacterized', 'upf\d', 'duf\d', 'unknown', 'ipr\d', 
                                              'putative'])])
                
        stats_table_df['bin'] = pandas.cut(stats_table_df['n_genomes'], bins=[0, 1, round(0.25*n, 0), round(0.75*n, 0), n-1, n], 
                                           labels=['1_unique', '2_rare', '3_intermediate', '4_near-core', '5_core'])
        
        for stats_row in stats_table_df.iterrows():
            ann = stats_row[1]['annotation']
            ann_exists = False
            for cat_key, cat_vals in cats_dict.iteritems():
                if [re.search(cat_val, ann) for cat_val in cat_vals]!=[None]*len(cat_vals):
                    stats_table_df.loc[stats_row[0], 'category'] = cat_key
                    ann_exists = True
                    break
            if not ann_exists:
                stats_table_df.loc[stats_row[0], 'category'] = 'other'
            #print stats_table_df.loc[stats_row[0]]
        plot_data_flag = stats_table_df.groupby(['viral_flag', 'bin', 'category']).size()
        
        plot_data = pandas.DataFrame([plot_data_flag[i].unstack().transpose()['5_core']/plot_data_flag[i].unstack().transpose()['5_core'].sum()
                                      for i in range(2)]).transpose()
        plot_data.columns = ['non viral', 'viral']
        print plot_data
        
        plot_data['non viral'].order(ascending=False)[:10].plot(ax=axes[0, ix], kind='pie', colormap='rainbow')
        plot_data['viral'].order(ascending=False)[:10].plot(ax=axes[1, ix], kind='pie', colormap='rainbow')
        axes[0, ix].set_title('%s\n(n_genomes=%d,\nn_pangenes=%d)' % (threshold_str_file, n, 
                                                                     sum([plot_data_flag[i].unstack().transpose()['5_core'].sum() for i in range(2)])))#, fontsize=18)
        
        patches_nonvir, labels_nonvir = axes[0, ix].get_legend_handles_labels()
        patches_vir, labels_vir = axes[1, ix].get_legend_handles_labels()
        for patch_vir, label_vir in zip(patches_vir, labels_vir):
            if label_vir not in labels_nonvir:
                print label_vir, (np.cos(np.radians((patch_vir.theta2-patch_vir.theta1)/2+patch_vir.theta1))/2,
                                  np.sin(np.radians((patch_vir.theta2-patch_vir.theta1)/2+patch_vir.theta1))/2)
                #return
                axes[1,ix].annotate('*', xy=(np.cos(np.radians((patch_vir.theta2-patch_vir.theta1)/2+patch_vir.theta1))/2,
                                             np.sin(np.radians((patch_vir.theta2-patch_vir.theta1)/2+patch_vir.theta1))/2))
        for patch_nonvir, label_nonvir in zip(patches_nonvir, labels_nonvir):
            if label_nonvir not in labels_vir:
                print label_nonvir, (np.cos(np.radians((patch_nonvir.theta2-patch_nonvir.theta1)/2+patch_nonvir.theta1))/2,
                                     np.sin(np.radians((patch_nonvir.theta2-patch_nonvir.theta1)/2+patch_nonvir.theta1))/2)
                #return
                axes[0,ix].annotate('*', xy=(np.cos(np.radians((patch_nonvir.theta2-patch_nonvir.theta1)/2+patch_nonvir.theta1))/2,
                                             np.sin(np.radians((patch_nonvir.theta2-patch_nonvir.theta1)/2+patch_nonvir.theta1))/2))
            
        #for i in range(2): axes[i, ix].axis('equal')
        
        '''
        chsq_calc_data = plot_data_flag.unstack(level=1)['5_core'].unstack()
        #print chsq_calc_data
        
        for v, val in enumerate(plot_data.index.values):
            if True not in [str(m)=='nan' for m in plot_data.transpose()[val].values] and False not in [m>5 for m in chsq_calc_data[val].values]:
                chsq_obs = [list(chsq_calc_data[val].values), 
                            [chsq_calc_data.transpose().sum()[i] - chsq_calc_data[val][i] for i in range(2)]]
                #print chsq_obs
                chsq_exp = [list((chsq_calc_data.transpose().sum()*chsq_calc_data[val].sum()/chsq_calc_data.transpose().sum().sum()).values), 
                            [chsq_calc_data.transpose().sum()[i] - chsq_calc_data.transpose().sum()[i]*chsq_calc_data[val].sum()/chsq_calc_data.transpose().sum().sum() for i in range(2)]]
                
                #print chsq_exp
                #print stats.chisquare(chsq_obs, f_exp=chsq_exp, ddof=1, axis=None)
                #return
                if stats.chisquare(chsq_obs, f_exp=chsq_exp, ddof=1, axis=None)[1]<p_val:
                    print val, plot_data.transpose()[val].values
                    axes[ix//3, ix%3].annotate('*', xy=(v, max(plot_data.transpose()[val].values)*1.005))
                    
        axes[ix//3, ix%3].set_xlabel('')'''
    
    #axes[-1, -1].axis('off')
    
    plt.show()
    return

def plot_codon_bias_indexes_by_origin_per_bin(db_fname):
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    matplotlib.style.use('ggplot')
    p_val = 0.05
    fig, axes = plt.subplots(2,5)
    plt.subplots_adjust(hspace=0.1)
    
    for ix, threshold in enumerate(THRESHOLDS):
        threshold_str_file = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        n = len(os.listdir(os.path.join(OUTPUTDIR, dir, 'fasta_out_pangenome', threshold_str_file, 'unified')))
        
        fig.suptitle('Average Codon Bias Index Values per Pangenomic Bin,\nAnnotated as Bacterial vs. Annotated as Viral', fontsize=22)
        
        threshold_str_sql = '%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        stats_table_name = 'pangenome_%s_stats' % threshold_str_sql
        stats_table_df = pandas.read_sql_query('SELECT pangene_id, n_genomes, viral_flag, encp, fop, origin FROM {} WHERE paralog_flag=0'.format(stats_table_name), 
                                               engine, index_col='pangene_id')
        
        stats_table_df['bin'] = pandas.cut(stats_table_df['n_genomes'], bins=[0, 1, round(0.25*n, 0), round(0.75*n, 0), n-1, n], 
                                           labels=['1_unique', '2_rare', '3_intermediate', '4_near-core', '5_core'])
                    
        plot_data = stats_table_df.groupby(['bin', 'origin']).agg({'fop': np.average, 'encp': np.average})
        tot_avg = stats_table_df.groupby(['bin']).agg({'fop': np.average, 'encp': np.average})
        
        encp_unstacked = plot_data['encp'].unstack()
        encp_unstacked['total'] = tot_avg['encp']
        
        fop_unstacked = plot_data['fop'].unstack()
        fop_unstacked['total'] = tot_avg['fop']
        '''
        plot_encp = pandas.DataFrame(data={'bacterial_nonviral': encp_unstacked.unstack().loc[:,'bacterial'][0],
                                           'phage_viral': encp_unstacked.unstack().loc[:,'viral'][1]},
                                    index=encp_unstacked.unstack().index)
        
        plot_fop = pandas.DataFrame(data={'bacterial_nonviral': fop_unstacked.unstack().loc[:,'bacterial'][0],
                                          'phage_viral': fop_unstacked.unstack().loc[:,'viral'][1]},
                                    index=fop_unstacked.unstack().index)
        
        print plot_encp
        print plot_fop
        #return
        
        plot_encp.plot(ax=axes[0,ix], kind='bar', ylim=(40,60), colormap='spring', rot=15, fontsize=10)
        plot_fop.plot(ax=axes[1,ix], kind='bar', ylim=(0,0.6), colormap='winter', rot=15, fontsize=10)
        '''
        encp_unstacked.loc[:,['bacterial','viral']].plot(ax=axes[0,ix], kind='bar', ylim=(40,60), colormap='spring', rot=15, fontsize=10)
        fop_unstacked.loc[:,['bacterial','viral']].plot(ax=axes[1,ix], kind='bar', ylim=(0,0.6), colormap='winter', rot=15, fontsize=10)
        
        '''tot_avg['encp'].plot(ax=axes[0,l], kind='bar', colormap='spring', rot=0, fontsize=12)
        tot_avg['fop'].plot(ax=axes[1,l], kind='bar', colormap='winter', rot=0, fontsize=12)'''
        
        axes[0, ix].set_title('Fluidity Threshold %s\n(n_genomes=%d, n_pangenes=%d)' % (threshold_str_file, n, stats_table_df.index.size), fontsize=14)
        for k in range(2): axes[k, ix].set_xlabel('')
        '''for k in range(2):
            patches, labels = axes[k, ix].get_legend_handles_labels()
            axes[k, ix].legend(patches,['non-viral', 'viral'], loc='upper right')
            axes[k, ix].set_xlabel('')'''
        '''
        #significance per origin
        for i, o in enumerate(['bacterial', 'unknown', 'viral']):
            if True not in [str(m)=='nan' for m in encp_unstacked[o].values]:
                if stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==label) & (stats_table_df['origin']==o)]['encp'], 
                                      stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==label) & (stats_table_df['origin']==o)]['encp'])[1]<p_val:
                    print stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==label) & (stats_table_df['origin']==o)]['encp'], 
                                             stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==label) & (stats_table_df['origin']==o)]['encp'])[1]
                    axes[0, l].annotate('*', (i, max(encp_unstacked[o].values)*1.005))
            
            if True not in [str(m)=='nan' for m in fop_unstacked[o].values]:
                if stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==label) & (stats_table_df['origin']==o)]['fop'], 
                                      stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==label) & (stats_table_df['origin']==o)]['fop'])[1]<p_val:
                    print stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==label) & (stats_table_df['origin']==o)]['fop'], 
                                             stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==label) & (stats_table_df['origin']==o)]['fop'])[1]
                    axes[1, l].annotate('*', (i, max(fop_unstacked[o].values)*1.005))
        '''
        '''significance of average value'''
        '''if True not in [str(m)=='nan' for m in encp_unstacked['total'].values]:
            if stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==label)]['encp'], 
                                  stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==label)]['encp'])[1]<p_val:
                print encp_unstacked['total'].values
                print stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==label)]['encp'], 
                                         stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==label)]['encp'])[1]
                axes[0, l].annotate('*', (0.5, max(encp_unstacked['total'].values)*1.005))
        
        if True not in [str(m)=='nan' for m in fop_unstacked['total'].values]:
            if stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==label)]['fop'], 
                                  stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==label)]['fop'])[1]<p_val:
                print fop_unstacked['total'].values
                print stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==label)]['fop'], 
                                         stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==label)]['fop'])[1]
                axes[1, l].annotate('*', (0.5, max(fop_unstacked['total'].values)*1.005)) '''
        '''
        #significance for phispy_flag=1, bacterial vs viral
        for l, label in enumerate(plot_encp.index):
            if True not in [str(m)=='nan' for m in plot_encp.loc[label].values]:
                if stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==label) & (stats_table_df['origin']=='bacterial')]['encp'], 
                                      stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==label) & (stats_table_df['origin']=='viral')]['encp'])[1]<p_val:
                    print stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==label) & (stats_table_df['origin']=='bacterial')]['encp'], 
                                             stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==label) & (stats_table_df['origin']=='viral')]['encp'])[1]
                    axes[0, ix].annotate('*', (l, max(plot_encp.loc[label].values)*1.005))
        
        for l, label in enumerate(plot_encp.index):
            if True not in [str(m)=='nan' for m in plot_fop.loc[label].values]:
                if stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==label) & (stats_table_df['origin']=='bacterial')]['fop'], 
                                      stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==label) & (stats_table_df['origin']=='viral')]['fop'])[1]<p_val:
                    print stats.mannwhitneyu(stats_table_df[(stats_table_df['viral_flag']==0) & (stats_table_df['bin']==label) & (stats_table_df['origin']=='bacterial')]['fop'], 
                                             stats_table_df[(stats_table_df['viral_flag']==1) & (stats_table_df['bin']==label) & (stats_table_df['origin']=='viral')]['fop'])[1]
                    axes[1, ix].annotate('*', (l, max(plot_fop.loc[label].values)*1.005))
    '''
        
        #significance for bacterial vs viral, regardless of phispy results
        for l, label in enumerate(encp_unstacked.index):
            if True not in [str(m)=='nan' for m in encp_unstacked.loc[label].values]:
                if stats.mannwhitneyu(stats_table_df[(stats_table_df['bin']==label) & (stats_table_df['origin']=='bacterial')]['encp'], 
                                      stats_table_df[(stats_table_df['bin']==label) & (stats_table_df['origin']=='viral')]['encp'])[1]<p_val:
                    print stats.mannwhitneyu(stats_table_df[(stats_table_df['bin']==label) & (stats_table_df['origin']=='bacterial')]['encp'], 
                                             stats_table_df[(stats_table_df['bin']==label) & (stats_table_df['origin']=='viral')]['encp'])[1]
                    axes[0, ix].annotate('*', (l, max(encp_unstacked.loc[label].values)*1.005))
        
        for l, label in enumerate(fop_unstacked.index):
            if True not in [str(m)=='nan' for m in fop_unstacked.loc[label].values]:
                if stats.mannwhitneyu(stats_table_df[(stats_table_df['bin']==label) & (stats_table_df['origin']=='bacterial')]['fop'], 
                                      stats_table_df[(stats_table_df['bin']==label) & (stats_table_df['origin']=='viral')]['fop'])[1]<p_val:
                    print stats.mannwhitneyu(stats_table_df[(stats_table_df['bin']==label) & (stats_table_df['origin']=='bacterial')]['fop'], 
                                             stats_table_df[(stats_table_df['bin']==label) & (stats_table_df['origin']=='viral')]['fop'])[1]
                    axes[1, ix].annotate('*', (l, max(fop_unstacked.loc[label].values)*1.005))
        
        fig.text(0.5, 0.05, 'bin', verticalalignment='bottom', fontsize=14)
        
        axes[0,0].set_ylabel('encp_avg')
        axes[1,0].set_ylabel('fop_avg')
        
        fig.text(0.04, 0.04, '* signifies p_value of less than %.2f in Mann-Whitney U test' % p_val, verticalalignment='bottom', fontsize=14)
        
    plt.show()
    return

def plot_codon_bias(db_fname):
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    matplotlib.style.use('ggplot')
    
    for ix, threshold in enumerate(THRESHOLDS):
        threshold_str_file = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        print threshold_str_file
                
        n = len(os.listdir(os.path.join(OUTPUTDIR, dir, 'fasta_out_pangenome', threshold_str_file, 'unified')))
        threshold_str_sql = '%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        stats_table_name = 'pangenome_%s_stats' % threshold_str_sql
        stats_table_df = pandas.read_sql_query('SELECT pangene_id, n_genomes, viral_flag, encp, fop, origin FROM {} WHERE paralog_flag=0'.format(stats_table_name), 
                                               engine, index_col='pangene_id')
        
        fig, axes = plt.subplots(2,5)
        plt.subplots_adjust(hspace=0.1)
        fig.suptitle('Distribution of Codon Bias Index Values per Pangenomic Bin for Fluidity Threshold %s\n(%d genomes, %d pangenes)' % 
                     (threshold_str_file, n, stats_table_df.index.size), fontsize=22)
        
        stats_table_df['bin'] = pandas.cut(stats_table_df['n_genomes'], bins=[0, 1, round(0.25*n, 0), round(0.75*n, 0), n-1, n], 
                                           #labels=['1_unique', '2_rare', '3_intermediate', '4_near-core', '5_core'])
                                           labels=[1, 2, 3, 4, 5])
        
        stats_table_df['origin_bin'] = stats_table_df['origin'].apply(lambda x: ['bacterial', 'viral', 'unknown'].index(x))
                
        for l, label in enumerate(['1_unique', '2_rare', '3_intermediate', '4_near-core', '5_core']):
            #plot_data = stats_table_df[stats_table_df['bin']=='1_unique']
            plot_data = stats_table_df[stats_table_df['bin']==l+1]
            
            plot_data.plot(ax=axes[0,l], kind='scatter', x='fop', y='encp', c='viral_flag', cmap='spring', edgecolor='k', colorbar=False, fontsize=12)        
            #plot_data[plot_data['viral_flag']==1].plot(ax=axes[ix], kind='scatter', x='fop', y='encp', color='DarkGreen', label='viral', fontsize=12)
                    
            axes[0,l].set_title('%s\n%d pangenes' % (label, plot_data.index.size), fontsize=14)
            axes[0,l].annotate('non viral: $\overline{FOP}=%.3f; \overline{ENC\'}=%.2f$\nviral: $\overline{FOP}=%.3f; \overline{ENC\'}=%.2f$' % 
                               (np.average(plot_data[plot_data['viral_flag']==0]['fop']),
                                np.average(plot_data[plot_data['viral_flag']==0]['encp']),
                                np.average(plot_data[plot_data['viral_flag']==1]['fop']),
                                np.average(plot_data[plot_data['viral_flag']==1]['encp'])),
                               (0.02, 0.04), size=10, xycoords="axes fraction", va="center", ha="left")
            
            plot_data.plot(ax=axes[1,l], kind='scatter', x='fop', y='encp', c='origin_bin', cmap='rainbow', edgecolor='k', colorbar=False, fontsize=12)
            axes[1,l].annotate('bacterial: $\overline{FOP}=%.3f; \overline{ENC\'}=%.2f$\nphage: $\overline{FOP}=%.3f; \overline{ENC\'}=%.2f$\nunknown: $\overline{FOP}=%.3f; \overline{ENC\'}=%.2f$' % 
                               (np.average(plot_data[plot_data['origin_bin']==0]['fop']),
                                np.average(plot_data[plot_data['origin_bin']==0]['encp']),
                                np.average(plot_data[plot_data['origin_bin']==1]['fop']),
                                np.average(plot_data[plot_data['origin_bin']==1]['encp']),
                                np.average(plot_data[plot_data['origin_bin']==2]['fop']),
                                np.average(plot_data[plot_data['origin_bin']==2]['encp'])),
                               (0.02, 0.06), size=10, xycoords="axes fraction", va="center", ha="left")
            
            '''print label
            print 'non viral: average fop: %.3f; average encp: %.2f' % (np.average(plot_data[plot_data['viral_flag']==0]['fop']), 
                                                                        np.average(plot_data[plot_data['viral_flag']==0]['encp']))
            print 'viral: average fop: %.3f; average encp: %.2f' % (np.average(plot_data[plot_data['viral_flag']==1]['fop']), 
                                                                    np.average(plot_data[plot_data['viral_flag']==1]['encp']))'''
        
        nonvirArtist = plt.Line2D((0,0),(0,0), color=cm.spring(0), marker='o', linestyle='')
        virArtist = plt.Line2D((0,0),(0,0), color=cm.spring(256), marker='o', linestyle='')
        fig.legend([nonvirArtist, virArtist], ['non viral', 'viral'], bbox_to_anchor=(0.1, 0.74), title='PhiSpy Category')
        
        bacArtist = plt.Line2D((0,0),(0,0), color=cm.rainbow(0), marker='o', linestyle='')
        phaArtist = plt.Line2D((0,0),(0,0), color=cm.rainbow(256//2), marker='o', linestyle='')
        unkArtist = plt.Line2D((0,0),(0,0), color=cm.rainbow(256), marker='o', linestyle='')
        fig.legend([bacArtist, phaArtist, unkArtist], ['bacteria', 'phage', 'unknown'], bbox_to_anchor=(0.1, 0.34), title='Origin by Annotation')
        
        plt.show()
    return
        
            
if __name__ == '__main__':
    for dir in os.listdir(OUTPUTDIR):
        if os.path.isdir(os.path.join(OUTPUTDIR, dir)):
            print dir
            plot_pangene_number_percentage(os.path.join(OUTPUTDIR, dir, 'output_datatables.db'))
            #plot_enc_encp(os.path.join(OUTPUTDIR, dir, 'output_datatables.db'))
            #plot_fop(os.path.join(OUTPUTDIR, dir, 'output_datatables.db'))
            #plot_viral_pcnt(os.path.join(OUTPUTDIR, dir, 'output_datatables.db'))
            #plot_origin_piechart(os.path.join(OUTPUTDIR, dir, 'output_datatables.db'))
            #plot_codon_bias_indexes_by_origin_per_bin(os.path.join(OUTPUTDIR, dir, 'output_datatables.db'))
            #plot_annotation_piechart(os.path.join(OUTPUTDIR, dir, 'output_datatables.db'))
            #plot_codon_bias(os.path.join(OUTPUTDIR, dir, 'output_datatables.db'))
            
