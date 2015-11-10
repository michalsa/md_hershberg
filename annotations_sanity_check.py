import pandas
import sqlite3
import numpy as np
import sqlalchemy
import os, re
from Bio import SeqIO

OUTPUTDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/output'
GBKDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/all_gbk'

THRESHOLDS = [0.1, 0.08, 0.05, 0.03, 0.02]

VIR_KEYWORDS = ['viral', 'virus', 'phage', 'capsid', 'tail']

def acquired_genes_by_annotations(db_fname, org_name):
    '''check how many pangenes flagged as viral vs. non-viral have viral keywords in their annotations'''
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    strains_df = pandas.read_sql_table('strains_plasmids', engine, index_col='strain_id')
    
    for threshold in THRESHOLDS:
        threshold_str_file = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        n = len(os.listdir(os.path.join(OUTPUTDIR, dir, 'fasta_out_pangenome', threshold_str_file, 'unified')))
        threshold_str_sql = '%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        pangenome_table_name, stats_table_name = 'pangenome_%s' % threshold_str_sql, 'pangenome_%s_stats' % threshold_str_sql
        pangenome_table_df = pandas.read_sql_table(pangenome_table_name, engine, index_col='pangene_id')          
        stats_table_df = pandas.read_sql_table(stats_table_name, engine, index_col='pangene_id')
        
        acquired_pangenes_stats = stats_table_df[(stats_table_df['n_genomes']==1) & (stats_table_df['paralog_flag']==0)]
        
        vir_count = [0,0]
        
        for stats_row in acquired_pangenes_stats.iterrows():
            pangene_row = pangenome_table_df[pangenome_table_df.index==stats_row[0]]
            orthologs = [pangene_row.index.values[0]]
            if pangene_row['ortholog_ids'].values[0]!=None:
                orthologs.extend(pangene_row['ortholog_ids'].values[0].split(';'))
            #annotations = []
            
            for ortholog in orthologs:
                stop_flag = False
                ortholog_loc = [int(l) for l in re.findall('\d+', ortholog.split('|')[-1])]
                strain_id = re.search('NC_\d{6}', ortholog.split('|')[3]).group(0)
                strain_row = strains_df[strains_df.index==strain_id]
                gb_filename = os.path.join(GBKDIR, '%s_%s' % (org_name, strain_row['strain_name_uid'].values[0]), '%s.gbk' % strain_id)
                with open(gb_filename, 'r') as gb_fh:
                    for record in SeqIO.parse(gb_fh, 'genbank'):
                        for feature in record.features:
                            if feature.type=='CDS' and feature.location.start in range(min(ortholog_loc)-1, min(ortholog_loc)+1) \
                            and feature.location.end in range(max(ortholog_loc)-1, max(ortholog_loc)+1):
                                ann_str = ''
                                if 'product' in feature.qualifiers.keys(): 
                                    ann_str = ';'.join(feature.qualifiers['product'])
                                for vir_keyword in VIR_KEYWORDS:
                                    if re.search(vir_keyword, ann_str, re.IGNORECASE):
                                        print strain_row['strain_name_uid'].values[0], ann_str, int(stats_row[1]['viral_flag'])
                                        vir_count[int(stats_row[1]['viral_flag'])]+=1
                                        print vir_count 
                                        stop_flag = True
                                        break
                                    '''if ann_str!='': 
                                    annotations.append(ann_str)
                                break'''
                                   
                        if stop_flag==True:
                            break
                if stop_flag==True:
                    break
            
            '''print annotations
            for ann in annotations:
                for vir_keyword in VIR_KEYWORDS:
                    if re.search(vir_keyword, ann, re.IGNORECASE):
                        print stats_row[1]['n_genomes'], ann
                        #strain_row['strain_name_uid'].values[0]
                        return'''
        print vir_count
        return


if __name__ == '__main__':
    for dir in os.listdir(OUTPUTDIR):
        if os.path.isdir(os.path.join(OUTPUTDIR, dir)):
            print dir
            acquired_genes_by_annotations(os.path.join(OUTPUTDIR, dir, 'output_datatables.db'), dir)
