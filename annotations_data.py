import pandas
import sqlite3
import subprocess, threading
from time import sleep
import numpy as np
import sqlalchemy, sqlalchemy.orm
import os, re
from Bio import SeqIO
from collections import defaultdict
from fuzzywuzzy import fuzz

BASEDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/'
OUTPUTDIR = os.path.join(BASEDIR, 'output')
GBKDIR = os.path.join(BASEDIR, 'all_gbk')
FFNDIR = os.path.join(BASEDIR, 'all_ffn')

THRESHOLDS = [0.1, 0.08, 0.05, 0.03, 0.02]
IDENTITYPCNT = 90
ALNIMLENGTH = 0.9

VIR_KEYWORDS = ['viral', 'virus', 'phage', 'capsid', 'tail']
UNKNOWN_KEYWORDS = ['hypothetical protein', 'uncharacterized protein']
ORIGIN_VALS = ['bacterial', 'unknown', 'viral']

def run_command_with_timeout(cmd, timeout_sec, output_path):
    """Execute `cmd` in a subprocess and enforce timeout `timeout_sec` seconds.
 
    Return subprocess exit code on natural completion of the subprocess.
    Raise an exception if timeout expires before subprocess completes."""
    proc = subprocess.Popen(cmd)
    print subprocess.list2cmdline(cmd)
    sleep(1)
    proc_thread = threading.Thread(target=proc.communicate)
    proc_thread.start()
    proc_thread.join(timeout_sec)
    '''if no progress has been made by timeout expiration, kill process'''
    if proc_thread.is_alive() and os.stat(output_path).st_size==0:
        # Process still running - kill it and raise timeout error
        try:
            proc.kill()
        except OSError, e:
            # The process finished between the `is_alive()` and `kill()`
            return proc.returncode
        '''# OK, the process was definitely killed
        raise SubprocessTimeoutError('Process #%d killed after %f seconds' % (proc.pid, timeout_sec))'''
    # Process completed naturally - return exit code
    return proc.returncode

def create_functional_annotations_dict(go_textfile):
    with open(go_textfile, 'r') as go_fh:
        func_dict = defaultdict(list)
        entry = []
        for line in go_fh.xreadlines():
            line = line.strip('\n')
            if re.match('//', line):
                entry = []
            if re.match('ID', line):
                entry.append(re.match('ID   (.*)\.$', line).group(1))
            if re.match('SY', line):
                entry.append(re.match('SY   ([^.]*)\.?$', line).group(1))
            if re.match('CA', line):
                entry.append(re.match('CA   (.*)\.$', line).group(1))
                if len(entry)>1:
                    if entry[-1] in func_dict.keys():
                        func_dict[entry[-1]].append('; '.join(entry[0:-1]))
                    else:
                        func_dict.update({entry[-1]: ['; '.join(entry[0:-1])]})
    
    return func_dict
 
def get_gene_annotations(db_fname, org_name):
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    
    for threshold in THRESHOLDS:
        threshold_str_file = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        print threshold_str_file
        
        threshold_str_sql = '%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        pangenome_table_name, stats_table_name = 'pangenome_%s' % threshold_str_sql, 'pangenome_%s_stats' % threshold_str_sql
        
        col_data = engine.execute('PRAGMA table_info(gene_annotations)').fetchall()
        
        if True not in [col[1]=='origin' for col in col_data]:
            engine.execute('ALTER TABLE {} ADD COLUMN origin TEXT'.format(stats_table_name))
                
        pangenome_table_df = pandas.read_sql_table(pangenome_table_name, engine, index_col='pangene_id')          
        stats_table_df = pandas.read_sql_table(stats_table_name, engine, index_col='pangene_id')
        
        if True in [val not in ORIGIN_VALS for val in stats_table_df['origin'].values]:
            for stats_row in stats_table_df.iterrows():
                if stats_row[1]['origin'] not in ORIGIN_VALS:
                    pangene_row = pangenome_table_df[pangenome_table_df.index==stats_row[0]]
                    orthologs = [pangene_row.index.values[0]]
                    if pangene_row['ortholog_ids'].values[0]!=None:
                        orthologs.extend(pangene_row['ortholog_ids'].values[0].split(';'))
                    annotations = []
                    
                    for ortholog in orthologs:
                        print ortholog
                        o_query = ortholog
                        if re.search(',', ortholog):
                            pos_arr = re.findall('\d+', ortholog.split('|')[-1])
                            final_sep = '|:'
                            pos_str = '%d-%d' % (min([int(p) for p in pos_arr]), max([int(p) for p in pos_arr]))
                            if re.search('c', ortholog):
                                final_sep='|:c'
                                pos_str = '%d-%d' % (max([int(p) for p in pos_arr]), min([int(p) for p in pos_arr]))
                            o_query = '|'.join(ortholog.split('|')[0:-1])+final_sep+pos_str
                        print o_query    
                        
                        gene_annotation = engine.execute('SELECT annotation FROM gene_annotations WHERE gene_id=(?)', (o_query,)).fetchone()
                        if gene_annotation==None:
                            pos_arr = re.findall('c?\d+', ortholog.split('|')[-1])
                            o_query = '|'.join(ortholog.split('|')[0:-1])+'|:'+pos_arr[0]+'->'+pos_arr[-1]
                            gene_annotation = engine.execute('SELECT annotation FROM gene_annotations WHERE gene_id=(?)', (o_query,)).fetchone()
                        
                        if gene_annotation==None:
                            pos_arr = re.findall('c?\d+', ortholog.split('|')[-1])
                            if re.search('c', pos_arr[0]): pos_arr[0] = 'c>%s' % re.search('(\d+)', pos_arr[0]).group(1)
                            o_query = '|'.join(ortholog.split('|')[0:-1])+'|:'+pos_arr[0]+'-'+pos_arr[-1]
                            gene_annotation = engine.execute('SELECT annotation FROM gene_annotations WHERE gene_id=(?)', (o_query,)).fetchone()
                        
                        #print o_query
                        if gene_annotation[0]=='viral':
                            stats_table_df.loc[stats_row[0],'origin'] = 'viral'
                            break
                        else:
                            annotations.append(gene_annotation[0])
                    
                    if len(annotations)==len(orthologs):
                        if annotations==['unknown']*len(annotations):
                            stats_table_df.loc[stats_row[0],'origin'] = 'unknown'
                        else:
                            stats_table_df.loc[stats_row[0],'origin'] = 'bacterial'
                    
                    print stats_table_df.loc[stats_row[0]]
                    if stats_table_df.loc[stats_row[0], 'origin'] not in ORIGIN_VALS: return
                    
                    conn = engine.connect()
                    with conn.begin() as trans: 
                        try: 
                            conn.execute('UPDATE {} SET origin=(?) WHERE pangene_id=(?)'.format(stats_table_name), (stats_table_df.loc[stats_row[0], 'origin'], stats_row[0]))
                            trans.commit()
                        except:
                            trans.rollback()
                            raise
                    conn.close()
        
    return
    
def assign_functional_category(db_fname):
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    
    for threshold in THRESHOLDS:
        threshold_str_file = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        print threshold_str_file
        
        threshold_str_sql = '%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        pangenome_table_name, stats_table_name = 'pangenome_%s' % threshold_str_sql, 'pangenome_%s_stats' % threshold_str_sql
        
        col_data = engine.execute('PRAGMA table_info({})'.format(stats_table_name)).fetchall()
        
        if [re.match('annotation', col[1]) for col in col_data]==[None]*len(col_data):
            engine.execute('ALTER TABLE {} ADD COLUMN annotation TEXT'.format(stats_table_name))
        
        pangenome_table_df = pandas.read_sql_table(pangenome_table_name, engine, index_col='pangene_id')          
        stats_table_df = pandas.read_sql_table(stats_table_name, engine, index_col='pangene_id')
        
        #if True in [val==None for val in stats_table_df['annotation'].values]:
        for stats_row in stats_table_df.iterrows():
            #if not stats_row[1]['annotation']:
            pangene_row = pangenome_table_df[pangenome_table_df.index==stats_row[0]]
            orthologs = [pangene_row.index.values[0]]
            if pangene_row['ortholog_ids'].values[0]!=None:
                orthologs.extend(pangene_row['ortholog_ids'].values[0].split(';'))
            annotations = {}
            
            for ortholog in orthologs:
                #print ortholog
                o_query = ortholog
                if re.search(',', ortholog):
                    pos_arr = re.findall('\d+', ortholog.split('|')[-1])
                    final_sep = '|:'
                    pos_str = '%d-%d' % (min([int(p) for p in pos_arr]), max([int(p) for p in pos_arr]))
                    if re.search('c', ortholog):
                        final_sep='|:c'
                        pos_str = '%d-%d' % (max([int(p) for p in pos_arr]), min([int(p) for p in pos_arr]))
                    o_query = '|'.join(ortholog.split('|')[0:-1])+final_sep+pos_str
                #print o_query    
                
                gene_annotation = engine.execute('SELECT annotation_full FROM gene_annotations WHERE gene_id=(?)', (o_query,)).fetchone()
                if gene_annotation==None:
                    pos_arr = re.findall('c?\d+', ortholog.split('|')[-1])
                    o_query = '|'.join(ortholog.split('|')[0:-1])+'|:'+pos_arr[0]+'->'+pos_arr[-1]
                    gene_annotation = engine.execute('SELECT annotation_full FROM gene_annotations WHERE gene_id=(?)', (o_query,)).fetchone()
                
                if gene_annotation==None:
                    pos_arr = re.findall('c?\d+', ortholog.split('|')[-1])
                    if re.search('c', pos_arr[0]): pos_arr[0] = 'c>%s' % re.search('(\d+)', pos_arr[0]).group(1)
                    o_query = '|'.join(ortholog.split('|')[0:-1])+'|:'+pos_arr[0]+'-'+pos_arr[-1]
                    gene_annotation = engine.execute('SELECT annotation_full FROM gene_annotations WHERE gene_id=(?)', (o_query,)).fetchone()
                
                #print o_query, gene_annotation
                annotations_str = gene_annotation[0].lower()
                                        
                if annotations_str not in annotations.keys():
                    annotations.update({annotations_str: 1})
                    '''if True not in [(fuzz.partial_ratio(annotations_str, key)>84 or fuzz.partial_ratio(key, annotations_str)>84) for key in annotations.keys()]:
                        annotations.update({annotations_str: 1})
                    else:
                        for key in annotations.keys():
                            if fuzz.partial_ratio(annotations_str, key)>84 or fuzz.partial_ratio(key, annotations_str)>84:
                                #print annotations_str, '|', key
                                if (len(key)<len(annotations_str) and not (re.match('putative', annotations_str, re.IGNORECASE) 
                                                                          or re.match('predict', annotations_str, re.IGNORECASE))) \
                                or (re.match('putative', key, re.IGNORECASE) or re.match('predict', key, re.IGNORECASE)):
                                    annotations.update({annotations_str: annotations[key]+1})
                                    try:
                                        annotations.pop(key)#, None)
                                    except KeyError:
                                        raise
                                else: annotations[key]+=1'''
                else:
                    annotations[annotations_str]+=1
            
            tot_counter = float(sum([val for val in annotations.values()]))
            largest_frac = ['hypothetical protein', 0]
            for key in annotations.keys():
                if key not in UNKNOWN_KEYWORDS:
                    key_frac = float(annotations[key])/tot_counter
                    if key_frac > largest_frac[1] or ([re.match(pred_word, largest_frac[0]) for pred_word in ['putative', 'predicted']]!=[None]*2 and
                                                      [re.match(pred_word, key) for pred_word in ['putative', 'predicted']]==[None]*2):
                        largest_frac = [key, key_frac]
                elif annotations[key]==tot_counter:
                    largest_frac = [key, 1]
            
            if stats_table_df.loc[stats_row[0], 'annotation']!=largest_frac[0]:   
                print stats_table_df.loc[stats_row[0], 'annotation'], '==>', largest_frac[0], '\n\t', annotations
                change_ann = raw_input('would you like to change the functional annotation for the gene? [y/n]') 
                update_str = ''
                if change_ann=='y':
                    update_str = largest_frac[0]
                    #stats_table_df.loc[stats_row[0], 'annotation'] = largest_frac[0]
                    #print stats_table_df.loc[stats_row[0]]
                    input_str = raw_input('alternative annotation [or d for default]:')
                    if input_str!='d':
                        update_str = input_str
                    '''if re.match('putative', largest_frac[0], re.IGNORECASE) or re.match('predict', largest_frac[0], re.IGNORECASE):
                        return
                    '''
                    conn = engine.connect()
                    with conn.begin() as trans: 
                        try: 
                            conn.execute('UPDATE {} SET annotation=(?) WHERE pangene_id=(?)'.format(stats_table_name), 
                                         (update_str, stats_row[0]))
                            trans.commit()
                        except:
                            trans.rollback()
                            raise
                    conn.close()
        
    return

def get_heg(db_fname, org_name):
    con = sqlite3.connect(db_fname)
    with con:
        def re_fn(expression, item):
            reg = re.compile(expression, re.IGNORECASE)
            return reg.search(item) is not None
        
        con.create_function("REGEXP", 2, re_fn)
        heg_table = pandas.read_sql_query("SELECT * FROM gene_annotations WHERE annotation_full REGEXP 'ribosom'", con, index_col='gene_id')
        #strains_table = pandas.read_sql_query("SELECT * FROM strains_plasmids", con, index_col='strain_id')
        unifiables_table = pandas.read_sql_query("SELECT * FROM unifiable_strains", con, index_col='strain_id')
        
        for i in range(1, len(THRESHOLDS)+1):
            print THRESHOLDS[i-1]
            threshold_str_file = '%se-%d' % (re.search('([1-9])', str(THRESHOLDS[i-1])).group(1), str(THRESHOLDS[i-1]).count('0'))
            if not os.path.isdir(os.path.join(OUTPUTDIR, org_name, 'cai', 'heg', threshold_str_file)):
                os.mkdir(os.path.join(OUTPUTDIR, org_name, 'cai', 'heg', threshold_str_file))
                
            uni_col = 'fluidity_%sEXPneg%d' % (re.search('([1-9])', str(THRESHOLDS[i-1])).group(1), str(THRESHOLDS[i-1]).count('0'))
            uni_series =  unifiables_table[uni_col]
            for j, uni_strain in enumerate(uni_series):
                if True not in [re.search(uni_series.index[j], f)!=None for f in os.listdir(os.path.join(OUTPUTDIR, org_name, 'cai', 'heg', threshold_str_file))]:
                    heg_records = []
                    strain_ids = [uni_series.index[j]]
                    if uni_strain!=None: strain_ids+=uni_strain.split(';')
                    print strain_ids
                    heg_not_in_unified_fasta = 0
                    ffn_filename = ''
                    for fname in os.listdir(os.path.join(OUTPUTDIR, org_name, 'fasta_out_pangenome', threshold_str_file, 'unified')):
                        #print fname
                        if re.search(uni_series.index[j], fname):
                            ffn_filename = os.path.join(OUTPUTDIR, org_name, 'fasta_out_pangenome', threshold_str_file, 'unified', fname)
                            #print ffn_filename
                            break
                    if ffn_filename=='':
                        raise IOError('no fasta file for %s in directory %s' % (uni_series.index[j], 
                                                                                os.path.join(OUTPUTDIR, org_name, 'fasta_out_pangenome', threshold_str_file, 'unified')))
                    for heg_id in heg_table.index:
                        if True in [re.search(st_id, heg_id)!=None for st_id in strain_ids]:
                            #print heg_table.loc[heg_id]
                            found_flag = False
                            with open(ffn_filename, 'r') as ffn_fh:
                                for record in SeqIO.parse(ffn_fh, 'fasta'):
                                    if record.id==heg_id and record not in heg_records:
                                        heg_records.append(record)
                                        found_flag = True
                                        break
                            if not found_flag:
                                heg_not_in_unified_fasta+=1
                    output_filename = os.path.join(OUTPUTDIR, org_name, 'cai', 'heg', threshold_str_file, '%s.ffn' %  '_'.join(strain_ids))
                    with open(output_filename, 'w') as output_fh:
                        SeqIO.write(heg_records, output_fh, 'fasta')
                    print '%d records written to %s;\n%d heg entries not found in unified fasta file' % (len(heg_records), output_filename, heg_not_in_unified_fasta)
                    #return
    return

def calculate_cai(org_name):
    for threshold in THRESHOLDS:
        threshold_str_file = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        if not os.path.isdir(os.path.join(OUTPUTDIR, org_name, 'cai', 'cut', threshold_str_file)):
            os.mkdir(os.path.join(OUTPUTDIR, org_name, 'cai', 'cut', threshold_str_file))
        if not os.path.isdir(os.path.join(OUTPUTDIR, org_name, 'cai', 'cai_calc', threshold_str_file)):
            os.mkdir(os.path.join(OUTPUTDIR, org_name, 'cai', 'cai_calc', threshold_str_file))
        for heg_file in os.listdir(os.path.join(OUTPUTDIR, org_name, 'cai', 'heg', threshold_str_file)):
            heg_fname = re.search('([^.]+)\.ffn', heg_file).group(1)
            '''print heg_fname
            return'''
            pangenome_fasta = os.path.join(OUTPUTDIR, org_name, 'fasta_out_pangenome', threshold_str_file, 'unified', heg_file)
            os.chdir(os.path.join(OUTPUTDIR, org_name, 'cai', 'heg', threshold_str_file))
            cut_file = os.path.join(OUTPUTDIR, org_name, 'cai', 'cut', threshold_str_file, '%s.cut' % heg_fname)
            cai_file = os.path.join(OUTPUTDIR, org_name, 'cai', 'cai_calc', threshold_str_file, '%s.cai' % heg_fname)
            run_command_with_timeout(['/bin/bash', '-i', '-c', 'cusp %s %s' % (heg_file, cut_file)], 30, cut_file)
            run_command_with_timeout(['/bin/bash', '-i', '-c', 'cai -seqall %s -cfile %s -outfile %s' % (pangenome_fasta, cut_file, cai_file)], 30, cai_file)
    return

if __name__ == '__main__':
    #f_dict = create_functional_annotations_dict(os.path.join(BASEDIR, 'go_keywlist.txt'))
    
    for dir in os.listdir(OUTPUTDIR):
        if os.path.isdir(os.path.join(OUTPUTDIR, dir)):
            print dir
            #get_gene_annotations(os.path.join(OUTPUTDIR, dir, 'output_datatables.db'), dir)
            assign_functional_category(os.path.join(OUTPUTDIR, dir, 'output_datatables.db'))
            #get_heg(os.path.join(OUTPUTDIR, dir, 'output_datatables.db'), dir)
            #calculate_cai(dir)
    