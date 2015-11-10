import os, csv, re, subprocess, threading, sqlalchemy
from Bio import SeqIO
from collections import defaultdict
from time import sleep
import numpy as np
import pandas

FFNDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/all_ffn'
OUTPUTDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/output'
PHISPYDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/phispy'

THRESHOLDS = [0.1, 0.08, 0.05, 0.03, 0.02]
IDENTITYPCNT = 90
ALMINLENGTH = 0.9

ORG_NAMES = ['Bacillus_cereus', 'Campylobacter_jejuni', 'Clostridium_botulinum', 'Corynebacterium_diphtheriae',
             'Corynebacterium_pseudotuberculosis', 'Escherichia_coli-Shigella', 'Francisella_tularensis', 
             'Haemophilus_influenzae', 'Listeria_monocytogenes', 'Mycobacterium_tuberculosis', 'Neisseria_meningitidis',
             'Pseudomonas_aeruginosa', 'Streptococcus_pneumoniae', 'Yersinia_pestis']

class fasta_hit:
    '''object representing a fasta result
        result indices: 
        0-query id, 1-subject id, 2-% identity, 3-alignment length, 4-mismatches, 5-gap opens, 
        6-q. start, 7-q. end, 8-s. start, 9-s. end, 10-evalue, 11-bit score'''
    def __init__(self, hit_arr):
        self.query_id = hit_arr[0]
        self.subject_id = hit_arr[1]
        self.identity = float(hit_arr[2])
        self.al_length = int(hit_arr[3])
        self.q_length = abs(int(hit_arr[7])-int(hit_arr[6]))+1
        self.s_length = abs(int(hit_arr[9])-int(hit_arr[8]))+1
            
    def __repr__(self):
        return '''<query id: %s; subject id: %s; identity percent: %.2f; 
                alignment length: %d; query length: %d; subject length: %d>''' % (self.query_id, self.subject_id, self.identity, 
                                                                                      self.al_length, self.q_length, self.s_length)

def hit2obj(fasta_line):
    '''convert result from fasta output file to object'''
    if not re.match('\#', fasta_line):
        hit_arr = fasta_line.strip().split('\t')
        return fasta_hit(hit_arr)
    else:
        return None

def isunique(hit):
    '''check whether hit object is unique, according to predetermined parameters IDENTITYPCNT and ALNIMLENGTH'''
    if hit and (hit.identity<IDENTITYPCNT or hit.al_length<ALMINLENGTH*hit.q_length):
        return True
    return False

def run_fasta(q_path, lib_path, output_path, timeout, b=1, d=1, use_stdin=False):
    '''execute fasta from commandline, using stdin for query is turned off by default'''
    params = 'fasta36 -b =%s -d %s -m 8C -m 10' % (b, d)
    io = [q_path, lib_path, '>', output_path]
    
    if use_stdin:
        params = 'cat %s | fasta36 -b =%s -d %s -m 8C -m 10' % (q_path, b, d)
        io = ['@' , lib_path, '>', output_path]
    
    run_command_with_timeout(['/bin/bash', '-i', '-c', '%s %s' % (params, ' '.join(io))], timeout, output_path)    

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

def create_unified_genome(org_name, unifiables, strain_row, uni_row, strains_df, unified_dir):
    '''create one ffn file from ffn files for all unifiables, using strain as reference'''
    print strain_row.index
    strain_fasta_fname = os.path.join(FFNDIR, '_'.join([org_name, strain_row['strain_name_uid']]), '%s.ffn' % strain_row.index)
    
    unified_fname = os.path.join(unified_dir, '%s.ffn' % strain_row.index)
    if unifiables: unified_fname = os.path.join(unified_dir, '%s.ffn' % '_'.join([strain_row.index] + unifiables))
    
    if not os.path.isfile(unified_fname):
        unified_fh = open(unified_fname, 'a')
        
        SeqIO.write(SeqIO.parse(open(strain_fasta_fname, 'r'), 'fasta'), unified_fh, 'fasta')
        print sum([1 for record in SeqIO.parse(open(strain_fasta_fname, 'r'), 'fasta')])
        
        if unifiables:
            unifiable_record_ids = defaultdict(list)
            for unifiable in unifiables:
                unifiable_name = strains_df[strains_df.index==unifiable]
                print unifiable_name                
                unifiable_record_ids.update({unifiable: []})
                
                with open(os.path.join('fasta_out_fluidity', '%s_vs_%s.out' % (unifiable_name, strain_row.index)), 'r') as fasta_comp_fh:
                    for line in fasta_comp_fh:
                        hit = hit2obj(line)
                        if hit:
                            unique_flag = isunique(hit)
                            
                            if unifiables.index(unifiable)>0:
                                for prev_unifiable in unifiables[:unifiables.index(unifiable)]:
                                    prev_unifiable_name = strains_df[strains_df.index==prev_unifiable]['strain_name_uid']
                                    with open(os.path.join('fasta_out_fluidity', '%s_vs_%s.out' % (unifiable_name, prev_unifiable_name)), 'r') as fasta_prev_comp_fh:
                                        for prev_line in fasta_prev_comp_fh:
                                            prev_hit = hit2obj(prev_line)
                                            if prev_hit!=None and hit.query_id==prev_hit.query_id and prev_hit.subject_id in unifiable_record_ids[prev_unifiable]:
                                                unique_flag = False
                                                break
                                        if unique_flag==False:
                                            break
                            
                            if unique_flag:
                                unifiable_record_ids[unifiable].append(hit.query_id)
                                unifiable_fname = os.path.join(FFNDIR, '_'.join([org_name, unifiable_name]), 
                                                           '%s.ffn' % re.search('\|ref\|(NC_\d{6})\.\d', hit.query_id).group(1))
                                with open(unifiable_fname, 'r') as unifiable_fh: 
                                    with SeqIO.parse(unifiable_fh, 'fasta') as strain_fasta: 
                                        for record in strain_fasta:
                                            if hit.query_id==record.id:
                                                SeqIO.write(record, unified_fh, 'fasta')
                                                break
            
                print len(unifiable_record_ids[unifiable])
    
    #return unified_fname
    
def create_pangenome(org_name):
    '''create pangenome for org_name 
    output: ffn file containing a single entry for each pangene and sqlite table containing all pangenes and orthologs'''
    db_fname = os.path.join(OUTPUTDIR, org_name, 'output_datatables_%s.db' % org_name)
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    
    if not os.path.isdir('fasta_out_pangenome'): os.mkdir('fasta_out_pangenome')
    
    insp = sqlalchemy.inspect(engine)
    cols = insp.get_columns('unifiable_strains')
    thresholds = [c['name'] for c in cols]
    unifiables_df = pandas.read_sql_table('unifiable_strains', engine, index_col='strain_id')
    strains_df = pandas.read_sql_table('strains_plasmids', engine, index_col='strain_id')
    
    for threshold in thresholds: 
        threshold_str_sql = "%sEXPneg%d" % (re.search('([1-9])', threshold), threshold.count('0'))
        threshold_str_file = "%se-%d" % (re.search('([1-9])', threshold), threshold.count('0'))
        
        ans = 'y' #change to '' for dialog to appear
        
        while ans not in ['y', 'n']:
            print "would you like to create a pangenome file for fluidity threshold %s? [y/n]: " % threshold_str_file
            ans = raw_input()
        
        if ans=='y':            
            threshold_output = os.path.join('fasta_out_pangenome', threshold_str_file)
            if not os.path.isdir(threshold_output): os.mkdir(threshold_output)
            
            pangenome_fasta_fname = '%s_pangenome_identity%d_fluidity%s.ffn' % (org_name, IDENTITYPCNT, threshold_str_file)
            pangenome_table_name = "pangenome_{}".format(threshold_str_sql)
            
            unified_dir = os.path.join('fasta_out_pangenome', threshold_str_file, 'unified')
            if not os.path.isdir(unified_dir): os.mkdir(unified_dir)
            
            print os.path.abspath(unified_dir)
            
            #unifiables_table = engine.execute("SELECT strain_id, fluidity_{} FROM unifiable_strains;".format(threshold_str_sql)).fetchall()
            
            for uni_row in unifiables_df.iterrows():
                #strain_row = engine.execute("SELECT strain_name_uid FROM strains_plasmids WHERE strain_id=(?);", [uni_row[0]]).fetchone()
                strain_row = strains_df[strains_df.index==uni_row.index]
                
                unifiable_strains = None
                if uni_row[threshold]: unifiable_strains = uni_row[threshold].split(';')
                
                unified_fname = os.path.join(unified_dir, '%s.ffn' % strain_row.index)
                if unifiable_strains: unified_fname = os.path.join(unified_dir, '%s.ffn' % '_'.join([strain_row.index] + unifiable_strains))
                
                if not os.path.isfile(unified_fname): create_unified_genome(org_name, unifiable_strains, strain_row, uni_row, strains_df, unified_dir)
                    
                print strain_row.index, unified_fname, sum([1 for record in SeqIO.parse(open(unified_fname, 'r'), 'fasta')])
        
            pangenome_dir = os.path.join('fasta_out_pangenome', threshold_str_file, 'pangenome')
            if not os.path.isdir(pangenome_dir): os.mkdir(pangenome_dir)            
            
            for fname in os.listdir(unified_dir):
                print fname

                genes_count, unique_genes_count = 0, 0
                
                with open(os.path.join(unified_dir, fname), 'r') as genome_fh:
                    genome_records = list(SeqIO.parse(genome_fh, 'fasta'))
                
                    pangenome_records = []
                    
                    if os.path.isfile(pangenome_fasta_fname):
                        pangenome_records = list(SeqIO.parse(open(pangenome_fasta_fname, 'r'), 'fasta'))
                    
                    if not os.path.isfile(pangenome_fasta_fname) or os.stat(pangenome_fasta_fname).st_size==0:
                        with open(pangenome_fasta_fname, 'a') as pangenome_fasta_fh:
                            SeqIO.write(genome_records, pangenome_fasta_fh, 'fasta')
                        genes_count = sum([1 for record in SeqIO.parse(open(os.path.join(unified_dir, fname), 'r'), 'fasta')])
                        unique_genes_count = genes_count
                        pangenome_records.to_sql(pangenome_table_name, engine, if_exists='replace', index_label='pangenome_id')
                        for record in genome_records:
                            conn = engine.connect()
                            with conn.begin() as trans: 
                                try: 
                                    conn.execute('INSERT INTO {} VALUES (?, ?)'.format(pangenome_table_name), [record.id, None])
                                    trans.commit()
                                except:
                                    trans.rollback()
                                    raise
                            conn.close()
                    
                    elif False not in [genome_rec.id in [pangenome_rec.id for pangenome_rec in pangenome_records] for genome_rec in genome_records]:
                        genes_count = len(genome_records)
                                                                
                    else:
                        output_fasta_fname = os.path.join(pangenome_dir, '%s.out' % re.search('(.*)\.ffn', fname).group(1))
                        
                        if not os.path.isfile(output_fasta_fname):
                            run_fasta(os.path.join(unified_dir, fname), pangenome_fasta_fname, output_fasta_fname, 600)
                        if os.stat(output_fasta_fname).st_size==0:
                            run_fasta(os.path.join(unified_dir, fname), pangenome_fasta_fname, output_fasta_fname, 600, use_stdin=True)
                        if os.stat(output_fasta_fname).st_size==0:
                            print 'could not process query'
                            return
                            
                        if os.path.isfile(output_fasta_fname) and os.stat(output_fasta_fname).st_size>0:
                            unique_records = []
                            with open(output_fasta_fname, 'r') as output_fasta_fh:
                                for line in output_fasta_fh:
                                    hit = hit2obj(line)
                                    if hit:
                                        #print hit
                                        genes_count+=1
                                        unique_flag = isunique(hit)
                                        pangenome_keys = [k[0] for k in engine.execute("SELECT pangene_id from {}".format(pangenome_table_name)).fetchall()]
                                        if unique_flag and hit.query_id not in [p_rec.id for p_rec in pangenome_records]:
                                            unique_genes_count+=1
                                            if hit.query_id not in pangenome_keys:
                                                conn = engine.connect()
                                                with conn.begin() as trans: 
                                                    try: 
                                                        conn.execute('INSERT INTO {} VALUES (?, ?)'.format(pangenome_table_name), [hit.query_id, None])
                                                        trans.commit()
                                                    except:
                                                        trans.rollback()
                                                        raise
                                                conn.close()
                                                for record in genome_records:
                                                    if hit.query_id==record.id:
                                                        unique_records.append(record)
                                                        if len(record.seq)==0:
                                                            print record
                                                            return
                                                        break
                                        else:
                                            pangenome_entry = engine.execute('SELECT ortholog_ids FROM {} WHERE pangene_id=(?)'.format(pangenome_table_name), 
                                                                          [hit.subject_id]).fetchone()
                                            
                                            if pangenome_entry not in [None, (None,)]:
                                                pangenome_entry = pangenome_entry[0].split(';')
                                            else:
                                                pangenome_entry = []
                                            if hit.query_id not in pangenome_entry:
                                                pangenome_entry.append(hit.query_id)
                                                conn = engine.connect()
                                                with conn.begin() as trans: 
                                                    try: 
                                                        conn.execute('UPDATE {} SET ortholog_ids=(?) WHERE pangene_id=(?)'.format(pangenome_table_name), 
                                                                    [';'.join(pangenome_entry), hit.subject_id])
                                                        trans.commit()
                                                    except:
                                                        trans.rollback()
                                                        raise
                                                conn.close()
                                
                                with open(pangenome_fasta_fname, 'a') as pangenome_fasta_fh:
                                    SeqIO.write(unique_records, pangenome_fasta_fh, 'fasta')
                            
                    print genes_count, unique_genes_count, len(engine.execute("SELECT pangene_id from {}".format(pangenome_table_name)).fetchall())
            
        elif ans=='n':
            continue
        
    return

def recover_lost_seqs(org_name, fasta_filename, strains_info_dict):
    fasta_updated = '%s_corrected.ffn' % fasta_filename.split('.')[0]
    if not os.path.isfile(fasta_updated):
        fasta_u_fh = open(fasta_updated, 'a')
        for record in SeqIO.parse(open(fasta_filename, 'r'), 'fasta'):
            if len(record.seq)==0:
                for strain_dict in strains_info_dict:
                    if re.search(strain_dict['keys'], record.id):
                        st_dir = '_'.join([org_name, strain_dict['strain'], 'uid%s' % strain_dict['uid']])
                        ffn_fname = os.path.join(FFNDIR, st_dir, '%s.ffn' % strain_dict['keys'])
                        for ffn_rec in SeqIO.parse(open(ffn_fname, 'r'), 'fasta'):
                            if record.id==ffn_rec.id:
                                record.seq = ffn_rec.seq
                                break
                        break
            SeqIO.write(record, fasta_u_fh, 'fasta')
    return fasta_updated

def self_compare(org_name):
    '''compare pangenome to itself, unite entries similar above alignment similarity threshold'''
    db_fname = os.path.join(OUTPUTDIR, org_name, 'output_datatables.db')
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
        
    for threshold in THRESHOLDS: 
        threshold_str_sql = '%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        threshold_str_file = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        
        print threshold_str_file
        output_dir = os.path.join('fasta_out_pangenome', threshold_str_file)
        if not os.path.isdir(output_dir): 
            raise OSError('output folder for fluidity threshold %s does not exist')
         
        output_fasta_fname = os.path.join(output_dir, 'pangenome_identity%d_fluidity%s.out' % (IDENTITYPCNT, threshold_str_file))
        
        pangenome_records = []
        
        fasta_fname = 'pangenome_identity%d_fluidity%s.ffn' % (IDENTITYPCNT, threshold_str_file)
        fasta_after_selfcmp_fname = 'pangenome_identity%d_fluidity%s_selfcmp.ffn' % (IDENTITYPCNT, threshold_str_file)
        
        pangenome_table_name = 'pangenome_{}'.format(threshold_str_sql)
        
        if os.path.isfile(fasta_fname):
            with open(fasta_fname, 'r') as pangenome_fh:
                for record in SeqIO.parse(pangenome_fh, 'fasta'):
                    pangenome_records.append(record)
        
        n_pangenes = engine.execute('SELECT Count(*) FROM {}'.format(pangenome_table_name)).fetchone()[0]
        print len(pangenome_records), n_pangenes
        
        if not os.path.isfile(output_fasta_fname):
            run_fasta(fasta_fname, fasta_fname, output_fasta_fname, 1800, b=2, d=2)
        if os.stat(output_fasta_fname).st_size==0:
            print 'could not process query'
            return
        
        if os.path.isfile(output_fasta_fname) and os.stat(output_fasta_fname).st_size>0: 
            with open(output_fasta_fname, 'r') as output_fasta_fh:    
                for line in output_fasta_fh:
                    hit = hit2obj(line)
                    if hit:
                        if hit.query_id!=hit.subject_id:
                            unique_flag = isunique(hit)
                            if not unique_flag:
                                query_entry = engine.execute('SELECT * FROM {} WHERE pangene_id=(?)'.format(pangenome_table_name), [hit.query_id]).fetchone()
                                subject_entry = engine.execute('SELECT * FROM {} WHERE pangene_id=(?)'.format(pangenome_table_name), [hit.subject_id]).fetchone()
                                if query_entry not in [None, (None,)] and subject_entry not in [None, (None,)]: 
                                    query_orthologs, subject_orthologs = [], []
                                    if query_entry[1] not in [None, (None,)]:
                                        query_orthologs = query_entry[1].split(';')
                                    if subject_entry[1] not in [None, (None,)]:
                                        subject_orthologs = subject_entry[1].split(';')
                                    entries = list(set(query_orthologs + [hit.query_id] + subject_orthologs))
                                    engine.execute('UPDATE {} SET ortholog_ids=(?) WHERE pangene_id=(?)'.format(pangenome_table_name), [';'.join(entries), 
                                                                                                                                    hit.subject_id])
                                    
                                    conn = engine.connect()
                                    with conn.begin() as trans: 
                                        try: 
                                            conn.execute('DELETE FROM {} WHERE pangene_id=(?)'.format(pangenome_table_name), [hit.query_id]) 
                                            trans.commit()
                                        except:
                                            trans.rollback()
                                            raise
                                    conn.close()
                                    
                                    rec_to_pop = None
                                    for rec in pangenome_records:
                                        if rec.id==hit.query_id:
                                            rec_to_pop = rec
                                            break
                                    if rec_to_pop!=None:
                                        pangenome_records.pop(pangenome_records.index(rec_to_pop))
            
            with open(fasta_after_selfcmp_fname, 'w') as pangenome_fh_updated:
                SeqIO.write(pangenome_records, pangenome_fh_updated, 'fasta')
            
            n_pangenes_updated = engine.execute('SELECT Count(*) FROM {}'.format(pangenome_table_name)).fetchone()[0]            
            print len(pangenome_records), n_pangenes_updated
        
    return

def find_paralogs(org_name):
    '''find paralogs in genome, export gene names to paralogs table'''
    strains_fh = open(os.path.join(OUTPUTDIR, org_name, 'strains_plasmids_ncbiFTP.csv'))
    strains_dict = [row for row in csv.DictReader(strains_fh)]
    strains_fh.seek(0)
    strains_headers = csv.DictReader(strains_fh).fieldnames
    
    fasta_output_dir = os.path.join(OUTPUTDIR, org_name, 'fasta_out_paralogs')
    if not os.path.isdir(fasta_output_dir): os.mkdir(fasta_output_dir)
    
    paralogs_dict = defaultdict(list)
    
    tsv_paralogs_fname = os.path.join(OUTPUTDIR, org_name, 'paralogs.tsv')
    if os.path.isfile(tsv_paralogs_fname):
            for row in csv.reader(open(tsv_paralogs_fname, 'r'), delimiter='\t'):
                if len(row)==2:
                    paralogs_dict.update({row[0]: row[1].split(';')})
        
    print len(paralogs_dict.keys())
    
    for strain_row in strains_dict:
        print strain_row['strain']
        strain_dir = os.path.join(FFNDIR, '_'.join([org_name, strain_row['strain'], 'uid%s' % strain_row['uid']]))
        strain_ffn_fname = '%s.ffn' % strain_row['keys']
        
        strain_fasta_output = os.path.join(fasta_output_dir, '%s.out' % strain_row['keys'])
        os.chdir(strain_dir)
        
        if not os.path.isfile(strain_fasta_output): 
            run_fasta(strain_ffn_fname, strain_ffn_fname, strain_fasta_output, 600, b=2, d=2)
        
        if os.stat(strain_fasta_output).st_size==0:
            print 'could not process query'
            return
        
        if os.stat(strain_fasta_output).st_size>0:
            genes_count, paralogs_count = 0, 0
            for line in open(strain_fasta_output, 'r'):
                hit = hit2obj(line)
                if hit:
                    if hit.query_id!=hit.subject_id:
                        unique_flag = isunique(hit)
                        if not unique_flag:
                            print hit
                            return
                            paralogs_count+=1
                            if hit.subject_id not in paralogs_dict.keys():
                                paralogs_dict.update({hit.subject_id: [hit.query_id]})
                            else:
                                paralogs_dict[hit.subject_id].append(hit.query_id)
                
                elif re.search('# FASTA processed \d+ queries', line):
                    genes_count = re.search('# FASTA processed (\d+) queries', line).group(1)
            
            print genes_count, paralogs_count, len(paralogs_dict.keys()) 
            
            tsv_selfcmp_writer = csv.writer(open(tsv_paralogs_fname, 'w'), delimiter='\t')
            for key in paralogs_dict.keys():
                tsv_selfcmp_writer.writerow([key, ';'.join(paralogs_dict[key])])
    
    return

def generate_pangenome_stats(org_name):
    '''create stats table for each pangenome, generated per fluidity threshold'''
    db_fname = os.path.join(OUTPUTDIR, org_name, 'output_datatables.db')
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    for threshold in THRESHOLDS: 
        threshold_str_sql = '%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        pangenome_table_name = 'pangenome_%s' % threshold_str_sql
        stats_table_name = 'pangenome_%s_stats' % threshold_str_sql
        
        print threshold_str_sql
        
        pangenes = engine.execute('SELECT * FROM {}'.format(pangenome_table_name)).fetchall()
        for pangene in pangenes:
            '''calculate number of orthologs'''
            orthologs = [pangene[0]]
            if pangene[1] not in [None, (None,)]:
                orthologs+=pangene[1].split(';')
            n_orthologs = len(orthologs)
            
            '''calculate number of genomes'''
            genomes = [re.search('\|(NC_\d{6})\.\d\|', o).group(1) for o in orthologs]
            unique_genomes = []
            for g in genomes:
                unifiables = engine.execute('SELECT fluidity_{} FROM unifiable_strains WHERE strain_id=(?)'.format(threshold_str_sql), [g]).fetchone()[0]
                if unifiables not in [None, (None,)]:
                    unifiables = unifiables.split(';')
                else:
                    unifiables = ''
                if g not in unique_genomes and True not in [u in unique_genomes for u in unifiables]:
                    unique_genomes.append(g)
            n_genomes = len(unique_genomes)
            if n_genomes==0:
                print pangene
            paralog_flag, viral_flag = False, False
            enc_arr, encprime_arr = [], []
            fop_arr = []
            
            for o in orthologs:
                genome = genomes[orthologs.index(o)]
                st_name = engine.execute('SELECT strain_name_uid FROM strains_plasmids WHERE strain_id=(?)', [genome]).fetchone()[0]
                
                '''check whether pangene is a paralog'''
                paralogs = engine.execute('SELECT * FROM paralogs WHERE pangene_id=(?)', [o]).fetchone()
                if paralogs not in [None, (None,)]:
                    paralog_flag = True
                
                '''check whether pangene is viral'''
                phispy_output_fname = os.path.join(PHISPYDIR, org_name, '%s_%s' % (org_name, st_name), genome, 'prophage_tbl.txt')
                if os.path.isfile(phispy_output_fname):
                    with open(phispy_output_fname, 'r') as phispy_output_fh:
                        for line in phispy_output_fh:
                            '''3-start, 4-stop, ..., 9-final_status'''
                            res = line.strip().split('\t')
                            o_arr = o.split('|')
                            pos = re.findall('\d+', o_arr[-1])
                            if res[3] in pos and res[4] in pos:
                                if res[9]=='1':
                                    viral_flag = True
                
                '''get ENC, ENCprime values'''
                encprime_fname = os.path.join(OUTPUTDIR, org_name, 'encprime', '%s_%s' % (org_name, st_name), '%s.ffn.results' % genome)
                with open(encprime_fname, 'r') as encprime_fh:
                    for ln in encprime_fh:
                        rw = ln.strip()
                        if o==rw.split(' ')[0]:
                            '''0-gene_name; 1-ENC; 2-ENCprime; ...'''
                            nums = re.findall('\d+\.\d+', rw)
                            enc_arr.append(float(nums[1]))
                            encprime_arr.append(float(nums[2]))
                            break
                
                '''get FOP values'''
                fop_fname = os.path.join(OUTPUTDIR, org_name, 'fop', '%s_%s' % (org_name, st_name), '%s.fop.out' % genome)
                with open(fop_fname, 'r') as fop_fh:
                    for line in fop_fh:
                        '''0-gene_name; 1-gene_length; 2-n_optimal_codons; 3-n_optimal_nonoptimal_codons; 5-fop'''
                        fop_row = line.strip().split('\t')
                        if fop_row[0]==o:
                            fop_arr.append(float(fop_row[4]))
                            break
            
            '''update stats table'''
            conn = engine.connect()
            with conn.begin() as trans: 
                try: 
                    conn.execute('INSERT INTO {} VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)'.format(stats_table_name),
                        [pangene[0], n_orthologs, n_genomes, int(paralog_flag), int(viral_flag),
                         np.average(enc_arr), np.std(enc_arr), np.average(encprime_arr), np.std(encprime_arr),
                         np.average(fop_arr), np.std(fop_arr)])
                    trans.commit()
                except:
                    trans.rollback()
                    raise
            conn.close()
                
    return

if __name__ == '__main__':
    for name in ORG_NAMES:
        if not os.path.isdir(os.path.join(OUTPUTDIR, name)):
            os.mkdir(os.path.join(OUTPUTDIR, name))
        print name
        os.chdir(os.path.join(OUTPUTDIR, name))
        create_pangenome(name)
        #self_compare(name)
        #find_paralogs(name)
        #generate_pangenome_stats(name)
            