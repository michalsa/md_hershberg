from Bio import SeqIO, SeqRecord, SeqFeature
import csv 
import re
import itertools
from scipy import stats
import subprocess, os, sys, threading
import pandas
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

IDENTITYPCNT = 90
ALMINLENGTH = 0.9

GBKDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/all_gbk'
FFNDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/all_ffn'
OUTPUTDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/output'

SEEDDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/seed'
PHISPYDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/phispy'

ENCPATH = '/home/michalsa/Downloads/ENCprime-master/bin'

AADICT = {'ala': ['gct', 'gcc', 'gca', 'gcg'],
          'arg': ['cgt', 'cgc', 'cga', 'cgg', 'aga', 'agg'],
          'asn': ['aat', 'aac'],
          'asp': ['gat', 'gac'],
          'cys': ['tgt', 'tgc'],
          'gln': ['caa', 'cag'],
          'glu': ['gaa', 'gag'],
          'gly': ['ggt', 'ggc', 'gga', 'ggg'],
          'his': ['cat', 'cac'],
          'ile': ['att', 'atc', 'ata'],
          'leu': ['tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg'],
          'lys': ['aaa', 'aag'],
          'phe': ['ttt', 'ttc'],
          'pro': ['cct', 'ccc', 'cca', 'ccg'],
          'ser': ['tct', 'tcc', 'tca', 'tcg', 'agt', 'agc'],
          'thr': ['act', 'acc', 'aca', 'acg'],
          'tyr': ['tat', 'tac'],
          'val': ['gtt', 'gtc', 'gta', 'gtg']}

OPTDICT = {'Escherichia_coli': {'ala': ['gcg'],
                                'arg': ['cgt'],
                                'asn': ['aac'],
                                'asp': ['gac'],
                                'cys': ['tgc'],
                                'gln': ['cag'],
                                'glu': ['gaa'],
                                'gly': ['ggc'],
                                'his': ['cac'],
                                'ile': ['atc'],
                                'leu': ['ctg'],
                                'lys': ['aaa'],
                                'phe': ['ttc'],
                                'pro': ['ccg'],
                                'ser': ['tcc', 'agc'],
                                'thr': ['acc'],
                                'tyr': ['tac'],
                                'val': ['gtg']}}

chunk_size = 1024
timeout_sec = 3000

csv.field_size_limit(sys.maxsize)

class SubprocessTimeoutError(RuntimeError):
    pass

def hit2dict(hit):
    '''result indices: 
       0-query id, 1-subject id, 2-% identity, 3-alignment length, 4-mismatches, 5-gap opens, 
       6-q. start, 7-q. end, 8-s. start, 9-s. end, 10-evalue, 11-bit score'''
    hit_arr = hit.strip().split('\t')
    hit_dict = {'query_id': hit_arr[0], 
                'subject_id': hit_arr[1],
                'identity%': float(hit_arr[2]), 
                'al_length': int(hit_arr[3]), 
                'q_length': abs(int(hit_arr[7])-int(hit_arr[6]))+1, 
                's_length': abs(int(hit_arr[9])-int(hit_arr[8]))+1}
    return hit_dict

def run_fasta(q_path, lib_path, output_path, use_stdin=False):
    params = 'fasta36 -b =1 -d 1 -m 8C -m 10'
    io = [q_path, lib_path, '>', output_path]
    
    if use_stdin==True:
        params = 'cat %s | ' % q_path + params + '-n'
        io[0] = '@'
    
    run_command_with_timeout(['/bin/bash', '-i', '-c', '%s %s' % (params, ' '.join(io))], 300, output_path)    
    '''
    sp = subprocess.Popen(['/bin/bash', '-i', '-c', '%s %s' % (params, ' '.join(io))])
    sp.communicate()
    '''

def run_command_with_timeout(cmd, timeout_sec, output_path):
    """Execute `cmd` in a subprocess and enforce timeout `timeout_sec` seconds.
 
    Return subprocess exit code on natural completion of the subprocess.
    Raise an exception if timeout expires before subprocess completes."""
    proc = subprocess.Popen(cmd)
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

def gb2seed(org_name, strains_fname):
    count = 0
    
    if not os.path.exists(os.path.join(SEEDDIR, org_name)):
        os.mkdir(os.path.join(SEEDDIR, org_name))
    
    strains_df = pandas.DataFrame.from_csv(strains_fname, index_col=None)    
    for st_ix, strain_row in strains_df.iterrows(): 
        strain_name = strain_row['strain']
        strain_uid = strain_row['uid']
        dir_name = '_'.join([org_name, strain_name, 'uid%s' % strain_uid])
        print dir_name
        
        if not os.path.exists(os.path.join(SEEDDIR, org_name, dir_name)):
            os.mkdir(os.path.join(SEEDDIR, org_name, dir_name))
        for filename in os.listdir(os.path.join(GBKDIR, dir_name)):
            print filename
            if strain_row['keys'] == re.search('^(.*).gbk$', filename).group(1):
                if not os.path.exists(os.path.join(SEEDDIR, org_name, dir_name, re.search('^(.*).gbk$', filename).group(1))) \
                or len(os.listdir(os.path.join(SEEDDIR, org_name, dir_name, re.search('^(.*).gbk$', filename).group(1)))) == 0:
                    os.chdir('/home/michalsa/Downloads/phiSpyNov11_v2.3/')
                    filename_abs = os.path.join(GBKDIR, dir_name, filename)
                    org_dir = os.path.join(SEEDDIR, org_name, dir_name, re.search('^(.*).gbk$', filename).group(1))
                    subprocess.call(['python', 'genbank_to_seed.py', filename_abs, org_dir])
                else:
                    print "\tseed conversion already exists"
                count+=1
                    
    print "%d seed conversions for %s" % (count, org_name)
 
def runPhiSpy(org_name, strains_fname):
    numOrg = 0
    if org_name == "Escherichia_coli":
        numOrg = 9
    
    if not os.path.exists(os.path.join(PHISPYDIR, org_name)):
        os.mkdir(os.path.join(PHISPYDIR, org_name))
              
    strains_df = pandas.DataFrame.from_csv(strains_fname, index_col=None)    
    for st_ix, strain_row in strains_df.iterrows(): 
        strain_name = strain_row['strain']
        strain_uid = strain_row['uid']
        strain_key = strain_row['keys']
        dir_name = '_'.join([org_name, strain_name, 'uid%s' % strain_uid])
        print dir_name
        if not os.path.exists(os.path.join(PHISPYDIR, org_name, dir_name)):
            os.mkdir(os.path.join(PHISPYDIR, org_name, dir_name))
        os.chdir('/home/michalsa/Downloads/phiSpyNov11_v2.3/') 
        dirname_abs = os.path.join(SEEDDIR, org_name, dir_name, strain_key)
        output_dir = os.path.join(PHISPYDIR, org_name, dir_name, strain_key)
        if not os.path.exists(output_dir):
            subprocess.call(['python', 'phiSpy.py', '-i', dirname_abs, '-o', output_dir, '-t', '%s' % numOrg])
        else:
            print '\tviral data file already exists'

def createPanGenome(org_name, strains_fname, unifiable_fname):  
    pairwise_output_path = os.path.join(OUTPUTDIR, org_name, 'fasta_out_fluidity')
       
    fasta_output_path = os.path.join(OUTPUTDIR, org_name, 'fasta_out_pangenome')
    if not os.path.isdir(fasta_output_path): os.mkdir(fasta_output_path)
    
    strains_df = pandas.DataFrame.from_csv(strains_fname, index_col=None)
    unifiable_df = pandas.DataFrame.from_csv(unifiable_fname)
    
    fluitidy_thresholds = list(unifiable_df.columns.values)
    
    for threshold in fluitidy_thresholds: 
        threshold_str = '%se-%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
        
        '''create directory for threshold, if it doesn't exist'''
        if not os.path.exists(os.path.join(fasta_output_path, threshold_str)):
            os.mkdir(os.path.join(fasta_output_path, threshold_str))
        
        '''create dict for pangenes'''
        pangenes = defaultdict(lambda: defaultdict (list)) 
        
        '''fasta filepath for pangenes'''
        pangenome_filepath = os.path.join(OUTPUTDIR, org_name, 'pangenome_identity%s_fluidity%s.ffn' % (IDENTITYPCNT, threshold_str))
        
        '''create csv summary file, if exists - add its contents to pangenes dict'''
        summary_path = os.path.join(OUTPUTDIR, org_name, 'pangenome_summary_identity%s_fluidity%s.tsv' % (IDENTITYPCNT, threshold_str))
        if os.path.isfile(summary_path) and os.stat(summary_path).st_size!=0:
            summary_fh = open(summary_path, 'r')
            #summary_reader = csv.reader(summary_fh, delimiter='\t')
            for line in summary_fh:
                row = line.split('\t')
                if len(row)==3:
                    pangenes.update({row[0]: {'paralogs': row[1].split(';'),
                                              'orthologs': row[2].split(';')}})
            print '%d records' % len(pangenes.keys())

        '''create directory for unified genomes and pangenome output, if they don't exist'''
        al_types = ['unified', 'pangenome']
        for al_type in al_types:
            if not os.path.exists(os.path.join(fasta_output_path, threshold_str, al_type)):
                os.mkdir(os.path.join(fasta_output_path, threshold_str, al_type))
        
        print threshold        
        for st_ix, strain_row in strains_df.iterrows():
            strain_name = strain_row['strain']
            uid = strain_row['uid']
            q_name = strain_row['keys']
                   
            strain_idx = '_'.join([strain_name, 'uid%s' % uid])
            dir_name = '_'.join([org_name, strain_idx])
            
            if type(q_name)!=str:
                raise ValueError('Incorrect key format for %s: %s. Should be str' % (dir_name, type(q_name)))
            
            q_path = os.path.join(FFNDIR, dir_name, '%s.ffn' % q_name)
            #print q_path
                            
            '''check whether strain is already accounted for in existing unified genomes - default is True, i.e., unified genome should be created for strain'''
            unifiable_flag = True
            for filename in os.listdir(os.path.join(fasta_output_path, threshold_str, 'unified')):
                if re.search(strain_idx, filename):
                    unifiable_flag = False
                    break
            
            if unifiable_flag == True:
                strain_unifiables = []
                if type(unifiable_df[threshold][strain_idx])==str: strain_unifiables = unifiable_df[threshold][strain_idx].split(';')
                
                unifiables_dir = os.path.join(fasta_output_path, threshold_str, 'unified')
                
                if len(strain_unifiables) > 0:
                    unified_filename = os.path.join(unifiables_dir, '%s_%s.ffn' % (strain_idx, '_'.join(strain_unifiables)))
                else:
                    unified_filename = os.path.join(unifiables_dir, '%s.ffn' % strain_idx)
                
                '''paste data from first entry in unified genome'''
                unified_fh = open(unified_filename, 'a')
                q_records = SeqIO.parse(open(q_path, 'r'), 'fasta')
                for q_record in q_records:
                    SeqIO.write(q_record, unified_fh, 'fasta')
                print strain_idx, '%d records' % sum([1 for record in SeqIO.parse(open(unified_filename, 'r'), 'fasta')])
                
                for unifiable_idx in strain_unifiables:
                    unique_count = 0
                    unifiable_name = None
                    unifiable_st, unifiable_uid = unifiable_idx.split('_uid')
                    for ix, rw in strains_df.iterrows():
                        if str(rw['strain']) == str(unifiable_st) and int(rw['uid']) == int(unifiable_uid):
                            unifiable_name = rw['keys']
                            #print rw['strain'], unifiable_st, rw['uid'], unifiable_uid

                    unifiable_fh = open(os.path.join(FFNDIR, '_'.join([org_name, unifiable_idx]), '%s.ffn' % unifiable_name), 'r') 
                    print unifiable_idx, '%d records' % sum([1 for record in SeqIO.parse(unifiable_fh, 'fasta')])

                    '''check whether pairwise alignment exists'''
                    pairwise_res = os.path.join(pairwise_output_path, '%s_vs_%s.txt' % (unifiable_idx, strain_idx))
                    if pairwise_res.split('/')[-1] not in os.listdir(pairwise_output_path) or os.stat(pairwise_res).st_size==0:
                        raise NameError('No such file: %s' % pairwise_res)
                    
                    res_fh = open(pairwise_res, 'r')
                    for line in res_fh:
                        if not re.search('^#', line):
                            al_dict = hit2dict(line)
                            unique_flag = False
                            '''if gene is unique, add it to unified genome'''
                            if al_dict['identity%'] < IDENTITYPCNT or al_dict['al_length'] < ALMINLENGTH*al_dict['q_length']:
                                unique_flag = True
                                '''parse all pairwise comparisons between unifiable and previous unifiables'''
                                if strain_unifiables.index(unifiable_idx)>0:
                                    for unifiable_comp in strain_unifiables[:strain_unifiables.index(unifiable_idx)]:
                                        unifiables_res = os.path.join(pairwise_output_path, '%s_vs_%s.txt' % (unifiable_idx, unifiable_comp))
                                        unifiables_fh = open(unifiables_res, 'r')
                                        for line_comp in unifiables_fh:
                                            if not re.search('^#', line_comp):
                                                '''unifiable is query'''
                                                al_comp_dict = hit2dict(line_comp)
                                                if al_comp_dict['query_id']==al_dict['query_id'] and \
                                                not (al_comp_dict['identity%'] < IDENTITYPCNT or al_comp_dict['al_length'] < ALMINLENGTH*al_comp_dict['q_length']):
                                                    unique_flag = False
                                                    break
                                        if unique_flag==False:
                                            break
                            '''add record to unified genome only if it didn't exist in any of the previous unifiables'''
                            if unique_flag==True:
                                unique_count+=1
                                unifiable_fh.seek(0)
                                unifiable_records = SeqIO.parse(unifiable_fh, 'fasta')
                                for record in unifiable_records:
                                    if record.id == al_dict['query_id']:
                                        SeqIO.write(record, unified_fh, 'fasta')
                                        break
                    print '%d unique' % unique_count
                
                print unified_filename.split('/')[-1], '%d records' % sum([1 for record in SeqIO.parse(open(unified_filename, 'r'), 'fasta')])
     
        '''iterate over unified genomes'''
        for unified_genome in sorted(os.listdir(os.path.join(fasta_output_path, threshold_str, 'unified'))):
            print unified_genome
            
            genes_count, unique_genes_count = 0, 0
            
            genome_filepath = os.path.join(fasta_output_path, threshold_str, 'unified', unified_genome)

            pangenome_alignment = os.path.join(fasta_output_path, threshold_str, 'pangenome', '%s.out' % re.search('(.*)\.ffn', unified_genome).group(1))
            
            if not os.path.isfile(pangenome_filepath) or os.stat(pangenome_filepath).st_size==0:
                SeqIO.write(SeqIO.parse(open(genome_filepath, 'r'), 'fasta'), open(pangenome_filepath, 'w'), 'fasta')
                pangenome_fh = open(pangenome_filepath, 'r')
                pangenome_fasta = SeqIO.parse(pangenome_fh, 'fasta')
                
                for record in pangenome_fasta:
                    genes_count+=1
                    pangenes.update({record.id: {'paralogs': [],
                                                 'orthologs': []}})
            
            else:
                if not os.path.isfile(pangenome_alignment): 
                    run_fasta(genome_filepath, pangenome_filepath, pangenome_alignment)

                if os.stat(pangenome_alignment).st_size==0:
                    run_fasta(genome_filepath, pangenome_filepath, pangenome_alignment, use_stdin=True)
            
            if os.path.isfile(pangenome_alignment) and os.stat(pangenome_alignment).st_size>0:
                al_outfh = open(pangenome_alignment, 'r')
                genome_fh = open(genome_filepath, 'r')
                pangenome_fh = open(pangenome_filepath, 'a')
                ids = []
                for line in al_outfh:
                    genome_fh.seek(0)
                    genome_fasta = SeqIO.parse(genome_fh, 'fasta')
                    
                    if not re.search('^#', line):
                        al_dict = hit2dict(line)
                        '''unified genome is query, pangenome is subject'''
                        if not al_dict['query_id'] in ids: 
                            ids.append(al_dict['query_id'])
                            genes_count+=1
                        else: continue
                        
                        '''if gene is unique, add it to pangenome; else, append it to the corresponding list'''
                        if al_dict['identity%'] < IDENTITYPCNT or al_dict['al_length']<ALMINLENGTH*al_dict['q_length']:
                            unique_genes_count+=1
                            pangenes.update({al_dict['query_id']: {'paralogs': [],
                                                                   'orthologs': []}})
                            #print 'query id: %s' % al_dict['query_id']
                            for record in genome_fasta:
                                if record.id == al_dict['query_id']:
                                    #print 'record id: %s' % record.id
                                    SeqIO.write(record, pangenome_fh, 'fasta')
                                    break
                        
                        else:
                            for pangene_key in pangenes.keys():
                                if al_dict['subject_id']==pangene_key:
                                    if al_dict['query_id'].split('|')[3]==pangene_key.split('|')[3]:
                                        pangenes[pangene_key]['paralogs'].append(al_dict['query_id'])
                                    else:
                                        pangenes[pangene_key]['orthologs'].append(al_dict['query_id'])
                                    break
                    
            print sum([1 for key in pangenes.keys()])
            pangenes_count = sum([1 for record in SeqIO.parse(open(pangenome_filepath, 'r'), 'fasta')])
            print '%d genes, %d unique; pangenome contains %d genes' % (genes_count, unique_genes_count, pangenes_count)
                
            '''write all information stored in pangenes dict to an output file'''
            summary_fh = open(summary_path, 'w')
            #genes_summary = csv.writer(summary_fh, delimiter='\t')
            for pangene in pangenes.keys():
                row = [pangene, ';'.join(pangenes[pangene]['paralogs']), ';'.join(pangenes[pangene]['orthologs'])]
                summary_fh.write('%s\n' % '\t'.join(row))
            
    return

def pangenome_self_compare(org_name):
    dir_path = os.path.join(OUTPUTDIR, org_name)
    list_dir = os.listdir(dir_path)
    
    for filename in list_dir:
        if re.search('^pangenome.*\.ffn$', filename) and not re.search('selfcmp', filename):
            identity, fluidity = re.search('identity(\d+)', filename).group(1), re.search('fluidity(\d+)', filename).group(1)
            pangenome_fname = os.path.join(dir_path, filename)
            self_cmp_fname = os.path.join(dir_path, 'fasta_out_pangenome', fluidity, 'pangenome_self_cmp.out')
            
            print filename, '%d genes' % sum([1 for record in SeqIO.parse(open(pangenome_fname, 'r'), 'fasta')])
            
            '''get pangenome df'''
            pangenes_fname = os.path.join(dir_path, 'pangenome_summary_identity%s_fluidity%s.csv' % (identity, fluidity))
            pangenes_fh = open(pangenes_fname, 'r')
            #pangenes_csv = csv.reader(open(pangenes_fname, 'r'), delimiter='\t')
            pangenes_dict = defaultdict(list)
            
            for line in pangenes_fh:
                row = line.split('\t')
                pangenes_dict.update({row[0]: row[1:]})
            
            '''perform self comparison on pangenome'''
            if not os.path.isfile(self_cmp_fname): 
                sp = subprocess.Popen(['/bin/bash', '-i', '-c', 'fasta36 -b =2 -d 2 -m 8C -m 10 %s %s > %s' % 
                                       (pangenome_fname, pangenome_fname, self_cmp_fname)])
                sp.communicate()
            if os.stat(self_cmp_fname).st_size==0:
                sp = subprocess.Popen(['/bin/bash', '-i', '-c', 'cat %s | fasta36 -b =2 -d 2 -m 8C -m 10 -n @ %s > %s' % 
                                       (pangenome_fname, pangenome_fname, self_cmp_fname)])
                sp.communicate()
            
            '''find all matchable genes'''
            matches_count = 0
            self_cmp_fh = open(self_cmp_fname, 'r')
            for line in self_cmp_fh:
                if not re.search('^#', line):
                    al_dict = hit2dict(line)
                    if al_dict['query_id']!=al_dict['subject_id']:
                        if not (al_dict['identity%']<int(identity) or al_dict['al_length']<ALMINLENGTH*al_dict['q_length']):
                            query, subject = None, None
                            for key in pangenes_dict.keys():
                                if re.search(re.escape(al_dict['query_id']), key):
                                    query = key
                                if re.search(re.escape(al_dict['subject_id']), key):
                                    subject = key
                            if query!=None and subject!=None:
                                matches_count+=1
                                append_vals = [query, ] + pangenes_dict.pop(query)
                                pangenes_dict[subject]+=append_vals
            print '%d matches' % matches_count
            
            '''copy unique entries to new ffn and write updated dict to new csv'''
            updated_pangenome_fname = os.path.join(dir_path, 'pangenome_identity%s_fluidity%s_selfcmp.ffn' % (identity, fluidity))
            updated_pangenome_fh = open(updated_pangenome_fname, 'w')
            
            updated_pangenes_fname = os.path.join(dir_path, 'pangenome_summary_identity%s_fluidity%s_selfcmp.csv' % (identity, fluidity))
            #updated_pangenes_csv = csv.writer(open(updated_pangenes_fname, 'w'), delimiter='\t')
            updated_pangenes_fh = open(updated_pangenes_fname, 'w')
            
            for record in SeqIO.parse(open(pangenome_fname, 'r'), 'fasta'):
                if record.id in pangenes_dict.keys():
                    print record
                    SeqIO.write(record, updated_pangenome_fh, 'fasta')
                    updated_row = [record.id, ] + pangenes_dict[record.id]
                    updated_pangenes_fh.write('%s\n' % '\t'.join(updated_row))
                
            print '%d records after self compare' % (len(pangenes_dict.keys()))
    return

def create_hist(org_name, strains_fname, unifiable_fname):
    strains_fh = open(os.path.join(OUTPUTDIR, org_name, 'strains_plasmids_ncbiFTP.csv'))
    unifiables_fh = open(os.path.join(OUTPUTDIR, org_name, 'unifiable_strains.csv'))
    paralogs_fh = open(os.path.join(OUTPUTDIR, org_name, 'paralogs.tsv'))
    
    strains = [row for row in csv.DictReader(strains_fh)]
    unifiables = [row for row in csv.DictReader(unifiables_fh)]
    
    unifiables_fh.seek(0)
    unifiables_headers = csv.DictReader(unifiables_fh).fieldnames
    thresholds = [header for header in unifiables_headers if re.search('0\.\d+', header)]
    
    paralogs = list(set([el for row in csv.reader(paralogs_fh, delimiter='\t') for el in row]))
    
    results_fh = open(os.path.join(OUTPUTDIR, org_name, 'hist_results.csv'), 'w')
    
    for threshold in thresholds:
        threshold_str = '%se-%d' % (re.search('([1-9])', threshold).group(1), threshold.count('0'))
        print threshold_str
        
        bins = []
        
        unified_genomes_dir = os.path.join(OUTPUTDIR, org_name, 'fasta_out_pangenome', threshold_str, 'unified')
        for unified_genome_fname in os.listdir(unified_genomes_dir):
            bins.append(re.search('(.*)\.ffn', unified_genome_fname).group(1))
        n = len(os.listdir(unified_genomes_dir))
        print n
        
        fname = os.path.join(OUTPUTDIR, org_name, 'pangenome_summary_identity%d_fluidity%s_corrected_afterselfcmp.tsv' % (IDENTITYPCNT, threshold_str))
        fh = open(fname, 'r')
        
        total = 0
        count_pangenes = [0]*n
        count_viral = [0]*n
        
        for line in fh:
            row = line.strip().split('\t')
            if len(row)>0:
                if row[0] not in paralogs and re.search('(NC_\d+)\.\d', row[0].split('|')[3]):
                #if re.search('(NC_\d+)\.\d', row[0].split('|')[3]):
                    viral_flag = False
                    g_refs = [row[0]]
                    g_id = re.search('(NC_\d+)\.\d', row[0].split('|')[3]).group(1)
                    g_ids = [g_id]
                    if len(row)>1:
                        g_refs+=[g for g in row[1].strip().split(';') if g!='']
                        g_ids+=[re.search('(NC_\d+)\.\d', g.split('|')[3]).group(1) for g in g_refs[1:] if re.search('(NC_\d+)\.\d', g.split('|')[3])]
                    
                    if len(g_refs)!=len(g_ids):
                        print g_refs
                        print g_ids
                        return
                    l_ids = []
                    for g_id in g_ids:
                        l_bins = ''
                        for b in bins:
                            if re.search(g_id, b):
                                l_bins = re.findall('NC_\d+', b) 
                        st_name = ''
                        for strain in strains:
                            if strain['keys']==g_id:
                                st_name = '%s_%s_uid%s' % (org_name, strain['strain'], strain['uid'])
                        phispy_output_fname = os.path.join(PHISPYDIR, org_name, st_name, g_id, 'prophage_tbl.txt')
                        if l_bins=='' or not (g_id in l_bins and True in [b in l_ids for b in l_bins]):
                            if g_id not in l_ids:
                                l_ids.append(g_id)
                        
                        '''add search for gene in prophage table generated by phispy'''
                        if os.path.isfile(phispy_output_fname):
                            phispy_output_fh = open(phispy_output_fname, 'r')
                            for line in phispy_output_fh:
                                '''3 - start, 4 - stop, 9 - final_status'''
                                res = line.strip().split('\t')
                                g_arr = g_refs[g_ids.index(g_id)].split('|')
                                pos = re.findall('\d+', g_arr[-1])
                                if res[3] in pos and res[4] in pos:
                                    if res[9]=='1':
                                        viral_flag = True
                                    break        
                                    
                    if len(row)==1 or False not in [l not in paralogs for l in l_ids]:
                        total+=1
                        for i in range(1, n+1):
                            if len(l_ids) in range(i, i+1):
                                count_pangenes[i-1]+=1
                                if viral_flag==True:
                                    count_viral[i-1]+=1
                        '''            
                        if len(l_ids)==n:
                            count_pangenes[-1]+=1
                            if viral_flag==True:
                                    count_viral[-1]+=1                        
                        '''        
        output_line = '%s\n' % '\t'.join([','.join([str(itm) for itm in count_pangenes]), 
                                          ','.join([str(itm) for itm in count_viral]), 
                                          ','. join([str(float(count_viral[i])/count_pangenes[i]) for i in range(n)])])
        print output_line
        results_fh.write(output_line)
    
    return results_fh.name

def plot_hist(y_arr_pangenes=[], y_arr_pangenes_viral=[], y_arr_percentages=[]):
    width = 0.4
    
    labels = ['1e-1', '8e-2', '5e-2', '3e-2', '2e-2']
    colors = ['r', 'g', 'b', 'y', 'c', 'm']

    fig, axarr = plt.subplots(2, 5, sharex='none')
        
    rects_pangenes = []
    rects_percentages = []
    for i in range(5):
        x_arr = range(len(y_arr_pangenes[i]))
        if y_arr_pangenes!=[] and y_arr_pangenes_viral!=[]:
            rects_pangenes.append([axarr[0][i].bar([width/2+j for j in x_arr], y_arr_pangenes[i], width, color=colors[i]),
                          axarr[0][i].bar([width/2+j for j in x_arr], y_arr_pangenes_viral[i], width, color=colors[-1])])
        if y_arr_percentages!=[]:
            rects_percentages.append(axarr[1][i].bar([width/2+j for j in x_arr], y_arr_percentages[i], width, color=colors[i]))
    
    for k in range(2): 
        for i in range(5):
            x_arr = range(len(y_arr_pangenes[i]))
            axarr[k][i].set_xlim(-0.5*width,len(x_arr)+0.5*width)
    
    for i in range(5):
        axarr[0][i].set_ylim(0,max(y_arr_pangenes[i])*1.05)
        #axarr[1][i].set_ylim(0,np.round(max(y_arr_percentages[i])+0.05, 2))
        axarr[0][i].set_yticks(np.arange(0, max(y_arr_pangenes[i]), 200))
        #axarr[1][i].set_yticks(np.arange(0, max(y_arr_pangenes[i]), 0.05))
        for k in range(2): 
            axarr[k][i].set_title('%s\n(n=%d)' % (labels[i], len(y_arr_pangenes[i])))
            axarr[k][i].grid(True)
        axarr[0][i].legend((rects_pangenes[i][0], rects_pangenes[i][1]), ('total', 'viral'))
    
    axarr[0][0].set_ylabel('# of pangenes')
    axarr[1][0].set_ylabel('viral genes ratio')
    
    for k in range(2):
        for j in range(5):
            x_arr = range(len(y_arr_pangenes[j]))
            axarr[k][j].set_xticks([i+width for i in x_arr])
            axarr[k][j].set_xticklabels(tuple([str(i+1) if (i==0 or (i+1)%5==0 or (i==len(x_arr)-1 and (i+1)%5>1)) else '' for i in x_arr]))
    
    def autolabel(rects, i, k):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            if height>1: axarr[k][i].text(rect.get_x()+rect.get_width()/2., 1.02*height, height, ha='center', va='bottom')
            elif height>0 and height<1: axarr[k][i].text(rect.get_x()+rect.get_width()/2., 1.02*height, '%.2f'%height, ha='center', va='bottom') 
    '''
    for i in range(5): 
        for j in range(2):
            autolabel(rects_pangenes[i][j], i, 0)
        autolabel(rects_percentages[i], i, 1)
    '''
    plt.show()
    return

def calculate_enc(org_name, prog_name):
    strains_fh = open(os.path.join(OUTPUTDIR, org_name, 'strains_plasmids_ncbiFTP.csv'))
    unifiables_fh = open(os.path.join(OUTPUTDIR, org_name, 'unifiable_strains.csv'))
    
    
    strains = [row for row in csv.DictReader(strains_fh)]
    #unifiables = [row for row in csv.DictReader(unifiables_fh)]
    
    for strain_row in strains:
        print '_'.join([org_name, strain_row['strain'], 'uid%s' % strain_row['uid']])
         
        strain_fname = os.path.join(FFNDIR, 
                                    '_'.join([org_name, strain_row['strain'], 'uid%s' % strain_row['uid']]),
                                    '%s.ffn' % strain_row['keys'])
        n_genes = sum([1 for record in SeqIO.parse(open(strain_fname, 'r'), 'fasta')])
        
        encprime_output_dir = os.path.join(OUTPUTDIR, org_name, 'encprime', '_'.join([org_name, strain_row['strain'], 'uid%s' % strain_row['uid']])) 
        
        if not os.path.isdir(encprime_output_dir):
            os.mkdir(encprime_output_dir)
        
        os.chdir(encprime_output_dir)
        if prog_name=='SeqCount_mod':
            for flag in ['-c', '-n']:
                opts = '%s %s %s %d' % (os.path.join(ENCPATH, prog_name), flag, strain_fname, n_genes)
                subprocess.call(['/bin/bash', '-i', '-c', opts])
        elif prog_name=='ENCprime':
            opts = '%s %s %s 11 %s 1 -q' % (os.path.join(ENCPATH, prog_name), 
                                         '%s.ffn.codcnt' % strain_row['keys'], 
                                         '%s.ffn.acgtfreq' % strain_row['keys'], 
                                         '%s.ffn.results' % strain_row['keys'])
            subprocess.call(['/bin/bash', '-i', '-c', opts])
        
    return

def encprime_per_pangene(org_name):
    strains_fh = open(os.path.join(OUTPUTDIR, org_name, 'strains_plasmids_ncbiFTP.csv'))
    unifiables_fh = open(os.path.join(OUTPUTDIR, org_name, 'unifiable_strains.csv'))
    paralogs_fh = open(os.path.join(OUTPUTDIR, org_name, 'paralogs.tsv'))
    
    strains = [row for row in csv.DictReader(strains_fh)]
    unifiables = [row for row in csv.DictReader(unifiables_fh)]
    unifiables_fh.seek(0)
    unifiables_headers = csv.DictReader(unifiables_fh).fieldnames
    thresholds = [header for header in unifiables_headers if re.search('0\.\d+', header)]
    
    paralogs = list(set([el for row in csv.reader(paralogs_fh, delimiter='\t') for el in row]))
    
    for threshold in thresholds:
        threshold_str = '%se-%d' % (re.search('([1-9])', threshold).group(1), threshold.count('0'))
        print threshold_str
        
        output_fname = os.path.join(OUTPUTDIR, org_name, 'pangenome_summary_identity%d_fluidity%s_corrected_stats.tsv' % (IDENTITYPCNT, threshold_str))
        output_fh = open(output_fname, 'a')
        
        bins = []
        
        unified_genomes_dir = os.path.join(OUTPUTDIR, org_name, 'fasta_out_pangenome', threshold_str, 'unified')
        for unified_genome_fname in os.listdir(unified_genomes_dir):
            bins.append(re.search('(.*)\.ffn', unified_genome_fname).group(1))
        n = len(os.listdir(unified_genomes_dir))
        print n

        enc = [[]]*n
        enc_viral = [[]]*n
        encprime = [[]]*n
        encprime_viral = [[]]*n
                
        pangenome_fname = os.path.join(OUTPUTDIR, org_name, 'pangenome_summary_identity%d_fluidity%s_corrected_afterselfcmp.tsv' % (IDENTITYPCNT, threshold_str))
        pangenome_fh = open(pangenome_fname, 'r')
        
        for line in pangenome_fh:
            row = line.strip().split('\t')
            if len(row)>0:
                if row[0] not in paralogs and re.search('(NC_\d+)\.\d', row[0].split('|')[3]):
                    viral_flag = False
                    enc_i, encprime_i = 0, 0
                    g_refs = [row[0]]
                    g_id = re.search('(NC_\d+)\.\d', row[0].split('|')[3]).group(1)
                    g_ids = [g_id]
                    if len(row)>1:
                        g_refs+=[g for g in row[1].strip().split(';') if g!='']
                        g_ids+=[re.search('(NC_\d+)\.\d', g.split('|')[3]).group(1) for g in g_refs[1:] if re.search('(NC_\d+)\.\d', g.split('|')[3])]
                    l_ids = []
                    for g_id in g_ids:
                        l_bins = ''
                        for b in bins:
                            if re.search(g_id, b):
                                l_bins = re.findall('NC_\d+', b) 
                        st_name = ''
                        for strain in strains:
                            if strain['keys']==g_id:
                                st_name = '%s_%s_uid%s' % (org_name, strain['strain'], strain['uid'])
                        phispy_output_fname = os.path.join(PHISPYDIR, org_name, st_name, g_id, 'prophage_tbl.txt')
                        if l_bins=='' or not (g_id in l_bins and True in [b in l_ids for b in l_bins]):
                            if g_id not in l_ids:
                                l_ids.append(g_id)
                        
                        '''add search for gene in prophage table generated by phispy'''
                        if os.path.isfile(phispy_output_fname):
                            phispy_output_fh = open(phispy_output_fname, 'r')
                            for line in phispy_output_fh:
                                '''3 - start, 4 - stop, 9 - final_status'''
                                res = line.strip().split('\t')
                                g_arr = g_refs[g_ids.index(g_id)].split('|')
                                pos = re.findall('\d+', g_arr[-1])
                                if res[3] in pos and res[4] in pos:
                                    if res[9]=='1':
                                        viral_flag = True
                                    break
                                
                        '''calculate average ENC, ENCprime'''
                        encprime_output_fname = os.path.join(OUTPUTDIR, org_name, 'encprime', st_name, 
                                                             '%s.ffn.results' % re.search('\|(NC_\d+)\.\d\|', g_refs[g_ids.index(g_id)]).group(1))
                        
                        for ln in open(encprime_output_fname, 'r'):
                            rw = ln.strip()
                            if re.search(g_refs[g_ids.index(g_id)].split('|')[-1], rw):
                                nums = re.findall('\d+\.\d+', rw)
                                enc_i+=float(nums[1])
                                encprime_i+=float(nums[2])
                    
                    enc_v = enc_i/len(g_refs)
                    encprime_v = encprime_i/len(g_refs)
                    
                    output_fh.write('%s\n' % '\t'.join([g_refs[0], str(len(g_refs)), str(len(l_ids)), str(viral_flag), str(enc_v), str(encprime_v)]))
    return

def create_enc_hist(org_name):
    unifiables_fh = open(os.path.join(OUTPUTDIR, org_name, 'unifiable_strains.csv'))
    unifiables_headers = csv.DictReader(unifiables_fh).fieldnames
    thresholds = [header for header in unifiables_headers if re.search('0\.\d+', header)]
    
    fig, axarr = plt.subplots(2, len(thresholds), sharex='none')
    width = 0.4
    
    vals_enc_encp = [[] for i in range(len(thresholds))]
    
    bins_enc_encp = [[] for i in range(len(thresholds))]
    bins_viral_enc_encp = [[] for i in range(len(thresholds))]
    
    for threshold in thresholds:
        threshold_str = '%se-%d' % (re.search('([1-9])', threshold).group(1), threshold.count('0'))
        print threshold_str

        unified_genomes_dir = os.path.join(OUTPUTDIR, org_name, 'fasta_out_pangenome', threshold_str, 'unified')
        n = len(os.listdir(unified_genomes_dir))
        print n
        
        bins_enc_encp[thresholds.index(threshold)] = [[[],[]] for i in range(n)]
        bins_viral_enc_encp[thresholds.index(threshold)] = [[[],[]] for i in range(n)]
        
        vals_enc_encp[thresholds.index(threshold)] = [[] for i in range(4)]
        
        output_fname = os.path.join(OUTPUTDIR, org_name, 'pangenome_summary_identity%d_fluidity%s_corrected_stats.tsv' % (IDENTITYPCNT, threshold_str))
        output_fh = open(output_fname, 'r')
        
        for line in output_fh:
            '''0- pangene name; 1- n_orthologs; 2-n_genomes; 3-viral_flag; 4-enc; 5-encprime'''
            row = line.strip().split('\t')
            if row[3]=='True':
                bins_viral_enc_encp[thresholds.index(threshold)][int(row[2])-1][0].append(float(row[4]))
                bins_viral_enc_encp[thresholds.index(threshold)][int(row[2])-1][1].append(float(row[5]))
            else:
                bins_enc_encp[thresholds.index(threshold)][int(row[2])-1][0].append(float(row[4]))
                bins_enc_encp[thresholds.index(threshold)][int(row[2])-1][1].append(float(row[5]))
        
        
        
        for i in range(n):
            vals_enc_encp[thresholds.index(threshold)][0].append(np.average(bins_enc_encp[thresholds.index(threshold)][i][0]))
            vals_enc_encp[thresholds.index(threshold)][1].append(np.average(bins_viral_enc_encp[thresholds.index(threshold)][i][0]))
            vals_enc_encp[thresholds.index(threshold)][2].append(np.average(bins_enc_encp[thresholds.index(threshold)][i][1]))
            vals_enc_encp[thresholds.index(threshold)][3].append(np.average(bins_viral_enc_encp[thresholds.index(threshold)][i][1]))
            '''
            print np.average(bins_enc_encp[i][0]), np.average(bins_viral_enc_encp[i][0]), len(bins_enc_encp[i][0])*len(bins_viral_enc_encp[i][0]), \
            stats.mannwhitneyu(bins_enc_encp[i][0], bins_viral_enc_encp[i][0]), \
            '\t', np.average(bins_enc_encp[i][1]), np.average(bins_viral_enc_encp[i][1]), len(bins_enc_encp[i][1])*len(bins_viral_enc_encp[i][1]), \
            stats.mannwhitneyu(bins_enc_encp[i][1], bins_viral_enc_encp[i][1])
            '''
    
    def autolabel(rects, t, k, u):
        # label significant values
        for rect in rects[0]:
            height = max(rect.get_height(), rects[1][rects[0].index(rect)].get_height())
            if u[rects[0].index(rect)][1]<=0.05:
                print t, k, rect.get_height(), rects[1][rects[0].index(rect)].get_height()
                axarr[k][t].text(rect.get_x()+rect.get_width(), 1.005*height, '*', ha='center', va='bottom')
    
    colors = ('r', 'g', 'b', 'y', 'c', 'm')
    series = ('non-viral', 'viral')
    for i in range(len(thresholds)):
        rects = [[], []]
        x_arr = range(len(vals_enc_encp[i][0]))
        axarr[0][i].set_title('%s\n(n=%d)' % (thresholds[i], len(vals_enc_encp[i][0])))
        for j in range(4):
            ind = np.arange(len(vals_enc_encp[i][j]))
            w = width*(j%2)
            c = colors[-1]
            if j%2==0:
                c = colors[i]
            rects[j//2].append(axarr[j//2][i].bar(ind+w, tuple([val for val in vals_enc_encp[i][j]]), width, color=c))
        
        for k in range(len(series)):
            axarr[k][i].set_xlim(-0.5*width,len(x_arr)+0.5*width)
            axarr[k][i].set_ylim(min([rect.get_height() for rect in rects[k][0]+rects[k][1]])*0.98, 
                                 max([rect.get_height() for rect in rects[k][0]+rects[k][1]])*1.02)
            axarr[k][i].grid(True)
            axarr[k][i].legend((rects[k][0], rects[k][1]), series, loc=3)
            u = [stats.mannwhitneyu(bins_enc_encp[i][m][k], bins_viral_enc_encp[i][m][k]) for m in range(len(x_arr))]
            autolabel(rects[k], i, k, u)
            axarr[k][i].set_xticks([x+width for x in x_arr])
            axarr[k][i].set_xticklabels(tuple([str(x+1) if (x==0 or (x+1)%5==0 or (x==len(x_arr)-1 and (x+1)%5>1)) else '' for x in x_arr]))
    
    axarr[0][0].set_ylabel('ENC')
    axarr[1][0].set_ylabel('ENCprime')
    
    plt.show()
    return

def calculate_fop(org_name):
    strains_fh = open(os.path.join(OUTPUTDIR, org_name, 'strains_plasmids_ncbiFTP.csv'))       
    strains = [row for row in csv.DictReader(strains_fh)]
    
    for strain_row in strains:
        print '_'.join([org_name, strain_row['strain'], 'uid%s' % strain_row['uid']])
         
        strain_fname = os.path.join(FFNDIR, 
                                    '_'.join([org_name, strain_row['strain'], 'uid%s' % strain_row['uid']]),
                                    '%s.ffn' % strain_row['keys'])
        
        fop_output_dir = os.path.join(OUTPUTDIR, org_name, 'fop', '_'.join([org_name, strain_row['strain'], 'uid%s' % strain_row['uid']])) 
        
        output_fname = os.path.join(fop_output_dir, '%s.fop.out' % strain_row['keys'])
        
        if not os.path.isdir(fop_output_dir):
            os.mkdir(fop_output_dir)
        
        output_fh = open(output_fname, 'a')
        
        for gene in SeqIO.parse(open(strain_fname, 'r'), 'fasta'):    
            print gene.name, len(gene.seq)
            count_opt, count_tot = 0, 0
            for i in range(0, len(gene.seq), 3):
                #print gene.seq[i:i+3], [key for key in AADICT.keys() if (str(gene.seq[i:i+3]) in AADICT[key] or str(gene.seq[i:i+3]).swapcase() in AADICT[key])]
                if len([key for key in AADICT.keys() if (str(gene.seq[i:i+3]) in AADICT[key] or str(gene.seq[i:i+3]).swapcase() in AADICT[key])])>0: 
                    count_tot+=1
                    '''print [key_opt for key_opt in OPTDICT[org_name].keys() if (str(gene.seq[i:i+3]) in OPTDICT[org_name][key_opt] or \
                                                                               str(gene.seq[i:i+3]).swapcase() in OPTDICT[org_name][key_opt])]'''
                    if len([key_opt for key_opt in OPTDICT[org_name].keys() if (str(gene.seq[i:i+3]) in OPTDICT[org_name][key_opt] or \
                                                                                str(gene.seq[i:i+3]).swapcase() in OPTDICT[org_name][key_opt])])>0:
                        count_opt+=1
            
            '''0- gene_name; 1-gene_length; 2-n_optimal_codons; 3- n_optimal_nonoptimal_codons; 5-fop'''
            output_fh.write('%s\n' % '\t'.join([str(itm) for itm in [gene.name, len(gene.seq), count_opt, count_tot, float(count_opt)/count_tot]]))
            #print count_opt, count_tot, float(count_opt)/count_tot
    return

def plot_fop(org_name):
    strains_fh = open(os.path.join(OUTPUTDIR, org_name, 'strains_plasmids_ncbiFTP.csv'))
    unifiables_fh = open(os.path.join(OUTPUTDIR, org_name, 'unifiable_strains.csv'))
    paralogs_fh = open(os.path.join(OUTPUTDIR, org_name, 'paralogs.tsv'))
    
    strains = [row for row in csv.DictReader(strains_fh)]
    unifiables = [row for row in csv.DictReader(unifiables_fh)]
    
    unifiables_fh.seek(0)
    unifiables_headers = csv.DictReader(unifiables_fh).fieldnames
    thresholds = [header for header in unifiables_headers if re.search('0\.\d+', header)]
    
    paralogs = list(set([el for row in csv.reader(paralogs_fh, delimiter='\t') for el in row]))
    
    for threshold in thresholds:
        threshold_str = '%se-%d' % (re.search('([1-9])', threshold).group(1), threshold.count('0'))
        print threshold_str
        
        bins = []
        
        unified_genomes_dir = os.path.join(OUTPUTDIR, org_name, 'fasta_out_pangenome', threshold_str, 'unified')
        for unified_genome_fname in os.listdir(unified_genomes_dir):
            bins.append(re.search('(.*)\.ffn', unified_genome_fname).group(1))
        n = len(os.listdir(unified_genomes_dir))
        print n
        
        pangenome_fname = os.path.join(OUTPUTDIR, org_name, 'pangenome_summary_identity%d_fluidity%s_afterselfcmp.tsv' % (IDENTITYPCNT, threshold_str))
        pangenome_fh = open(pangenome_fname, 'r')
        
        stats_fname = os.path.join(OUTPUTDIR, org_name, 'pangenome_summary_identity%d_fluidity%s_stats.tsv' % (IDENTITYPCNT, threshold_str))
        stats_fh = open(stats_fname, 'r') 
        
        for line_p in pangenome_fh:
            row_p = line_p.strip().split('\t')
            if len(row_p)>0:
                if row_p[0] not in paralogs and re.search('(NC_\d+)\.\d', row_p[0].split('|')[3]):
                    line_s = next(stats_fh)
                    row_s = line_s.strip().split('\t')
                    g_refs = [row_p[0]]
                    g_id = re.search('(NC_\d+)\.\d', row_p[0].split('|')[3]).group(1)
                    g_ids = [g_id]
                    if len(row_p)>1:
                        g_refs+=[g for g in row_p[1].strip().split(';') if g!='']
                        g_ids+=[re.search('(NC_\d+)\.\d', g.split('|')[3]).group(1) for g in g_refs[1:] if re.search('(NC_\d+)\.\d', g.split('|')[3])]
                    
                    fop_arr = []     
                    #print len(g_ids)               
                    for g_ix, g_id in enumerate(g_ids):
                        g_ref = g_refs[g_ix]
                        st_name = ''
                        for strain in strains:
                            if strain['keys']==g_id:
                                st_name = '_'.join([org_name, strain['strain'], 'uid%s' % strain['uid']])
                                fop_output_dir = os.path.join(OUTPUTDIR, org_name, 'fop', st_name) 
                                fop_fname = os.path.join(fop_output_dir, '%s.fop.out' % g_id)
                                fop_fh = open(fop_fname, 'r')
                                for line in fop_fh:
                                    fop_row = line.strip().split('\t')
                                    if fop_row[0]==g_ref:
                                        fop_arr.append(float(fop_row[4]))
                    
                    print int(row_s[2]), np.average(fop_arr), np.std(fop_arr), row_s[0]==row_p[0]
                    if int(row_s[2])==1 and np.std(fop_arr)>0:
                        print row_p
                        return
                    if row_s[0]!=row_p[0]:
                        print row_s[0], row_p[0]
                        return
        return

        
def get_highly_expressed_genes(org_name):
    strains_fh = open(os.path.join(OUTPUTDIR, org_name, 'strains_plasmids_ncbiFTP.csv'))
    unifiables_fh = open(os.path.join(OUTPUTDIR, org_name, 'unifiable_strains.csv'))
    
    strains = [row for row in csv.DictReader(strains_fh)]
    unifiables = [row for row in csv.DictReader(unifiables_fh)]
    
    unifiables_fh.seek(0)
    unifiables_headers = csv.DictReader(unifiables_fh).fieldnames
    thresholds = [header for header in unifiables_headers if re.search('0\.\d+', header)]
    
    for threshold in thresholds:
        threshold_str = '%se-%d' % (re.search('([1-9])', threshold).group(1), threshold.count('0'))
        print threshold_str
        
        pangenome_fname = os.path.join(OUTPUTDIR, org_name, 'pangenome_summary_identity%d_fluidity%s_afterselfcmp.tsv' % (IDENTITYPCNT, threshold_str))
        pangenome_fh = open(pangenome_fname, 'r')
        
        pangenome_heg_fname = os.path.join(OUTPUTDIR, org_name, 'pangenome_identity%d_fluidity%s_heg.ffn' % (IDENTITYPCNT, threshold_str))
        if not os.path.isfile(pangenome_heg_fname):
            pangenome_heg_fh = open(pangenome_heg_fname, 'a')
            heg_records = []
            found_rows = []
            for strain in strains:
                print strain['strain']
                gbk_fname = os.path.join(GBKDIR, '%s_%s_uid%s' % (org_name, strain['strain'], strain['uid']), '%s.gbk' % strain['keys'])
                gbk_fh = SeqIO.parse(open(gbk_fname, 'r'), 'genbank')
                for record in gbk_fh:
                    for feature in record.features:
                        for k in feature.qualifiers.keys():
                            for term in ['ribosomal', 'elongation factor']:
                                if [re.search(term, feature.qualifiers[k][j]) for j in range(len(feature.qualifiers[k]))]!=[None]*len(feature.qualifiers[k]):
                                    pangenome_fh.seek(0)
                                    found_flag = False
                                    for line in pangenome_fh:
                                        row = line.strip().split('\t')
                                        if len(row)>0 and row not in found_rows and len(row[0].split('|'))==4:
                                            g_refs = [row[0]]
                                            #print row[0].split('|')
                                            g_id = re.search('(NC_\d+)\.\d', row[0].split('|')[3]).group(1)
                                            g_ids = [g_id]
                                            if len(row)>1:
                                                g_refs+=[g for g in row[1].strip().split(';') if g!='']
                                                g_ids+=[re.search('(NC_\d+)\.\d', g.split('|')[3]).group(1) for g in g_refs[1:] if re.search('(NC_\d+)\.\d', g.split('|')[3])]
                                                
                                            for g_id in g_ids:
                                                g_arr = g_refs[g_ids.index(g_id)].split('|')
                                                pos = [int(p) for p in re.findall('\d+', g_arr[-1])]
                                                #print pos
                                                
                                                if min(pos) in range(feature.location.start-1, feature.location.start+2) \
                                                and max(pos) in range(feature.location.end-1, feature.location.end+2):
                                                    heg_record = SeqRecord.SeqRecord(feature.extract(record.seq), 
                                                                                     id=g_refs[g_ids.index(g_id)],
                                                                                     name=g_refs[g_ids.index(g_id)], 
                                                                                     description=g_refs[g_ids.index(g_id)])
                                                    if heg_record.id not in [rec.id for rec in heg_records]:
                                                        print feature.location, '\n', feature.qualifiers
                                                        print pos, '\n', feature.extract(record.seq)
                                                        heg_records.append(heg_record)
                                                        found_rows.append(row)
                                                        SeqIO.write(heg_record, pangenome_heg_fh, 'fasta')
                                                    found_flag = True
                                                    break
                                        if found_flag==True:
                                            break
                                        
            print sum([1 for record in heg_records])
    return        

def get_annotations(directory):
    org_name = directory.split('/')[-1]
    csv.field_size_limit(sys.maxsize)
    #input_csv = csv.reader(open(os.path.join(directory, 'temp', 'genes_summary_withViralFlag.txt'), 'r'))
    #summary_df = pandas.DataFrame.from_csv(os.path.join(directory, 'temp', 'genes_summary_withViralFlag.txt'), index_col=None)
    viral_df = pandas.DataFrame.from_csv(os.path.join(directory, 'temp', 'genes_summary_hist_cutoff90_noplasmids_nophages.txt'), index_col=None)

    core_output_csv = csv.writer(open(os.path.join(directory, 'temp', 'core_genes_annotation_cutoff90_noplasmids_nophages.txt'), 'w'))
    acquired_output_csv = csv.writer(open(os.path.join(directory, 'temp', 'acquired_genes_annotation_cutoff90_noplasmids_nophages.txt'), 'w'))
    
    for idx, pcnt in enumerate(viral_df['%unique_genomes']):
        output_flag = False
        if viral_df['viral_flag'][idx]==1:
            if pcnt < 0.11:
                print viral_df['id'][idx], viral_df['name'][idx], pcnt
                output_csv = acquired_output_csv
                output_flag = True
            elif pcnt > 0.89:
                print viral_df['id'][idx], viral_df['name'][idx], pcnt
                output_csv = core_output_csv
                output_flag = True
        
            if output_flag==True:
                output_csv.writerow([viral_df['id'][idx], viral_df['name'][idx], '{:.2f}'.format(pcnt)])
                pangene_id = viral_df['id'][idx]
                pangene_name = viral_df['name'][idx]
                
                filename = None
                filepath = os.path.join(directory, '%s_%s_refGenome.gb' % (org_name, pangene_id))
                if os.path.isfile(filepath):
                    filename = filepath
                
                else:
                    if re.search('[A-Z]{1,2}\w+\.\d', pangene_id):
                        filepath = os.path.join(directory, '%s_%s_%s_refGenome.gb' % (org_name, re.search('([A-Z]{1,2})\w+\.\d', pangene_id).group(1), 
                                                                                      re.search('[A-Z]{1,2}(\w+\.\d)', pangene_id).group(1)))
                    if os.path.isfile(filepath):
                        filename = filepath
                    else:
                        print 'no such file: %s' % filepath
                        return 
                
                pangene_gbfile = SeqIO.parse(open(filename, 'r'), 'genbank')
                for record in pangene_gbfile:
                    #print record.features[3064]
                    #return
                    counter = 0
                    for feature in record.features:
                        if feature.type=='CDS':
                            counter+=1
                            if ('gene' in feature.qualifiers.keys() and feature.qualifiers['gene'][0]==re.search('(\w+\d?)_\d+', pangene_name).group(1)) \
                             or counter==int(re.search('\w+\d?_(\d+)', pangene_name).group(1)):
                                #print ';'.join(feature.qualifiers['gene']), '\t', ';'.join(feature.qualifiers['product']), '\t', ';'.join(feature.qualifiers['note'])
                                output = [';'.join(feature.qualifiers[key]) for key in ['gene', 'product', 'note'] if key in feature.qualifiers.keys()]
                                if output:
                                    print '\t'.join([';'.join(feature.qualifiers[key]) for key in ['gene', 'product'] if key in feature.qualifiers.keys()])
                                    output_csv.writerow(output)
                                else:
                                    output_csv.writerow(feature)
                                    return
        
if __name__=="__main__":
    org_dirs = {'Escherichia_coli': []}
    org_folders = sorted([name for name in os.listdir(GBKDIR) if os.path.isdir(os.path.join(GBKDIR, name))])
    for query_org in org_dirs.keys():
        print query_org
        strains_csv = os.path.join(OUTPUTDIR, query_org, 'strains_plasmids_ncbiFTP.csv')
        unifiable_csv = os.path.join(OUTPUTDIR, query_org, 'unifiable_strains.csv')
        
        for org_folder in sorted(org_folders):
            if re.search(query_org, org_folder):
                org_dirs[query_org].append(org_folder)
        '''
        gb2seed(query_org, strains_csv)
        runPhiSpy(query_org, strains_csv)
        '''
        #createPanGenome(query_org, strains_csv, unifiable_csv)
        
        #pangenome_self_compare(query_org)
        '''
        hist_specs = os.path.join(OUTPUTDIR, query_org, 'hist_results.csv')
        if not os.path.isfile(hist_specs):
            hist_specs = create_hist(query_org, strains_csv, unifiable_csv)
            
        pangenes, pangenes_viral, percentages = [], [], []
        for line in open(hist_specs, 'r'):
            row = line.strip().split('\t')
            pangenes.append([int(el) for el in row[0].split(',')])
            pangenes_viral.append([int(el) for el in row[1].split(',')])
            percentages.append([float(el) for el in row[2].split(',')])
        
        print len(pangenes), pangenes
        print len(pangenes_viral), pangenes_viral
        print len(percentages), percentages
        '''
        #plot_hist(pangenes, pangenes_viral, percentages)
        
        #calculate_enc(query_org, 'SeqCount_mod')
        #calculate_enc(query_org, 'ENCprime')
        #encprime_per_pangene(query_org)
        create_enc_hist(query_org)
        #calculate_fop(query_org)
        #plot_fop(query_org)
        