from Bio import SeqIO
import pandas
import numpy
import os, re
import sqlalchemy
import matplotlib
from matplotlib import pyplot as plt

import pangenome
from numpy import NaN

GBKDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/all_gbk/'
FFNDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/all_ffn/'
OUTPUTDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/output/'

ORG_NAMES = ['Bacillus_cereus', 'Campylobacter_jejuni', 'Clostridium_botulinum', 'Corynebacterium_diphtheriae',
             'Corynebacterium_pseudotuberculosis', 'Escherichia_coli-Shigella', 'Francisella_tularensis', 
             'Haemophilus_influenzae', 'Listeria_monocytogenes', 'Mycobacterium_tuberculosis', 'Neisseria_meningitidis',
             'Pseudomonas_aeruginosa', 'Streptococcus_pneumoniae', 'Yersinia_pestis']

CLONAL_ORGS = ['Corynebacterium_pseudotuberculosis', 'Francisella_tularensis', 
               'Mycobacterium_tuberculosis', 'Yersinia_pestis']

def populate_strains_table(org_name, dir_names):
    if not os.path.isdir(os.path.join(OUTPUTDIR, org_name)): os.mkdir(os.path.join(OUTPUTDIR, org_name))
    
    db_fname = os.path.join(OUTPUTDIR, org_name, 'output_datatables_%s.db' % org_name)
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    st_pl = sqlalchemy.Table('strains_plasmids', sqlalchemy.MetaData(engine), autoload=True)
    
    for dir_name in dir_names:
        #print dir_name
        group_names = dir_name.split('-')
        for group_name in group_names:
            print group_name
            strain_dir = os.path.join(GBKDIR, group_name)
            genome_files = [name for name in os.listdir(strain_dir) if re.search('.gbk$', name)]
            print '%d gbk files in directory' % len(genome_files)
            for fname in sorted(genome_files):
                fh = open(os.path.join(strain_dir, fname))
                records = SeqIO.parse(fh, 'genbank')            
                for record in records:
                    #print record.description
                    rec_key = fname.replace('.gbk', '')
                    if not re.search('plasmid', record.description, re.IGNORECASE):   
                        conn = engine.connect()
                        with conn.begin() as trans:
                            print rec_key, group_name
                            try:
                                conn.execute('INSERT INTO strains_plasmids VALUES (?, ?, ?, ?)', (rec_key, group_name, None, None))
                                trans.commit()
                            except:
                                trans.rollback()
                                #print traceback.format_exc()
                        conn.close()
                    else:
                        for feature in record.features:
                            #print record.description
                            if feature.type=='source':
                                pl_name = feature.qualifiers['plasmid'][0]
                                #st_row = engine.execute('SELECT (plasmid_ids, plasmid_names) FROM strains_plasmids WHERE strain_name_uid=(?)', (group_name)).fetchone()
                                st_row = sqlalchemy.select([st_pl.c.plasmid_ids, st_pl.c.plasmid_names]).where(st_pl.c.strain_name_uid==group_name).execute().fetchone()
                                #print st_row
                                conn = engine.connect()
                                with conn.begin() as trans:
                                    try:
                                        if st_row==(None,None):
                                            conn.execute('UPDATE strains_plasmids SET plasmid_ids=(?), plasmid_names=(?) WHERE strain_name_uid=(?)', 
                                                          (rec_key, pl_name, group_name))   
                                        else:
                                            conn.execute('UPDATE strains_plasmids SET plasmid_ids=(?), plasmid_names=(?) WHERE strain_name_uid=(?)', 
                                                          (';'.join([st_row[0], rec_key]), ';'.join([st_row[1], pl_name]), group_name))
                                        trans.commit()
                                    except:
                                        trans.rollback()
                                        #print traceback.format_exc()
                                conn.close()
                                    
    print 'done, table now contains %d strain entries' % engine.execute('SELECT Count(*) FROM strains_plasmids').fetchone()[0]
    return

def get_16s_pos(gbk_filename):
    with open(gbk_filename, 'r') as gbk_fh:
        gbk_records = SeqIO.parse(gbk_fh, 'genbank')
        for gbk_record in gbk_records:
            locations = []
            for feature in gbk_record.features:
                if feature.type=='rRNA':
                    if 'product' in feature.qualifiers.keys() and re.search('16s', feature.qualifiers['product'][0], re.IGNORECASE):
                        locations.append(feature.location)
    if len(locations)>0: 
        return locations
    else:
        return None

def calculate_fluidity(fasta_output_file, q_path, s_path, total_count):           
    if not os.path.isfile(fasta_output_file):
        pangenome.run_fasta(q_path, s_path, fasta_output_file, 600)
    if os.stat(fasta_output_file).st_size==0:
        pangenome.run_fasta(q_path, s_path, fasta_output_file, 600, use_stdin=True)
    if os.stat(fasta_output_file).st_size==0:
        print 'could not process query'
        return
           
    unique_count = 0
    
    with open(fasta_output_file, 'r') as fasta_fh:
    # result indices: 0-query id, 1-subject id, 2-% identity, 3-alignment length, 4-mismatches, 5-gap opens, 
    #                 6-q. start, 7-q. end, 8-s. start, 9-s. end, 10-evalue, 11-bit score
        while True:
            try:
                line = next(fasta_fh)
                hitarr = []
                while not re.search('^#', line):
                    hit = pangenome.hit2obj(line)
                    if hit:
                        unique_flag = pangenome.isunique(hit)
                        hitarr.append(unique_flag)
                    line = next(fasta_fh)
                if hitarr!=[] and False not in hitarr:
                    unique_count+=1
                if re.match('# FASTA processed (\d+) queries', line):
                    matches_count = int(re.match('# FASTA processed (\d+) queries', line).group(1))
                    unique_count += total_count - matches_count
            except StopIteration:
                break
        
    return unique_count 

def compare_strains(org_name):          
    db_fname = os.path.join(OUTPUTDIR, org_name, 'output_datatables_%s.db' % org_name)
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
      
    fasta_output_path = os.path.join(OUTPUTDIR, org_name, 'fasta_out_fluidity')
    if not os.path.isdir(fasta_output_path): os.mkdir(fasta_output_path)
    
    strains_df = pandas.read_sql_table('strains_plasmids', engine, index_col='strain_id')
    counts_df = pandas.DataFrame(index=strains_df.index.values, columns=strains_df.index.values)

    os.chdir(FFNDIR)
    
    for st_ix, strain_row in strains_df.iterrows():
        q_path = os.path.join(strain_row['strain_name_uid'], '%s.ffn' % st_ix)
        if os.path.isfile(q_path): 
            q_records = SeqIO.parse(open(q_path, 'r'), 'fasta')
            total_count = sum([1 for q_record in q_records])
                   
            for comp_st_ix, comp_strain_row in strains_df.iterrows():
                if comp_st_ix!=st_ix:
                    s_path = os.path.join(comp_strain_row['strain_name_uid'], '%s.ffn' % comp_st_ix)
                    output_fname = os.path.join(fasta_output_path, '%s_vs_%s.out' % (st_ix, comp_st_ix))
                    #unique_count = compare_pool.apply_async(calculate_fluidity, (output_fname, q_path, s_path, total_count)).get()
                    if os.path.isfile(s_path):
                        unique_count = calculate_fluidity(output_fname, q_path, s_path, total_count)
                        print '\t'.join([str(x) for x in [st_ix, comp_st_ix, total_count, unique_count]])
                        counts_df[st_ix][comp_st_ix] = (unique_count, total_count)
                    else:
                        strains_df = strains_df[strains_df.index!=comp_st_ix]
                        conn = engine.connect()
                        with conn.begin() as trans:
                            try:
                                conn.execute('DELETE FROM strains_plasmids WHERE strain_id=(?)', (comp_st_ix))
                                trans.commit()
                                print 'no FFN was found for %s (id %s), row deleted from strains_plasmids table' % (comp_strain_row['strain_name_uid'], comp_st_ix) 
                            except:
                                trans.rollback()
                                #print traceback.format_exc()
                        conn.close()
        else:
            strains_df = strains_df[strains_df.index!=st_ix]
            conn = engine.connect()
            with conn.begin() as trans:
                try:
                    conn.execute('DELETE FROM strains_plasmids WHERE strain_id=(?)', (st_ix))
                    trans.commit()
                    print 'no FFN was found for %s (id %s), row deleted from strains_plasmids table' % (strain_row['strain_name_uid'], st_ix) 
                except:
                    trans.rollback()
                    #print traceback.format_exc()
            conn.close()
    
    fluidity_df = pandas.DataFrame(index=strains_df.index.values, columns=strains_df.index.values)
    for ix in counts_df.index.values:
        for comp_ix in counts_df.index.values:
            if ix!=comp_ix:
                print tuple(counts_df[ix][comp_ix]), tuple(counts_df[comp_ix][ix])
                u1 = float(tuple(counts_df[ix][comp_ix])[0])
                u2 = float(tuple(counts_df[comp_ix][ix])[0])
                t1 = float(tuple(counts_df[ix][comp_ix])[1])
                t2 = float(tuple(counts_df[comp_ix][ix])[1])
                fluidity_df[ix][comp_ix] = (u1+u2)/(t1+t2)
    
    fluidity_df.to_sql('fluidity_scores', engine, if_exists='replace', index_label='strain_index')
    
    print 'done, %d strains were compared' % fluidity_df.index.size  
           
    return 

def plot_fluidity_vals(org_name, axes, n, ix):
    db_fname = os.path.join(OUTPUTDIR, org_name, 'output_datatables_%s.db' % org_name)
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    
    fluidity_df = pandas.read_sql_table('fluidity_scores', engine, index_col='strain_index', coerce_float=True)
    axes[ix%n, ix%2].set_title(org_name)
    fluidity_df.astype(float).plot(ax=axes[ix%n, ix%2], kind='box', rot=45)
    print fluidity_df.astype(float).max().max()
    axes[ix%n, ix%2].set_ylim([0,fluidity_df.astype(float).max().max()*1.05])
    return

def group_by_fluidity(org_name):
    '''for a given threshold, append strain to existing group only if its fluidity score with all existing strains in that group is <= threshold'''
    db_fname = os.path.join(OUTPUTDIR, org_name, 'output_datatables_%s.db' % org_name)
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    
    fluidity_df = pandas.read_sql_table('fluidity_scores', engine, index_col='strain_index').astype(float)
    print '%d genome files' % len(fluidity_df.index)
    if query_org not in CLONAL_ORGS: 
        thresholds = [0.05, 0.03, 0.01]
        unifiable_df = pandas.DataFrame(index=fluidity_df.index.values, columns=thresholds)
        
        for threshold in thresholds:
            print threshold
            count_unique = fluidity_df.index.size
            for strain_ix in fluidity_df.index.values:
                unifiable_df.loc[strain_ix,threshold] = []
                for comp_ix in fluidity_df.index.values:
                    if strain_ix!=comp_ix:
                        if fluidity_df.loc[strain_ix,comp_ix] < threshold:
                            unifiable_df.loc[strain_ix,threshold].append(comp_ix)
            
            for strain in unifiable_df.index.values:
                for comp_strain in unifiable_df.loc[strain,threshold]:
                    if strain not in unifiable_df.loc[comp_strain,threshold]:
                        unifiable_df.loc[strain,threshold].remove(comp_strain)
                for prev_strain in unifiable_df.index.values[:fluidity_df.index.get_loc(strain)]:
                    if prev_strain in unifiable_df.loc[strain, threshold]:
                        count_unique-=1
                        break
                unifiable_df.loc[strain,threshold] = ';'.join(unifiable_df.loc[strain,threshold])
            
            print '%d unique genomes for threshold %.2f' % (count_unique, threshold)
        #print unifiable_df
        unifiable_df.to_sql('unifiable_strains', engine, if_exists='replace', index_label='strain_id')
        print 'added unifiable_strains table to %s db' % org_name
    else:
        meta = sqlalchemy.MetaData()
        if 'unifiable_strains' in meta.tables.keys():
            unifiable_sql = sqlalchemy.Table('unifiable_strains', meta, autoload=True, autoload_with=engine)
            unifiable_sql.drop(engine, checkfirst=True)
            print 'removed unifiable_strains table from %s db' % org_name
    return

if __name__ == "__main__":
    org_dirs = {}
    for org_name in ORG_NAMES:
        org_dirs.update({org_name: []})
    '''
    matplotlib.style.use('ggplot')
    n = len(ORG_NAMES)-len(CLONAL_ORGS)
    fig, axes = plt.subplots(n//2, 2)
    fig.suptitle('Fluidity Thresholds per Strain (Non-Clonal Species Only)', fontsize=22)
    plt.subplots_adjust(hspace=0.8)
    '''      
    org_folders = sorted([name for name in os.listdir(GBKDIR) if os.path.isdir(os.path.join(GBKDIR, name))])
    i=0
    for query_org in org_dirs.keys():
        print query_org
        '''group_names = query_org.split('-')
        for group_name in group_names:
            for org_folder in sorted(org_folders):
                if re.search(group_name, org_folder):
                    org_dirs[query_org].append(org_folder)
        populate_strains_table(query_org, org_dirs[query_org])'''
        #compare_strains(query_org)
        group_by_fluidity(query_org)
        #if query_org not in CLONAL_ORGS:
            #i+=1
            #plot_fluidity_vals(query_org, axes, n//2, i)
    
    #plt.show()
        
