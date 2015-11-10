import sqlite3, os, re, numpy, pandas, sqlalchemy
from Bio import SeqIO

OUTPUTDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/output'
FFNDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/all_ffn'
GBKDIR = '/home/michalsa/Documents/research/codon_bias/ncbi_ftp_genomes/all_gbk'

ORG_NAMES = ['Bacillus_cereus', 'Campylobacter_jejuni', 'Clostridium_botulinum', 'Corynebacterium_diphtheriae',
             'Corynebacterium_pseudotuberculosis', 'Escherichia_coli-Shigella', 'Francisella_tularensis', 
             'Haemophilus_influenzae', 'Listeria_monocytogenes', 'Mycobacterium_tuberculosis', 'Neisseria_meningitidis',
             'Pseudomonas_aeruginosa', 'Streptococcus_pneumoniae', 'Yersinia_pestis']

THRESHOLDS = [0.1, 0.08, 0.05, 0.03, 0.02]
IDENTITY = 90

VIR_KEYWORDS = ['viral', 'virus', 'phage', 'capsid', 'tail']
UNKNOWN_KEYWORDS = ['hypothetical protein', 'uncharacterized protein']

def db_init(filename):
    con = sqlite3.connect(filename)
    with con:
        cur = con.cursor()
        
        con.execute("PRAGMA foreign_keys=ON")
        
        cur.execute("CREATE TABLE IF NOT EXISTS strains_plasmids(strain_id TEXT PRIMARY KEY, strain_name_uid TEXT, plasmid_ids TEXT, plasmid_names TEXT) WITHOUT ROWID;")
        
        cur.execute("CREATE TABLE IF NOT EXISTS paralogs(pangene_id TEXT PRIMARY KEY NOT NULL, paralog_ids TEXT) WITHOUT ROWID;")
        
        threshold_strs = []
        for threshold in THRESHOLDS:
            threshold_str = "%sEXPneg%d" % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
            threshold_strs.append("fluidity_%s TEXT" % threshold_str)
            cur.execute("CREATE TABLE IF NOT EXISTS pangenome_%s(pangene_id TEXT PRIMARY KEY, ortholog_ids TEXT) WITHOUT ROWID;" % threshold_str)
            cur.execute("CREATE TABLE IF NOT EXISTS pangenome_%s_stats(pangene_id TEXT PRIMARY KEY NOT NULL REFERENCES pangenome_%s(pangene_id), n_orthologs INT, n_genomes INT, \
            paralog_flag INT, viral_flag INT, enc REAL, enc_std REAL, encp REAL, encp_std REAL, fop REAL, fop_std REAL) WITHOUT ROWID;" 
            % (threshold_str, threshold_str))
        
        cur.execute("CREATE TABLE IF NOT EXISTS unifiable_strains(strain_id TEXT PRIMARY KEY NOT NULL REFERENCES strains_plasmids(strain_id), %s) WITHOUT ROWID;" 
                    % ", ".join(threshold_strs))
    return

def db_show(filename):
    con = sqlite3.connect(filename)
    with con:
        cur = con.cursor()
        '''show existing data'''
        cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
        #print [d[0] for d in cur.description]
        for row in cur.fetchall():
            print row
            cur.execute("PRAGMA table_info(%s)" % row)
            #cur.execute("SELECT * FROM %s" % row)
            '''print column info'''
            for r in cur.fetchall():
                print r
    return

def insert_data_from_file(db_filename, table_name, data_filename):
    con = sqlite3.connect(db_filename)
    with con:
        cur = con.cursor()
        input_fh = open(data_filename, 'r')
        count = 0
        with input_fh:
            for line in input_fh:
                count+=1
                row = line.strip().split('\t')
                insert_str = '(%s)' % ','.join([el for el in row])
                cur.execute('INSERT INTO ? VALUES ?', table_name, insert_str)
        
        print 'table %s successfully updated, %d values inserted' % (table_name, count)
    return

def insert_strain_info(db_filename, data_filename):
    con = sqlite3.connect(db_filename)
    with con:
        cur = con.cursor()
        input_fh = open(data_filename, 'r')
        count = 0
        with input_fh:
            for line in input_fh:
                row = line.strip().split(',')
                if row[0]!='strain':
                    count+=1
                    st_key = row[2]
                    st_name_uid = '%s_uid%s' % (row[0], row[1])
                    pl_key = None
                    pl_name = None
                    if len(row)>3:
                        pl_key = row[4]
                        pl_name = row[3]
                    insert_lst = [st_key, st_name_uid, pl_key, pl_name]
                    cur.execute('INSERT INTO strains_plasmids VALUES (?, ?, ?, ?)', insert_lst)
                    con.commit()
        
        print 'table strains_plasmids successfully updated, %d entries inserted' % count
    return

def insert_unifiables_info(db_filename, data_filename):
    con = sqlite3.connect(db_filename)
    with con:
        cur = con.cursor()
        cur.execute("SELECT * FROM unifiable_strains;")
        cur.execute("DELETE FROM unifiable_strains;")
        input_fh = open(data_filename, 'r')
        count = 0
        with input_fh:
            for line in input_fh:
                row = line.strip().split(',')
                if row[0]:
                    st_name_uid = row[0]
                    #print unicode(st_name_uid, "utf-8")
                    cur.execute("SELECT strain_id FROM strains_plasmids WHERE strain_name_uid=(?)", [unicode(st_name_uid, "utf-8")])
                    #print cur.fetchone()[0]
                    st_id = cur.fetchone()[0]
                    if st_id==None:
                        print 'entry unsuccessful for strain: %s - id does not exist in strains_plasmids table' % st_name_uid
                    else:
                        uni_st_ids = [None]*len(row[1:])
                        for ix, el in enumerate(row[1:]):
                            if el!='':
                                uni_st_ids[ix] = []
                                strains_fluidity = el.split(';')
                                for strain in strains_fluidity:
                                    print strain
                                    cur.execute("SELECT strain_id FROM strains_plasmids WHERE strain_name_uid=(?)", [unicode(strain, "utf-8")])
                                    uni_st_id = cur.fetchone()[0]
                                    if uni_st_id==None:
                                        print 'entry unsuccessful for strain: %s - id does not exist in strains_plasmids table' % uni_st_id
                                    else:
                                        uni_st_ids[ix].append(uni_st_id)
                        insert_lst = [st_id]+[el if el!=[] else None for el in uni_st_ids]
                        insert_lst = [';'.join(el) if type(el)==list else el for el in insert_lst]
                        cur.execute('INSERT INTO unifiable_strains VALUES (?, ?, ?, ?, ?, ?)', insert_lst)
                        count+=1
                        
        print 'table unifiable_strains successfully updated, %d entries inserted' % count
    return

def insert_paralogs_info(db_filename, data_filename):
    con = sqlite3.connect(db_filename)
    with con:
        cur = con.cursor()
        input_fh = open(data_filename, 'r')
        count = 0
        with input_fh:
            for line in input_fh:
                row = line.strip().split('\t')
                count+=1
                cur.execute('INSERT INTO paralogs VALUES (?, ?)', row)
                con.commit()
        
        print 'table paralogs successfully updated, %d entries inserted' % count
    return

def insert_pangenome_info(db_filename):
    con = sqlite3.connect(db_filename)
    with con:
        for threshold in THRESHOLDS:
            print threshold
            threshold_str_sql = "%sEXPneg%d" % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
            threshold_str_file = "%se-%d" % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
            
            data_filename = os.path.join(OUTPUTDIR, ORGNAME, 'pangenome_summary_identity%d_fluidity%s_corrected_afterselfcmp.tsv' % (IDENTITY, threshold_str_file))
            
            cur = con.cursor()
            input_fh = open(data_filename, 'r')
            count = 0
            with input_fh:
                for line in input_fh:
                    row = line.strip().split('\t')
                    if len(row)>1:
                        orthologs = row[1].split(';')
                        orthologs_unique = []
                        for ortholog in orthologs:
                            if ortholog!=row[0] and ortholog not in orthologs_unique:
                                orthologs_unique.append(ortholog)
                    count+=1
                    insert_lst = [row[0], ';'.join(orthologs_unique)]
                    
                    cur.execute('INSERT INTO pangenome_{} VALUES (?, ?)'.format(threshold_str_sql), insert_lst)
    
            print 'table pangenome_{thr} successfully updated, {cnt} entries inserted'.format(thr=threshold_str_sql, 
                                                                                              cnt=count)
    return

def populate_stats_tables(db_filename):
    con = sqlite3.connect(db_filename)
    with con:
        for threshold in THRESHOLDS:
            print threshold
            threshold_str = "%sEXPneg%d" % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
            
            cur = con.cursor()
            pangenes = cur.execute("SELECT * FROM pangenome_{}".format(threshold_str))
            count = 0
            for row in pangenes:
                pangene_id = row[0]
                orthologs = row[1].split(';')
                
                genes = [pangene_id] + orthologs
                
                n_orthologs = len(genes)+1
                
                genes_lengths = []
                genomes_unique = []
                enc_vals = []
                encp_vals = []
                for gene in genes:
                    genome_id = re.search('\|(NC_\d+)\.\d\|', gene).group(1)
                    
                    '''get unifiable strains for gene using its genome_id'''
                    q_uni = cur.execute("SELECT fluidity_{} FROM unifiable_strains WHERE strain_id=(?)".format(threshold_str), [genome_id]).fetchall()
                    unifiables = []
                    for unifiable in q_uni:
                        for u in unifiable:
                            unifiables.append(u)
                    
                    '''get strain name & uid for gene using its genome_id'''
                    q_strain = cur.execute("SELECT strain_name_uid FROM strains_plasmids WHERE strain_id=(?)".format(threshold_str), [genome_id]).fetchall()
                    st_name_uid = ''
                    if len(q_strain)==1: st_name_uid = q_strain[0][0]
                    else:
                        raise ValueError('error: more than one strain found with id %s' % genome_id)
                        return
                    
                    '''get gene length'''
                    if genome_id not in genomes_unique and [unifiable not in genomes_unique for unifiable in unifiables]==[True]*len(unifiables):
                        genomes_unique.append(genome_id)
                    else:
                        print gene
                    pos_str = re.search('\|:c?(\d+-\d+)', gene).group(1)
                    pos_arr = [int(p) for p in pos_str.split('-')]
                    genes_lengths.append(abs(pos_arr[0]-pos_arr[1]))
                    
                    '''get enc, encprime values for gene'''
                    encprime_output_fname = os.path.join(OUTPUTDIR, ORGNAME, 'encprime', '_'.join([ORGNAME, st_name_uid]), '%s.ffn.results' % genome_id)
                    encprime_output_fh = open(encprime_output_fname, 'r')
                    for encprime_line in encprime_output_fh:
                        if re.search('^(gi[^\s]+) ', encprime_line):
                            encprime_gene_name = re.search('^(gi[^\s]+) ', encprime_line).group(0)
                            if re.search(gene.split('|')[-1], encprime_gene_name):
                                encprime_trimmed_line = re.search(': (.*)$', encprime_line).group(1)
                                encprime_nums = [float(n) if re.search('\.', n) else int(n) for n in encprime_trimmed_line.split(' ')]
                                enc_vals.append(encprime_nums[0])
                                encp_vals.append(encprime_nums[1])
                
                n_genomes = len(genomes_unique)
                
                l_gene = numpy.mean(genes_lengths)
                l_gene_std = numpy.std(genes_lengths)
                
                enc = numpy.mean(enc_vals)
                enc_std = numpy.std(enc_vals)
                
                encp = numpy.mean(encp_vals)
                encp_std = numpy.std(encp_vals)
                
                print pangene_id, l_gene, l_gene_std, n_orthologs, n_genomes, enc, enc_std, encp, encp_std
                return
                
                count+=1
                insert_lst = [row[0], ';'.join(orthologs_unique)]
                
                #cur.execute('INSERT INTO pangenome_stats_{} VALUES (?, ?)'.format(threshold_str), insert_lst)
    
            print 'table pangenome_stats_{thr} successfully updated, {cnt} entries inserted'.format(thr=threshold_str, 
                                                                                              cnt=count)
    return

def populate_annotations_table(db_filename, org_name):
    engine = sqlalchemy.create_engine('sqlite:///{}'.format(db_fname))
    strains_df = pandas.read_sql_table('strains_plasmids', engine, index_col='strain_id')
    
    engine.execute('CREATE TABLE IF NOT EXISTS gene_annotations(gene_id TEXT PRIMARY KEY NOT NULL,\
                                                                annotation TEXT NOT NULL,\
                                                                annotation_full TEXT NOT NULL) WITHOUT ROWID;')
    
    col_data = engine.execute('PRAGMA table_info(gene_annotations)').fetchall()
    
    if True not in [col[1]=='annotation_full' for col in col_data]:
        engine.execute('ALTER TABLE gene_annotations ADD COLUMN annotation_full TEXT')
    
    counter = 0
        
    for strain_row in strains_df.iterrows():
        #print strain_row
        gb_filename = os.path.join(GBKDIR, '%s_%s' % (org_name, strain_row[1]['strain_name_uid']), '%s.gbk' % strain_row[0])
        with open(gb_filename, 'r') as gb_fh:
            for record in SeqIO.parse(gb_fh, 'genbank'):
                for feature in record.features:
                    if feature.type=='CDS':
                        loc_components = [':', '', '%s-%s' % (str(feature.location.start+1), str(feature.location.end))]
                        if feature.location.strand==-1:
                                loc_components[1]='c'
                                loc_components[2] = '%s-%s' % (str(feature.location.end), str(feature.location.start+1))
                        loc_str = ''.join(loc_components)
                        gene_id = '|'.join(['gi', record.annotations['gi'], 'ref', record.id, loc_str])
                        
                        if engine.execute('SELECT annotation_full FROM gene_annotations WHERE gene_id=(?)', (gene_id,)).fetchone() == (None,):
                            print gene_id
                            
                            ann_str = ''
                            
                            if 'product' in feature.qualifiers.keys(): 
                                ann_str = ';'.join(feature.qualifiers['product'])
                            
                            elif 'pseudo' in feature.qualifiers.keys():
                                ann_str = ';'.join(feature.qualifiers['note'])
                            
                            origin = 'bacterial'
                            for vir_keyword in VIR_KEYWORDS:
                                if re.search(vir_keyword, ann_str, re.IGNORECASE):
                                    origin = 'viral'
                                    break
                            for unknown_keyword in UNKNOWN_KEYWORDS:
                                if re.search(unknown_keyword, ann_str, re.IGNORECASE):
                                    origin = 'unknown'
                                    break
                            
                            print origin, '|', ann_str
                                
                            #gene absent from annotations table
                            if len(engine.execute('SELECT * FROM gene_annotations WHERE gene_id=(?)', (gene_id,)).fetchall())==0:
                                conn = engine.connect()
                                with conn.begin() as trans: 
                                    try: 
                                        conn.execute('INSERT INTO gene_annotations VALUES (?,?,?)', (gene_id, origin, ann_str))
                                        trans.commit()
                                    except:
                                        trans.rollback()
                                        raise
                                conn.close()
                            
                            #gene exists in annotations table, add full annotation text and update origin
                            else: 
                                conn = engine.connect()
                                with conn.begin() as trans: 
                                    try: 
                                        conn.execute('UPDATE gene_annotations SET annotation=(?), annotation_full=(?) WHERE gene_id=(?)', (origin,ann_str,gene_id))
                                        trans.commit()
                                        counter+=1
                                    except:
                                        trans.rollback()
                                        raise
                                conn.close()
                                
    #print '%d gene entries updated' % engine.execute('SELECT Count(*) FROM gene_annotations').fetchone()[0]
    print '%d gene entries updated' % counter
    return

def empty_table(db_filename, table_name):
    con = sqlite3.connect(db_filename)
    with con:
        cur = con.cursor()
        cur.execute("SELECT * FROM {};".format(table_name))
        cur.execute("DELETE FROM {};".format(table_name))
        print 'table %s successfully emptied' % table_name

if __name__ == '__main__':
    for name in ORG_NAMES:
        db_fname = os.path.join(OUTPUTDIR, name, 'output_datatables_%s.db' % name)
        #db_init(db_fname)
        #db_show(db_fname)
        #insert_strain_info(db_fname, os.path.join(OUTPUTDIR, name, 'strains_plasmids_ncbiFTP.csv'))
        #insert_unifiables_info(db_fname, os.path.join(OUTPUTDIR, name, 'unifiable_strains_revised.csv'))
        #insert_paralogs_info(db_fname, os.path.join(OUTPUTDIR, name, 'paralogs.tsv'))
        #insert_pangenome_info(db_fname)
        #populate_stats_tables(db_fname)
        #populate_annotations_table(db_fname, name)
        '''for threshold in THRESHOLDS:
            t_name = 'pangenome_%sEXPneg%d' % (re.search('([1-9])', str(threshold)).group(1), str(threshold).count('0'))
            empty_table(db_fname, t_name)'''
    
    