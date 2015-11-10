import csv, os, re
from collections import defaultdict

import ggplot as gplt
import pandas as pd

import matplotlib.pyplot as mplt

ROOTDIR = '/home/michalsa/Documents/research/metagenomes/humanGut'

def phylo_stats(filesDict, outfile=os.path.join(ROOTDIR, 'phylogeneticDist.csv')):
    
    phyloDict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for key, val in filesDict.items():
        print key
        phyloDict[key] = defaultdict(lambda: defaultdict(int))
        with open(val, 'r') as input_fh:
            phyloReader = csv.reader(input_fh, delimiter='\t')
            for row in phyloReader:
                lineage = row[-1]
                linArr = lineage.split(';')
                #if float(row[3]) > 50:
                if linArr[0] == 'Viruses':
                    #print linArr[4]
                    phyloDict[key][linArr[0]][linArr[4]]+=1
                else:
                    phyloDict[key][linArr[0]][linArr[1]]+=1
                    
    with open(outfile, 'w') as output_fh:
        phyloWriter = csv.writer(output_fh)
        phyloWriter.writerow(['metagenome', 'domain', 'phylum', '#contigs', '%contigs'])
        for kPhylo in sorted(phyloDict.keys()):
            #phyloWriter.writerow([kPhylo])
            for kDom in sorted(phyloDict[kPhylo].keys()):
                #phyloWriter.writerow(['', kDom])
                totPhy = sum([vPhy for vPhy in phyloDict[kPhylo][kDom].values()])
                #print totPhy
                for kPhy, vPhy in sorted(phyloDict[kPhylo][kDom].items()):
                    fracPhy = float(vPhy)/totPhy
                    #print '%d/%d = %f' % (vPhy, totPhy, fracPhy)
                    phyloWriter.writerow([kPhylo, kDom, kPhy, vPhy, fracPhy])
       
    #phyloDF = pd.DataFrame.from_csv(outfile, index_col=None)
    #print phyloDF
        
def phylo_barcharts(keys, outfile=os.path.join(ROOTDIR, 'phylogeneticDist.csv')):
    statsDF = pd.DataFrame.from_csv(outfile, index_col=None)
    orgs = ['Bacteria', 'Archaea', 'Viruses']
    '''
    df = statsDF[statsDF['metagenome'].isin(keys) & statsDF['domain'].isin(orgs)]
    print plt.ggplot(plt.aes(x=u'phylum', y=u'#contigs'), data=df) + \
    plt.ggtitle(', '.join([key for key in keys])+'_Phylogenetic Distribution') + \
    plt.geom_bar(stat='bar')
    '''
    for org in orgs:
        df = pd.DataFrame.from_records(statsDF[(statsDF['metagenome'].isin(keys)) & (statsDF['domain'] == org) & (statsDF['%contigs'] > 0.01)])
        print df
        
        print gplt.ggplot(gplt.aes(x='phylum', y='%contigs', fill='metagenome'), data=df) + \
        gplt.ggtitle(org+'_'+';'.join([key for key in keys])+'_Phylogenetic Distribution') + \
        gplt.geom_bar(stat='bar', position='dodge')

def phylo_piecharts(keys=None, outfile=os.path.join(ROOTDIR, 'phylogeneticDist_corr.csv')):   
    statsDF = pd.DataFrame.from_csv(outfile, index_col=None) 
    
    if keys==None:
        keys = list(set(statsDF['metagenome']))
        
    doms = ['Bacteria', 'Archaea', 'Viruses']  
        
    for dom in doms:
        phyDict = defaultdict(lambda: defaultdict(int)) 
        for key in keys:
            phys = list(set(statsDF[(statsDF['metagenome'] == key) & (statsDF['domain'] == dom)]['phylum']))
            for phy in phys:
                df = statsDF[(statsDF['metagenome'] == key) & (statsDF['domain'] == dom) & (statsDF['phylum'] == phy)]      
                #print df          
                phyDict[key][phy] = sum(df['#contigs'])
            
        phyDF = pd.DataFrame.from_dict(phyDict)
        #print phyDF
                
        df_pie = phyDF.sum(1)
        other_sum = 0
        #print [df_pie[key] for key in df_pie.keys()]
        b_arr = [df_pie[key]>0.01*df_pie.sum() for key in df_pie.keys()]
        #print b_arr
        for ix, b in enumerate(b_arr):
            if not b:
                #print ix
                #print df_pie[ix]
                other_sum+=df_pie[ix]
                
        df_pie_new = {}
        for key in df_pie.keys():
            if df_pie[key]>0.01*df_pie.sum():
                df_pie_new[key] = df_pie[key]
        if other_sum > 0:
            df_pie_new['other'] = other_sum
        
        print df_pie_new
        
        fig = mplt.figure()
        fig.suptitle('%s_%s_phyla\n(total: %d)' % (ROOTDIR.split('/')[-1], dom, df_pie.sum()), fontsize=14, fontweight='bold')
        pieplt = fig.add_subplot(111) 
        pieplt.pie(df_pie_new.values(), labels=df_pie_new.keys(), autopct='%1.2f%%', pctdistance=1.0)
        mplt.legend(loc=2)
        # Set aspect ratio to be equal so that pie is drawn as a circle.
        pieplt.axis('equal')
        mplt.show()
        #fig.savefig(os.path.join(ROOTDIR, 'lineP_%s_phylums.png' % dom), dpi=fig.dpi)  
             
def dom_piecharts(keys=None, outfile=os.path.join(ROOTDIR, 'phylogeneticDist_corr.csv')):
    statsDF = pd.DataFrame.from_csv(outfile, index_col=None) 
    
    if keys==None:
        keys = list(set(statsDF['metagenome']))
        #print list(set(keys))
    
    #orgs = ['Bacteria', 'Archaea', 'Viruses']
    
    domDict = defaultdict(lambda: defaultdict(int))  
    for key in keys:
        doms = list(set(statsDF[statsDF['metagenome'] == key]['domain']))
        #print orgs
        #tot = sum(statsDF[statsDF['metagenome'] == key]['#contigs'])
        #print tot
        for dom in doms:
            df = statsDF[(statsDF['metagenome'] == key) & (statsDF['domain'] == dom)]
            dom_tot = sum(df['#contigs'])
            #if dom_tot > 0.01*tot:
            domDict[key][dom] = dom_tot
            '''
            else:
                domDict[key]['other'] += dom_tot
            '''
        #print df['phylum'], df['%contigs'] 
    
    domDF = pd.DataFrame.from_dict(domDict)
    
    df_pie = domDF.sum(1)
    other_sum = 0
    #print [df_pie[key] for key in df_pie.keys()]
    b_arr = [df_pie[key]>0.01*df_pie.sum() for key in df_pie.keys()]
    #print b_arr
    for ix, b in enumerate(b_arr):
        if not b:
            #print ix
            #print df_pie[ix]
            other_sum+=df_pie[ix]
    
    df_pie_new = {}
    for key in df_pie.keys():
        if df_pie[key]>0.01*df_pie.sum():
            df_pie_new[key] = df_pie[key]
    if other_sum > 0:
        df_pie_new['other'] = other_sum
    print df_pie_new
    
    fig = mplt.figure()
    fig.suptitle('%s_domains\n(total: %d)' % (ROOTDIR.split('/')[-1], df_pie.sum()), fontsize=14, fontweight='bold')
    pieplt = fig.add_subplot(111) 
    pieplt.pie(df_pie_new.values(), labels=df_pie_new.keys(), autopct='%1.2f%%')
    mplt.legend(loc=2)
    # Set aspect ratio to be equal so that pie is drawn as a circle.
    pieplt.axis('equal')
    mplt.show()
    #fig.savefig(os.path.join(ROOTDIR, 'lineP_domains.png'), dpi=fig.dpi)     
    
    
if __name__ == '__main__':
    '''
    phyloFiles = {}
    for fname in os.listdir(ROOTDIR):
        if os.path.isdir(os.path.join(ROOTDIR, fname)):
            MGDIR = fname
            patt = re.compile('(%s\.[au])\.phylodist\.txt' % MGDIR)
            for mgfname in os.listdir(os.path.join(ROOTDIR, MGDIR)):
                phylofname = patt.match(mgfname)
                if phylofname:
                    phyloFiles['%s' % phylofname.group(1)] = os.path.join(ROOTDIR, MGDIR, mgfname)
    
    phylo_stats(phyloFiles)
    '''
    #keysList = ['3300000149.a', '3300001472.a']
    #phylo_barcharts(keysList)
    phylo_piecharts()
    #dom_piecharts()
    
    