from Bio import SeqIO
import os, re
def convert_divide():
    filename = '/home/michalsa/Downloads/ncbi-blast-2.2.30+/db/2219.1.1822.CCCATG.fastq'
    handle = open(filename, 'r')
    count = 1
    
    for record in SeqIO.parse(handle, 'fastq'):
        out_filename = '/home/michalsa/Downloads/ncbi-blast-2.2.30+/db/2219.1.1822.CCCATG_%d.fasta' % count
        if not os.path.isfile(out_filename) or os.stat(out_filename).st_size < 1000000001:
            out_handle = open(out_filename, 'a')
            SeqIO.write(record, out_handle, 'fasta')
        else:
            count+=1

def locate():    
    search_terms = ['2219:1:2204:20169:54082/2', '2219:1:2204:20169:54082/1', '2219:1:2202:16762:51632/1', '2219:1:1203:4975:91982/2', '2219:1:1203:4975:91982/1', 
                    '2219:1:2204:20169:54082/2', '2219:1:2204:20169:54082/1', '2219:1:2202:16762:51632/1', '2219:1:1203:4975:91982/2', '2219:1:1203:4975:91982/1']
    
    for term in search_terms:
        for num in range(1, 11):
            filename = '/home/michalsa/Downloads/ncbi-blast-2.2.30+/db/2219.1.1822.CCCATG_%d.fasta' % num
            handle = open(filename, 'r')
            for record in SeqIO.parse(handle, 'fasta'):
                if re.search(term, record.id):
                    print '%s: %s' % (term, filename)
                    break
                
if __name__ == "__main__":
    locate()