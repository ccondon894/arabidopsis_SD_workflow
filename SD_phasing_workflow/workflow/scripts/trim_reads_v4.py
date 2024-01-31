import sys
import gzip
from xopen import xopen
from collections import defaultdict
# My own version of Russ' trim_reads script. It seems to work more correctly?
# Updates:
# Trying xopen to see if it is faster than gzip.open
# Now prints retrimmed fastq and redirects to correct file using redirect option within print()
# Fixed the error that would skip reads of certain sizes  

pollen = sys.argv[1]
leaf = sys.argv[2]

pollen_sample_name = sys.argv[3]
leaf_sample_name = sys.argv[4]

path = sys.argv[5]

# get read length dist for smaller sample 
def collect_read_dist(sample, lengths):
    
    with xopen(sample) as f:
        while True:
            header = f.readline().strip()
            sequence = f.readline().strip()
            trash = f.readline().strip()
            qual = f.readline().strip()
            
            read_len = len(sequence)
            lengths[read_len] += 1

            if not header: 
                break
    
    return lengths

# 
def trim_other(pollen_lengths, leaf, outfile):
    with xopen(leaf) as f, open(outfile, 'w') as o:
        while True:
            header = f.readline().strip()
            sequence = f.readline().strip()
            trash = f.readline().strip()
            qual = f.readline().strip()

            read_len = len(sequence)
            for i in reversed(range(20, read_len + 1)):
                if i in pollen_lengths and pollen_lengths[i] > 0:

                    print(header, file=o)
                    print(sequence[0:i], file=o)
                    print(trash, file=o)
                    print(qual[0:i], file=o)
                    pollen_lengths[i] -= 1
                    if pollen_lengths[i] == 0:
                        pollen_lengths.pop(i)
                    break
            
            if len(pollen_lengths.keys()) == 0:
                break
            if not header: 
                break

        o.close()

def run_it(leaf, pollen, lname, pname, path):
    
    l_outfile = path + lname + ".retrimmed.fastq" 
    p_outfile = path + pname + ".retrimmed.fastq"

    pollen_lengths = defaultdict(int)
    leaf_lengths = defaultdict(int)

    pollen_lengths = collect_read_dist(pollen, pollen_lengths)
    leaf_lengths = collect_read_dist(leaf, leaf_lengths)

    l_counts = 0
    p_counts = 0
    for i in reversed(range(130, 152)):
        l_counts += leaf_lengths[i]
        p_counts += pollen_lengths[i]

    if l_counts >= p_counts:
        trim_other(pollen_lengths, leaf, l_outfile)
    elif l_counts < p_counts:
        trim_other(leaf_lengths, pollen, p_outfile)

run_it(leaf, pollen, leaf_sample_name, pollen_sample_name, path)







