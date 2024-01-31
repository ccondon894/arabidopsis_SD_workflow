# This program do the following:
# Open sample vcf and collect everything in various lists into dict
# Open each chr mpileup vcf with retrimmed + remapped variant counts
# Update the positions in the dict to new counts
# make new vcf with new counts
import sys
import os
import gzip

vcf = sys.argv[1]
sample_name = sys.argv[2]
mpileup_path = sys.argv[3]

mpileup_files = os.listdir(mpileup_path)
mp_vcfs = [file for file in mpileup_files if (file.endswith("vcf.gz"))] # list of mpileup vcfs

def read_vcf(vcf): # get vcf lines, one key per line
    sample_dict = {}
    with open(vcf, 'r') as f:
        for line in f:
            split = line.strip().split()
            chrom = split[0].split("_")[-1] # this should work with all species (I hope)
            pos = split[1]
            info = split[9].split(':')
            # split[0:9] is everything up until the sample gentoype info
            # info is the list of genotype info split by :
            # split[10:] is everything else
            # last spot is boolean to keep track of whether position was updated or not
            sample_dict[chrom+":"+pos] = [split[0:9], info, split[10:], False]
            
    return sample_dict


def read_mpileup(chrom_vcf, sample_dict, sample): # read mpileup for each chr

    with gzip.open(chrom_vcf, 'rt') as f:
        count = 0
        for line in f:
            if line.startswith('#CHROM'):
                # get index for current sample
                format_line = line.split("\t")
                samples = [i for i in format_line[9: len(format_line)]]
                samples = [i.split("/")[1] for i in samples]
                samples = [i[0:5] for i in samples]
                sample_idx = samples.index(sample)

            elif line.startswith("scaffold_") or line.startswith("Ahal_AUBY1_scaffold"):
                split = line.strip().split()
                chrom = split[0].split("_")[-1]
                pos = split[1]
                key = chrom+":"+pos
                if key in sample_dict:
                    # update counts
                    info = split[9 + sample_idx].split(":") # type: ignore
                    # get mpileup line genotype
                    mp_gt = info[0]
                    # get mpileup line allele depths
                    mp_ad = info[5].split(',')
                    ad1 = mp_ad[0]
                    ad2 = mp_ad[1]
                    # store sample line genotype in easy to compare format
                    samp_gt = (sample_dict[key][1][0][0], sample_dict[key][1][0][2])
                    # make sure genotypes are consistent between files 
                    # if genotypes are not consistent, skip that site
                    if (mp_gt[0], mp_gt[2]) == samp_gt and (mp_gt[0], mp_gt[2]) == ('0','0'):
                        sample_dict[key][1][1] = ad1 + ',' + ad2
                        sample_dict[key][3] = True
                        count += 1
                    elif (mp_gt[0], mp_gt[2]) == samp_gt and (mp_gt[0], mp_gt[2]) == ('1','1'):
                        sample_dict[key][1][1] = ad1 + ',' + ad2
                        sample_dict[key][3] = True
                        count += 1
                    # also check phasing for het genotypes and flip allele depths if need be
                    elif (mp_gt[0], mp_gt[2]) == samp_gt and (mp_gt[0], mp_gt[2]) == ('0','1'):
                        sample_dict[key][1][1] = ad1 + ',' + ad2
                        sample_dict[key][3] = True
                        count += 1
                    elif (mp_gt[2], mp_gt[0]) == samp_gt and (mp_gt[0], mp_gt[2]) == ('0','1'):
                        sample_dict[key][1][1] = ad2 + ',' + ad1
                        sample_dict[key][3] = True
                        count += 1
    

def write_vcf(sample_dict): # write new vcf with updated allele depths

    for key in sample_dict:
        # if site was updated, convert lists back into string and write to file
        if sample_dict[key][3] == True:
            chunk_1 = "\t".join(sample_dict[key][0])
            chunk_2 = ":".join(sample_dict[key][1])
            chunk_3 = "\t".join(sample_dict[key][2])
            print(chunk_1 + "\t" + chunk_2 + "\t" + chunk_3)

sample = read_vcf(vcf)
for mp in mp_vcfs:
    read_mpileup(mpileup_path + mp, sample, sample_name)

write_vcf(sample)