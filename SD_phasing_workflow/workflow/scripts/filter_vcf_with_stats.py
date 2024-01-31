# Input files are phased sample vcfs from phase_vcfs.smk
# This program do the following:
# count num of unphased sample sites in leaf vcf and pollen vcf
# count num of identically phased leaf vs. pollen sites
# count num of non-identically phased leaf vs. pollen sites
# count num of non-identically genotyped leaf vs. pollen sites
# filter out all sites that are not identically genotyped
# filter out all sites that are non Hardy-Weinberg between sample vs. parents
# calc avg. site depth, 95th, 5th percentiles
# report stats (.log file)
# rewrite subsetted vcf 
import sys
import numpy as np

l_vcf = sys.argv[1]
p_vcf = sys.argv[2]
l_sample_name = sys.argv[3]
p_sample_name = sys.argv[4]
final_lvcf = sys.argv[5]
final_pvcf = sys.argv[6]

# extract relevant stuff from vcf
def read_vcf(vcf, sample_name, filter):
    sample_dict = {}
    inconsistent = 0
    unphased = 0
    impossible_gts = {('00', '00', '01'):0, 
                        ('00', '00', '11'):0, 
                        ('00', '01', '11'):0,
                        ('00', '11', '00'):0,
                        ('00', '11', '11'):0,
                        ('01', '00', '11'):0,
                        ('01', '11', '00'):0,
                        ('11', '00', '00'):0,
                        ('11', '00', '11'):0,
                        ('11', '01', '00'):0,
                        ('11', '11', '00'):0,
                        ('11', '11', '01'):0}
    with open(vcf, 'r') as f:
        for line in f:
            if line.startswith("scaffold_") or line.startswith("Ahal_AUBY1_scaffold"):
                split = line.strip().split()
                #print(split)
                chrom = split[0].split("_")[-1]
                pos = split[1]
                sample_info = split[9].split(':')
                #print(sample_info)
                p01_info = split[10].split(':')
                p02_info = split[11].split(':')
                sample_gt = sample_info[0]
                # child genotype could not be determined? Weird.
                if sample_gt != "./.":
                    p01_gt = p01_info[0]
                    p02_gt = p02_info[0]
                    sample_dp = int(sample_info[2])
                    # if true, then we now filter out all inconsistent genotypes
                    if filter == True:
                        gt_set = (p01_gt[0]+p01_gt[2], p02_gt[0]+p02_gt[2], sample_gt[0]+sample_gt[2])
                        if "/" in sample_gt:
                            unphased += 1
                        elif gt_set in impossible_gts:
                            print(gt_set)
                            inconsistent += 1
                        else:
                            sample_dict[chrom + ':' + pos] = [sample_gt, sample_dp]

    sys.stderr.write(sample_name + "\t" + "number of unphased sites:" + "\t" + str(unphased) + "\n")
    sys.stderr.write(sample_name + "\t" + "number of consistent to parents sites:" + "\t" + str(len(sample_dict)) + "\n")
    sys.stderr.write(sample_name + "\t" + "number of inconsistent to parents sites:" + "\t" + str(inconsistent) + "\n")
    return sample_dict

# 
def filter_sites(l, p):

    identical_l = {}
    identical_p = {}
    diff_phase = 0
    diff_gt = 0
    unique_l = 0
    unique_p = 0
    l_depths = []
    p_depths = []

    for site in l:
        l_depths.append(l[site][1])
        if site in p:
            l_gt = l[site][0]
            p_gt = p[site][0]
            # genotype identical and phased identically
            if l_gt == p_gt:
                identical_l[site] = l[site]
                identical_p[site] = p[site]
            # genotype identical but opposite phasing
            elif {l_gt[0], l_gt[2]} == {p_gt[0], p_gt[2]}:
                diff_phase += 1
            # different alleles at this site           
            else:
                diff_gt += 1
        else:
            unique_l += 1

    for site in p:
        p_depths.append(p[site][1])

        if site not in l:
            unique_p += 1
    
    sys.stderr.write("number of identically phased sites:" + "\t" + str(len(identical_l)) + "\n")
    sys.stderr.write("number of non-identically phased sites:" + "\t" +  str(diff_phase) + "\n")
    sys.stderr.write("number of non-identically genotyped sites:" + "\t" + str(diff_gt) + "\n")

    return identical_l, identical_p, l_depths, p_depths

def percentile_filter(l, p, l_depths, p_depths):

    l_avg = np.average(l_depths)
    p_avg = np.average(p_depths)
    l_5 = np.percentile(l_depths, 5)
    l_95 = np.percentile(l_depths, 95)
    p_5 = np.percentile(p_depths, 5)
    p_95 = np.percentile(p_depths, 95)
    final_sites = {}

    # get rid of sites below 5th percentile or above 95th percentile
    for site in list(l):
        if l[site][1] > l_95 or l[site][1] < l_5:
            l.pop(site)
    for site in list(p):
        if p[site][1] > p_95 or p[site][1] < p_5:
            p.pop(site)
    # grab final 
    for site in l:
        if site in p:
            final_sites[site] = 0
        
    sys.stderr.write("leaf avg depth:" + "\t" + str(l_avg) + "\n")
    sys.stderr.write("leaf 5th and 95th percentiles:" + "\t" + str(l_5) + "\t" + str(l_95) + "\n")
    sys.stderr.write("pollen avg depth:" + "\t" + str(p_avg) + "\n")
    sys.stderr.write("pollen 5th and 95th percentiles:" + "\t" + str(p_5) + "\t" + str(p_95) + "\n")
    sys.stderr.write("final number of indentical sites post-depth filtering:" + "\t" + str(len(final_sites)) + "\n")
    return final_sites

# write new vcf
def make_final_vcf(final_sites, vcf, final_vcf):
    with open(vcf, 'r') as f, open(final_vcf, 'w') as o:
        for line in f:
            if line.startswith("scaffold_") or line.startswith("Ahal_AUBY1_scaffold"):
                split = line.strip().split()
                maybe_key = split[0].split("_")[-1] + ":" + split[1]
            
                if maybe_key in final_sites:
                    o.write(line)



sys.stderr.write("phasing stats for cross " + "\t" + l_sample_name[0:3] + "\n")
l_sites = read_vcf(l_vcf, l_sample_name, True)
p_sites = read_vcf(p_vcf, p_sample_name, True)

id_l_sites, id_p_sites, l_depths, p_depths = filter_sites(l_sites, p_sites)
final_sites = percentile_filter(id_l_sites, id_p_sites, l_depths, p_depths)

make_final_vcf(final_sites, l_vcf, final_lvcf)
make_final_vcf(final_sites, p_vcf, final_pvcf)

