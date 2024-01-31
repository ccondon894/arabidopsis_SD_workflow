import sys
# Simple program grabs phasing info from Whatshap phased vcf files
# Output is a .haplotypes file in the following format:
# (chr, pos, ref/alt, ref/alt, ref_depth/alt_depth, ref_depth/alt_depth)


l_file = sys.argv[1]
p_file = sys.argv[2]

def extract_vcf_info(file):
    prev_pos = 0
    sample_dict = {}
    with open(file, 'r') as f:
        for line in f:
            if line.startswith("scaffold_") or line.startswith("Ahal_AUBY1_scaffold"):
                line_info = line.strip().split("\t")
                chr = line_info[0].split("_")[-1]
                pos = line_info[1]
                phase_info = line_info[9]
                if abs(int(pos) - prev_pos) > 100:
                    prev_pos = int(pos)
                    if phase_info.startswith('0|1'):
                        
                        split_phase_info = phase_info.split(':')
                        dp = split_phase_info[1].split(',')
                        ref_dp = dp[0]
                        alt_dp = dp[1]
                        sample_dict[chr+":"+pos] = [chr, int(pos), int(ref_dp), int(alt_dp)]
                        
                    elif phase_info.startswith('1|0'):
                        
                        split_phase_info = phase_info.split(':')
                        dp = split_phase_info[1].split(',')
                        ref_dp = dp[1]
                        alt_dp = dp[0]
                        sample_dict[chr+":"+pos] = [chr, int(pos), int(alt_dp), int(ref_dp)]
                        
    return sample_dict

def make_file(l_dict, p_dict):

    l_sites = set(list(l_dict.keys()))
    p_sites = set(list(p_dict.keys()))

    shared_sites = l_sites.intersection(p_sites)

    for key in l_dict:
        if key in shared_sites:
            l_info = l_dict[key]
            p_info = p_dict[key]
            print(l_info[0], l_info[1], str(l_info[2]), str(l_info[3]), str(p_info[2]), str(p_info[3]))

l_info = extract_vcf_info(l_file)
p_info = extract_vcf_info(p_file)
make_file(l_info, p_info)