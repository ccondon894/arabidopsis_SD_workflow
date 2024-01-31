import sys
# Program grabs phasing info from Whatshap phased vcf files
# Extracts both leaf and pollen phasing info into one file
# Uses genetic map and linear interpolation to fill in centimorgan positions
# Output is a .haplotypes file in the following format:
# (chr, postion in bp, position in cM, somatic allele A counts, somatic allele a counts, gametic allele A counts, gametic allele a counts)


leaf = sys.argv[1]
pollen = sys.argv[2]
gen_map = sys.argv[3]

def extract_vcf_info(file):
    prev_pos = 0
    sample_dict = {}
    with open(file, 'r') as f:
        for line in f:
            if line.startswith("scaffold"):
                line_info = line.strip().split("\t")
                chr = line_info[0][9:]
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

def read_map(file):
    map = []
    scaffold = '1'
    scaffold_map = []
    with open(file, 'r') as f:
        for line in f:
            if not line.startswith("Chr"):
                split = line.strip().split()
                chr = split[0]
                bp = split[1]
                cm = split[2]
                if chr == scaffold:    
                    scaffold_map.append([chr, int(bp), float(cm)])
                else:
                    map.append(scaffold_map)
                    scaffold = chr
                    scaffold_map = []
                    scaffold_map.append([chr, int(bp), float(cm)])

        map.append(scaffold_map)

    return map


def linterpol(l0, l1, b0, b1, b):

    l = (((l1 - l0) * (b - b0)) / (b1 - b0)) + l0
    return l


def make_file(leaf, pollen, map):

    scaffold = '1'
    # i used to increment map positions 
    i = 0
    b0 = map[int(scaffold)-1][i][1]
    b1 = map[int(scaffold)-1][i+1][1]
    l0 = map[int(scaffold)-1][i][2]
    l1 = map[int(scaffold)-1][i+1][2]

    for key in leaf:
        if key in pollen:

            # check if pos is on new scaffold, update scaffold and move to next scaffold map
            if leaf[key][0] != scaffold:
                
                scaffold = leaf[key][0]
                # print("next scaffold", scaffold)
                i = 0
                b0 = map[int(scaffold)-1][i][1]
                b1 = map[int(scaffold)-1][i+1][1]
                l0 = map[int(scaffold)-1][i][2]
                l1 = map[int(scaffold)-1][i+1][2]

            # if map bounds are not on the right scaffold, increment forward until they bound the current position or we reach end of the map
            while leaf[key][1] > b1 and i < len(map[int(scaffold)-1]) - 2: 
                i += 1
                b0 = map[int(scaffold)-1][i][1]
                b1 = map[int(scaffold)-1][i+1][1]
                l0 = map[int(scaffold)-1][i][2]
                l1 = map[int(scaffold)-1][i+1][2]

            # print("bounds", b0, b1)
            # print("chromsome", leaf[key][0])
            # print("position", leaf[key][1])
            # only two scenarios left here:
            # scenario 1: the position is within the new map bounds. In this case we enter the conditional statement below
            # scenario 2: the position is still outside the bounds and we have reached the edge of our scaffold map. In this case we just skip all positions until we run off the scaffold
            if leaf[key][1] <= b1:
                # if position falls within current map bounds
                if leaf[key][1] > b0 and leaf[key][1] < b1:
                    chrom = leaf[key][0]
                    b = leaf[key][1]
                    l_A = leaf[key][2]
                    l_a = leaf[key][3]
                    p_A = pollen[key][2]
                    p_a = pollen[key][3]
                    l = linterpol(l0, l1, b0, b1, b)
                    print(chrom, b, l, l_A, l_a, p_A, p_a, sep="\t")
                # if position equals lower map bound
                elif leaf[key][1] == b0:
                    chrom = leaf[key][0]
                    b = leaf[key][1]
                    l_A = leaf[key][2]
                    l_a = leaf[key][3]
                    p_A = pollen[key][2]
                    p_a = pollen[key][3]
                    l = l0
                    print(chrom, b, l, l_A, l_a, p_A, p_a, sep="\t")
                # if position equals upper map bound
                elif leaf[key][1] == b1:
                    chrom = leaf[key][0]
                    b = leaf[key][1]
                    l_A = leaf[key][2]
                    l_a = leaf[key][3]
                    p_A = pollen[key][2]
                    p_a = pollen[key][3]
                    l = l1
                    print(chrom, b, l, l_A, l_a, p_A, p_a, sep="\t")


leaf_dict = extract_vcf_info(leaf)
pollen_dict = extract_vcf_info(pollen)
map = read_map(gen_map)
# for scaffold in map:
#     print(len(scaffold))

make_file(leaf_dict, pollen_dict, map)


                    




