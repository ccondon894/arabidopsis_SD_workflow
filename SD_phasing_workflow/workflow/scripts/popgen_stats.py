import sys
from collections import defaultdict
pi_file = sys.argv[1]

def read_file(file):
    pi_dict = defaultdict(lambda : defaultdict(float))
    with open(file, 'r') as f:
        for line in f:
            info = line.strip().split(" ")
            # grab compared samples from filename
            compar = info[0].split('.')[0]
            pi = info[1]

            if 'vs' in compar:
                compar = compar.split("_vs_")
                c1 = compar[0]
                c2 = compar[1]
                pi_dict[c1][c2] = float(pi)
            else:
                pi_dict[compar][compar] = float(pi)
    
    return pi_dict


def read_fam_file(file):
    fam_dict = {}
    with open(file, 'r') as f:
            for line in f:
                info = line.strip().split(" ")
                sample = info[0].split(".")[0]
                pi = info[1]
                fam_dict[sample] = float(pi)
    
    return fam_dict


def pi_between(pi1, pi2):

    return (pi1 * 0.5) + (pi2 * 0.5)


def calc_pi_between(pi_dict):
    pi_dict.pop('all_samples')
    fam_list = list(pi_dict.keys())
    between_dict = defaultdict(lambda : defaultdict(float))
    for i in range(len(fam_list)):
        for j in range(i, len(fam_list)):
            cross1 = fam_list[i]
            cross2 = fam_list[j]
            between_dict[cross1][cross2] = pi_between(pi_dict[cross1], pi_dict[cross2])

    return between_dict

def get_fst(pi_dict, c1, c2):
    c11_pi = pi_dict[c1][c1]
    c12_pi = pi_dict[c1][c2]
    return (c12_pi - ((c11_pi + c12_pi) / 2)) / c12_pi


def all_fst(pi_dict):

    fst_dict = defaultdict(lambda : defaultdict(float))
    for s1 in pi_dict:
        for s2 in pi_dict[s1]:
            fst = get_fst(pi_dict, s1, s2)
            fst_dict[s1][s2] = fst

    return fst_dict

def fam_fst(between_dict):
    fst_dict = defaultdict(lambda : defaultdict(float))
    for s1 in between_dict:
        for s2 in between_dict[s1]:
            pi_within = between_dict[s1][s1]
            pi_between = between_dict[s1][s2]
            fst = (pi_between - pi_within) / pi_between
            fst_dict[s1][s2] = fst

    return fst_dict

pi_dict = read_fam_file(pi_file)
print(pi_dict)

between_dict = calc_pi_between(pi_dict)
between_dict = dict(between_dict)
for key in between_dict:
    between_dict[key] = dict(between_dict[key])
# print(between_dict)

fst_dict = fam_fst(between_dict)
for s1 in fst_dict:
    for s2 in fst_dict[s1]:
        print(s1, s2, between_dict[s1][s1], between_dict[s1][s2], fst_dict[s1][s2])



# pi_dict = read_file(pi_file)
# fst_dict = all_fst(pi_dict)

# fst_dict = dict(fst_dict)
# for key in fst_dict:
#     fst_dict[key] = dict(fst_dict[key])


