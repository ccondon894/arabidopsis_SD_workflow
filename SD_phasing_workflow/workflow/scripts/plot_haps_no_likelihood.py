import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import sys
# This program uses leaf and pollen haplotypes input files, and plots data per chr

hap_file = sys.argv[1]
png = sys.argv[2]

def get_hap_data(file):
    l_arr = []
    p_arr = []
    diff_arr = []
    l_scaff = []
    p_scaff = []
    diff_scaff = []
    scaffold = 1
    
    with open(file, 'r') as f:
        for line in f:
            split = line.strip().split()
            if scaffold == int(split[0]):
                pos = split[1]
                l_count_a = int(split[2])
                l_count_b = int(split[3])
                p_count_a = int(split[4])
                p_count_b = int(split[5])
                
                l_ratio = l_count_a / (l_count_a + l_count_b)
                p_ratio = p_count_a / (p_count_a + p_count_b)
                diff = p_ratio - l_ratio

                l_scaff.append([int(pos), l_ratio])
                p_scaff.append([int(pos), p_ratio])
                diff_scaff.append([int(pos), diff])

            elif scaffold != int(split[0]) and scaffold < 8:
                scaffold += 1
                l_arr.append(l_scaff)
                p_arr.append(p_scaff)
                diff_arr.append(diff_scaff)

                l_scaff = []
                p_scaff = []
                diff_scaff = []

                pos = split[1]
                l_count_a = int(split[2])
                l_count_b = int(split[3])
                p_count_a = int(split[4])
                p_count_b = int(split[5])
                
                l_ratio = l_count_a / (l_count_a + l_count_b)
                p_ratio = p_count_a / (p_count_a + p_count_b)
                diff = p_ratio - l_ratio

                l_scaff.append([int(pos), l_ratio])
                p_scaff.append([int(pos), p_ratio])
                diff_scaff.append([int(pos), diff])

        l_arr.append(l_scaff)
        p_arr.append(p_scaff)
        diff_arr.append(diff_scaff)

    return l_arr, p_arr, diff_arr                 


l_arr, p_arr, diff_arr = get_hap_data(hap_file)


def read_likelihood_file(file):
    scaffold = 1
    scaffold_arr = []
    ll_arr = []
    with open(file, 'r') as f:
        
        for line in f:
            split = line.strip().split()
            if scaffold == int(split[0]):
                pos = split[1]
                log_likelihood = split[6]
                scaffold_arr.append([int(pos), abs(float((log_likelihood)))])

            elif scaffold != int(split[0]) and scaffold < 8:
                scaffold += 1
                ll_arr.append(scaffold_arr)
                scaffold_arr = []
                
                pos = split[1]
                log_likelihood = split[6]
                scaffold_arr.append([int(pos), abs(float((log_likelihood)))])

        ll_arr.append(scaffold_arr)
    
    return ll_arr



# subsample hets to a specific window size and take the average ratio of that window
def bin_ratios(hap_arr, window_size):

    ratios = []
    positions = []
    window_ratio = 0
    i = 0
    for read in hap_arr:
        if i <= window_size:
            window_ratio += read[1]
            i += 1
        elif i > window_size:
            window_ratio = window_ratio / window_size        
            ratios.append(window_ratio)
            positions.append(read[0])
            i = 0

    return ratios, positions


def mask_aneuploidy(ratios, positions):
    filt_ratios = []
    filt_pos = []
    abs_ratios = [abs(i) for i in ratios]
    percentile = np.percentile(abs_ratios, 95)
    for i in range(len(abs_ratios)):
        if abs_ratios[i] < percentile:
            filt_ratios.append(ratios[i])
            filt_pos.append(positions[i])

    return filt_ratios, filt_pos


plt.rcParams.update({'font.size': 20})
fig, axes = plt.subplots(ncols=8, nrows=2, figsize=(60,20))

# create subplots for each chromosome and plot the data
boxes = []

for i, ax in enumerate(axes.flatten()): # type: ignore
    
    boxes.append(ax.get_position())
    if i < 8:
        L_ratios, L_positions = bin_ratios(l_arr[i], 1000)
        L_ratios, L_positions = mask_aneuploidy(L_ratios, L_positions)
        P_ratios,  P_positions = bin_ratios(p_arr[i], 1000)
        P_ratios, P_positions = mask_aneuploidy(P_ratios, P_positions)

        ax.plot(L_positions, L_ratios, color='red', label='leaf')
        ax.plot(P_positions,P_ratios, color='blue', label='pollen')
        ax.set_ylim(0.30, 0.70)
        ax.set_title("scaffold_" + str(i+1))
        ax.legend()
        ax.set_xlabel('bp')
        #ax.set_xticks([])
        #ax.set_xticklabels([])
        if i == 0:
            ax.set_ylabel("Raw Ancestry Ratios")
        else:
            ax.set_yticks([])
            ax.set_yticklabels([])

    elif i >= 8:
        diff_ratios, diff_positions = bin_ratios(diff_arr[i-8], 1000)
        diff_ratios, diff_positions = mask_aneuploidy(diff_ratios, diff_positions)
        #print(i-8)
        #print(diff_positions, diff_ratios)
        ax.plot(diff_positions, diff_ratios, color='magenta')
        ax.axhline(y=0, color='grey',linestyle='dashed', linewidth=3)
        ax.set_ylim(-0.15, 0.15)    
        ax.set_xlabel('bp')
        if i == 8:
            ax.set_ylabel("Ancestry Ratio Difference")
        else:
            ax.set_yticks([])
            ax.set_yticklabels([])

    
    
#print(boxes)
fig.subplots_adjust(wspace=0.25)
fig.suptitle(hap_file[0:3], fontsize=50)
plt.savefig(png, dpi=300)