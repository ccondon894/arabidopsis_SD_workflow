import matplotlib.pyplot as plt
import numpy as np
import sys
# This program uses leaf and pollen haplotypes input files, and plots data per chr
# input haplotypes file: XXX.final.haplotypes
# input likelihoods file: XXX.window_10k.likelihoods
# 9/5/23 Update: Make all likelihood Y axes scale to the largest chr Y axis

hap_file = sys.argv[1]
likelihood_file = sys.argv[2]
png = sys.argv[3]

# plt.rcParams.update({'font.size': 22})

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
                l_count_a = int(split[3])
                l_count_b = int(split[4])
                p_count_a = int(split[5])
                p_count_b = int(split[6])
                
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
                l_count_a = int(split[3])
                l_count_b = int(split[4])
                p_count_a = int(split[5])
                p_count_b = int(split[6])
                
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
    lllist = []
    with open(file, 'r') as f:
        
        for line in f:
            split = line.strip().split()
            if scaffold == int(split[0]):
                pos = split[1]
                log_likelihood = abs(float(split[6]))
                lllist.append(log_likelihood)
                scaffold_arr.append([int(pos), log_likelihood])

            elif scaffold != int(split[0]) and scaffold < 8:
                scaffold += 1
                ll_arr.append(scaffold_arr)
                scaffold_arr = []
                pos = split[1]
                log_likelihood = abs(float(split[6]))
                lllist.append(log_likelihood)
                scaffold_arr.append([int(pos), log_likelihood])

        ll_arr.append(scaffold_arr)
    
    return ll_arr, lllist

ll_arr, lllist = read_likelihood_file(likelihood_file)

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

max_ll = max(lllist)


fig, axes = plt.subplots(ncols=8, nrows=3, figsize=(15,5))

# create subplots for each chromosome and plot the data

for i, ax in enumerate(axes.flatten()): # type: ignore

    # plot ancestry ratios 
    if i < 8:
        L_ratios, L_positions = bin_ratios(l_arr[i], 1000)
        L_ratios, L_positions = mask_aneuploidy(L_ratios, L_positions)
        P_ratios,  P_positions = bin_ratios(p_arr[i], 1000)
        P_ratios, P_positions = mask_aneuploidy(P_ratios, P_positions)
        min_pos = min(L_positions)
        max_pos = max(L_positions)
        ax.plot(L_positions, L_ratios, color='red', label='leaf')
        ax.plot(P_positions,P_ratios, color='blue', label='pollen')
        ax.set_xlim(min_pos, max_pos)
        ax.set_ylim(0.30, 0.70)
        ax.tick_params(axis='y', labelsize=10)
        ax.set_yticks([0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65])
        ax.set_title("scaffold_" + str(i+1), fontsize=10)
        # ax.legend()
        ax.set_xticks([])
        #ax.set_xticks([])
        #ax.set_xticklabels([])
        if i == 0:
            ax.set_ylabel("Raw Ancestry\n Ratios", fontsize=10)
        else:
            ax.set_yticks([])
            ax.set_yticklabels([])

    # plot ancestry ratio differences 
    elif i >= 8 and i < 16:
        diff_ratios, diff_positions = bin_ratios(diff_arr[i-8], 1000)
        diff_ratios, diff_positions = mask_aneuploidy(diff_ratios, diff_positions)
        max_pos = max(diff_positions)
        min_pos = min(diff_positions)
        ax.plot(diff_positions, diff_ratios, color='magenta')
        ax.axhline(y=0, color='grey',linestyle='dashed', linewidth=3)
        ax.set_xlim(min_pos, max_pos)
        ax.set_ylim(-0.2, 0.2)
        ax.tick_params(axis='y', labelsize=10)
        ax.set_yticks([-0.15, -0.10, -0.05, 0, 0.05, 0.1, 0.15])   
        ax.set_xticks([])
        if i == 8:
            ax.set_ylabel("Ancestry Ratio\n Difference", fontsize=10)
        else:
            ax.set_yticks([])
            ax.set_yticklabels([])

    # plot likelihoods
    elif i >= 16:
        scaffold_ll = ll_arr[i-16]
        ll_ratios = [j[1] for j in scaffold_ll]
        ll_positions = [j[0] for j in scaffold_ll]
        max_pos = max(ll_positions)
        # pos_25 = round(np.percentile(ll_positions, 25)
        # pos_50 = np.percentile(ll_positions, 50)
        # pos_75 = np.percentile(ll_positions, 75)
        ax.plot(ll_positions, ll_ratios, color='red')
        ax.set_xlim([0, max_pos])
        ax.set_ylim([0, max_ll + (max_ll * 0.1)])
        if max_pos < 20000000:
            ax.set_xticks([5000000, 10000000, 15000000])
        elif max_pos < 25000000:
            ax.set_xticks([5000000, 10000000, 15000000, 20000000])
        else:
            ax.set_xticks([5000000, 10000000, 15000000, 20000000, 25000000])
        ax.set_xlabel('bp', fontsize=10)
        ax.tick_params(axis='both', labelsize=10)
        ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
        ax.xaxis.get_offset_text().set_fontsize(10)
        if i == 16:
            ax.set_ylabel("Log Likelihood\n Ratio", fontsize=10)
        else:
            ax.set_yticks([])
            ax.set_yticklabels([])
            
    
# fig.subplots_adjust(wspace=0.25)
fig.suptitle(hap_file[0:3], fontsize=13)
fig.subplots_adjust(hspace=0.1, wspace=0)
fig.align_ylabels()

plt.savefig(png, dpi=600)