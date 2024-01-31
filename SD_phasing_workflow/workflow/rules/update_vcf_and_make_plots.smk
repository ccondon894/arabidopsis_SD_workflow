# Update phased vcfs with remapped depths and do basic filtering 

def get_cross(wildcards):
    # print(config["crosses"][wildcards.cross])
    return config["crosses"][wildcards.cross]

rule filter_vcf: # filters the vcf of sites that do not match genotype or are not identically phased.
    input:
        l = config["species_folder"] + "phased_vcfs/" + config["vcf_file"] + ".{cross}_L.phased.vcf",
        p = config["species_folder"] + "phased_vcfs/" + config["vcf_file"] + ".{cross}_P.phased.vcf",
    output:
        l = temp(config["species_folder"] + "phased_vcfs/" + config["vcf_file"] + ".{cross}_L.filtered.vcf"),
        p = temp(config["species_folder"] + "phased_vcfs/" + config["vcf_file"] + ".{cross}_P.filtered.vcf"),
    params:
        cross = get_cross,
        logfile = config["vcf_file"] + "{cross}.filtering.log",
    conda:
        config["conda_env"],
    shell:
        '''
        python workflow/scripts/filter_vcf_with_stats.py {input.l} {input.p} {params.cross}_L {params.cross}_P {output.l} {output.p} 2>logs/{params.logfile}
        '''

rule update_depths_l: # Updates leaf allele depths inside the filtered vcfs to the ones obtained after trimming and remapping reads.
    input:
        l = config["species_folder"] + "phased_vcfs/" + config["vcf_file"] + ".{cross}_L.filtered.vcf",
    output:
        l = config["species_folder"] + "phased_vcfs/" + config["vcf_file"] + ".{cross}_L.updated.vcf",
    params:
        cross = get_cross,
        mp_path = config["species_folder"] + "mpileup/",
    conda:
        config["conda_env"],
    shell:
        '''
        python workflow/scripts/update_filt_vcf_depths.py {input.l} {params.cross}_L {params.mp_path} > {output.l}
        '''

rule update_depths_p: # Updates pollen allele depths inside the filtered vcfs to the ones obtained after trimming and remapping reads.
    input:
        p = config["species_folder"] + "phased_vcfs/" + config["vcf_file"] + ".{cross}_P.filtered.vcf",
    output:
        p = config["species_folder"] + "phased_vcfs/" + config["vcf_file"] + ".{cross}_P.updated.vcf",
    params:
        cross = get_cross,
        mp_path = config["species_folder"] + "mpileup/",
    conda:
        config["conda_env"],
    shell:
        '''
        python workflow/scripts/update_filt_vcf_depths.py {input.p} {params.cross}_P > {output.p}
        '''


rule make_haps_file: # makes haplotypes file for ancestry ratio plotting
    input:
        l = config["species_folder"] + "phased_vcfs/" + config["vcf_file"] + ".{cross}_L.updated.vcf",
        p = config["species_folder"] + "phased_vcfs/" + config["vcf_file"] + ".{cross}_P.updated.vcf",
    output:
        haps = config["species_folder"] + "haplotypes/{cross}.final.haplotypes",
    params:
        cross = get_cross,
    conda:
        config["conda_env"],
    shell:
        '''
        python workflow/scripts/extract_haps_v2.py {input.l} {input.p} > {output.haps}
        '''

rule map_sd: # runs the max likelihood estimation framework MAP_SD. Window size is 10k sites.
    input: 
        haps = config["species_folder"] + "haplotypes/{cross}.final.haplotypes",
    output:
        likelihoods = config["species_folder"] + "likelihoods/{cross}.window_10K.likelihoods",
    params:
        cross = get_cross,
    conda:
        config["conda_env"],
    shell:
        '''
        python workflow/MAP_SD/src/map_sd -d {input.haps} -w 10000 > {output.likelihoods} &
        '''

rule plot_ratios: # make ratio plots
    input:
        haps = config["species_folder"] + "haplotypes/{cross}.final.haplotypes",
        likelihoods = config["species_folder"] + "likelihoods/{cross}.window_10K.likelihoods",
    output:
        png = config["species_folder"] + "ratio_plots/{cross}.ratios.png",
    params:
        cross = get_cross,
    conda:
        config["conda_env"],
    shell:
        '''
        python workflow/scripts/plot_haps_v3.py {input.haps} {input.likelihoods} {output.png}
        '''

