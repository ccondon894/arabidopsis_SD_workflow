# Run trio phasing for vcf using Whatshap
# 1/8/23: Everything is in working order with this .smk file and conifg + Snakefile
def get_cross(wildcards):
    # print(config["crosses"][wildcards.cross])
    return config["crosses"][wildcards.cross]

# this rule splits the main vcf into family-specific smaller vcfs 
rule split_vcf:
    output:
        L_vcf = temp(config["species_folder"] + "unphased_vcfs/" + config["vcf_file"] + ".{cross}_L.vcf.gz"),
        P_vcf = temp(config["species_folder"] + "unphased_vcfs/" + config["vcf_file"] + ".{cross}_P.vcf.gz"),
    params:
        big_vcf = config["vcf_path"] + config["vcf_file"] + ".vcf.gz",
        cross = get_cross,
    conda:
        config["conda_env"],
    shell:
        '''
        vcf-subset --exclude-ref -c {params.cross}_L,{params.cross}_P01,{params.cross}_P02 {params.big_vcf} > {output.L_vcf}
        vcf-subset --exclude-ref -c {params.cross}_P,{params.cross}_P01,{params.cross}_P02 {params.big_vcf} > {output.P_vcf}
        '''        
# this rule generates the pedigree files for trio phasing
rule make_peds:
    output:
        L_ped = config["species_folder"] + "peds/{cross}_L.ped",
        P_ped = config["species_folder"] + "peds/{cross}_P.ped",
    params:
        outdir= config["species_folder"] + "peds/",
        cross=get_cross
    conda:
        config["conda_env"],
    shell:
        '''
        echo -e {params.cross}'\t'{params.cross}_P'\t'{params.cross}_P01'\t'{params.cross}_P02'\t'0'\t'1'\n' > {output.P_ped}
        echo -e {params.cross}'\t'{params.cross}_L'\t'{params.cross}_P01'\t'{params.cross}_P02'\t'0'\t'1'\n' > {output.L_ped}
        '''
# phase pollen samples
rule phase_p_sample:
    input:
        P01_bam = config["bam_folder"] + "{cross}_P01_final.bam",
        P02_bam = config["bam_folder"] + "{cross}_P02_final.bam", 
        P_bam = config["bam_folder"] + "{cross}_P_final.bam", 
        P_ped = config["species_folder"] + "peds/{cross}_P.ped",
        P_vcf = config["species_folder"] + "unphased_vcfs/" + config["vcf_file"] + ".{cross}_P.vcf.gz",
    output:
        phase_P_vcf = config["species_folder"] + "phased_vcfs/" + config["vcf_file"] + ".{cross}_P.phased.vcf",
    params:
        outdir = config["species_folder"] + "phased_vcfs/",
        ref = config["ref_genome"],
        cross = get_cross,
    conda:
        config["conda_env"], 
    shell:
        '''
        whatshap phase --ped {input.P_ped} --ref={params.ref} --recombrate 4.8 -o {output.phase_P_vcf} {input.P_vcf} {input.P_bam} {input.P01_bam} {input.P02_bam}
        '''
# phase leaf samples
rule phase_l_sample:
    input:
        P01_bam = config["bam_folder"] + "{cross}_P01_final.bam",
        P02_bam = config["bam_folder"] + "{cross}_P02_final.bam",
        L_bam = config["bam_folder"] + "{cross}_L_final.bam",
        L_ped = config["species_folder"] + "peds/{cross}_L.ped",
        L_vcf = config["species_folder"] + "unphased_vcfs/" + config["vcf_file"] + ".{cross}_L.vcf.gz",
    output:
        phase_L_vcf = config["species_folder"] + "phased_vcfs/" + config["vcf_file"] + ".{cross}_L.phased.vcf",
    params:
        outdir = config["species_folder"] + "phased_vcfs/",
        ref = config["ref_genome"],
        cross = get_cross,
    conda:
        config["conda_env"],
    shell:
        '''
        whatshap phase --ped {input.L_ped} --ref={params.ref} --recombrate 4.8 -o {output.phase_L_vcf} {input.L_vcf} {input.L_bam} {input.P01_bam} {input.P02_bam}
        '''
