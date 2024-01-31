# remove mapping bias by remapping only forward reads and matching read length distributions in leaf-pollen pairs. 
# 1/10/23: Everything is working in the dry run. I still haven't done a real run but I expect the python scripts should all work as intended.

def get_cross(wildcards):
    # print(config["crosses"][wildcards.cross])
    return config["crosses"][wildcards.cross]

def get_scaffold(wildcards):
    # print(config["scaffolds"][wildcards.scaffold])
    return config["scaffolds"][wildcards.scaffold]

# trim only forward reads since we are mapping only those
rule trim_r1:
    input:
        lzip = config["fastq_folder"] + "{cross}_L_R1.fastq.gz",
        pzip = config["fastq_folder"] + "{cross}_P_R1.fastq.gz",
        
    output:
        lzip = config["fastq_folder"] + "{cross}_L_R1.retrimmed.fastq.gz", # trimmed + zipped leaf
        pzip = config["fastq_folder"] + "{cross}_P_R1.retrimmed.fastq.gz", # trimmed + zipped pollen
    params:
        cross = get_cross,
        path = config["fastq_folder"],
        lunzip =config["fastq_folder"] + "{cross}_L_R1.retrimmed.fastq", # trimmed + unzipped leaf
        punzip = config["fastq_folder"] + "{cross}_P_R1.retrimmed.fastq" # trimmed + unzipped pollen
    conda:
        config["conda_env"],
    shell: # python script trim_reads_v4.py automatically determines which file needs to be trimmed
        '''
        python workflow/scripts/trim_reads_v4.py {input.pzip} {input.lzip} {params.cross}_P_R1 {params.cross}_L_R1 {params.path} && 
        if [ -f {params.lunzip} ]
        then
            cp {input.pzip} {output.pzip} &&
            pigz {params.lunzip} &&
        elif [ -f {params.punzip} ]
        then
            cp {input.lzip} {output.lzip} &&
            pigz {params.punzip} &&
        fi
        '''

rule map_leaf: # remap retrimmed leaf reads 
    output:
        bam = config["species_folder"] + "retrim_bams/" + "{cross}_L.retrim.sorted.bam",
        bai = config["species_folder"] + "retrim_bams/" + "{cross}_L.retrim.sorted.bam.bai",
    params:
        cross = get_cross,
        ref = config["ref_genome"],
        fq = config["fastq_folder"] + "{cross}_L_R1.retrimmed.fastq",
        fqz = config["fastq_folder"] + "{cross}_L_R1.retrimmed.fastq.gz",
    threads:
        config['threads'],
    conda:
        config["conda_env"],
    shell:
        '''
        if test -f "{params.fqz}"
        then
            gzip -d {params.fqz}
        fi
        bwa mem -t {threads} {params.ref} {params.fq} | samtools sort -o {output.bam}
        samtools index {output.bam}
        pigz {params.fq}
        '''

rule map_pollen: # remap retrimmed pollen reads
    output:
        bam = "retrim_bams/{cross}_P.Ahal_ref.retrim.sorted.bam",
        bai = "retrim_bams/{cross}_P.Ahal_ref.retrim.sorted.bam.bai",
    params:
        cross = get_cross,
        ref = config["ref_genome"],
        fq = config["fastq_folder"] + "{cross}_P_R1.retrimmed.fastq",
        fqz = config["fastq_folder"] + "{cross}_P_R1.retrimmed.fastq.gz",
    threads: 
        config['threads'],
    conda:
        config["conda_env"],
    shell:
        '''
        if test -f "{params.fqz}"
        then
            gzip -d {params.fqz}
        fi
        bwa mem -t {threads} {params.ref} {params.fq} | samtools sort -o {output.bam}
        samtools index {output.bam}
        pigz {params.fq}
        '''

rule call: # call variants per scaffold for all produced bam files. Include per sample genotype info. Only call for the 8 main scaffolds.
    input:
        ref = config["ref_genome"],
        p_bam = expand(config["species_folder"] + "retrim_bams/{cross}_P.Ahal_ref.retrim.sorted.bam",cross=config["crosses"]),
        p_bai = expand(config["species_folder"] + "retrim_bams/{cross}_P.Ahal_ref.retrim.sorted.bam.bai",cross=config["crosses"]),
        l_bam = expand(config["species_folder"] + "retrim_bams/{cross}_L.Ahal_ref.retrim.sorted.bam",cross=config["crosses"]),
        l_bai = expand(config["species_folder"] + "retrim_bams/{cross}_L.Ahal_ref.retrim.sorted.bam.bai",cross=config["crosses"])
    output:
        bcf = config["species_folder"] + "mpileup/{scaffold}.retrim_calls.bcf"
    params:
        scaffold = get_scaffold,
    conda:
        config["conda_env"],
    shell:
        '''
        bcftools mpileup -a 'FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP' -Ou -r {params.scaffold} -f {input.ref} {input.p_bam} {input.l_bam} | bcftools call -mv -Ob -o {output.bcf}
        '''

rule index: # index bcf files for filtering.
    input:
        bcf = config["species_folder"] + "mpileup/{scaffold}.Ahal_ref.retrim_calls.bcf"
    output:
        csi = config["species_folder"] + "mpileup/{scaffold}.Ahal_ref.retrim_calls.bcf.csi"
    params:
        scaffold = get_scaffold,
    conda:
        config["conda_env"],
    shell:
        '''
        bcftools index {input.bcf}
        '''

rule filter: # quality filter and remove indels/multiallelic sites from bcf, and convert to vcf.
    input:
        bcf = config["species_folder"] + "mpileup/{scaffold}.Ahal_ref.retrim_calls.bcf"
    output:
        bcf = config["species_folder"] + "mpileup/{scaffold}.Ahal_ref.retrim_calls.final.vcf.gz"
    params:
        scaffold = get_scaffold,
    conda:
        config["conda_env"],
    shell:
        '''
        bcftools view {input.bcf} -i '%QUAL>=20' --max-alleles 2 --exclude-types indels -Oz -o {output.bcf}
        '''



