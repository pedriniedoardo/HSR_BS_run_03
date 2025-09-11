# rule test:
#     '''
#     This is a test rule
#     '''
#     input:
#         test = config['test']
#     conda:
#         'bioinfo'
#     output:
#         test = 'results/test/test.txt'
#     log:
#         'logs/test/test.log'
#     benchmark:
#         'benchmarks/test/test.txt'
#     resources:
#         mem_mb = 500,
#         cpus = 1
#     threads: 1
#     params:
#         annotations = config['REGION1']
#     shell:
#         '''
#         echo "run the test" > {log}
#         cat {input.test} > {output.test}
#         zcat {params.annotations} | awk 'NR<=10' >> {output.test}
#         echo "test done" >> {log}
#         '''

rule unzip_barcodes:
    '''
    This rule is to unzip the barcodes file in the cellranger output folder. the conversion is adding the unzipped version of the file in the folder, keeping the zipped one. 
    '''
    input:
        barcodes = "../../data/cellranger61/out/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    conda:
        'cellsnp_vireo'
    output:
        barcodes = "../../data/cellranger61/out/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv"
    log:
        'logs/{sample}/00_barcode_unzipping.log'
    benchmark:
        'benchmarks/{sample}/00_barcode_unzipping.txt'
    resources:
        mem_mb = lambda wildcards, attempt: 500 * attempt,
        # mem_mb = 500,
        cpus = 1
    threads: 1
    params:
    shell:
        '''
        echo "unzipping <{wildcards.sample}>" > {log}
        gunzip < {input.barcodes} > {output.barcodes} 2>> {log}
        echo "completed unzipping of <{wildcards.sample}>" >> {log}
        '''

# implementation version 01
# rule cellsnp_mode1a:
#     '''
#     This is the rule to run cell snp lite in mode 1a (for droplet based dataset).
#     '''
#     input:
#         barcodes = rules.unzip_barcodes.output.barcodes,
#         bam = "data/cellranger/{sample}/outs/possorted_genome_bam.bam"
#     conda:
#         'cellsnp_vireo'
#     output:
#         base_vcf = "results/01_cellsnplite/{sample}/cellSNP.base.vcf",
#         samples_tsv = "results/01_cellsnplite/{sample}/cellSNP.samples.tsv",
#         AD_mtx = "results/01_cellsnplite/{sample}/cellSNP.tag.AD.mtx",
#         DP_mtx = "results/01_cellsnplite/{sample}/cellSNP.tag.DP.mtx",
#         OTH_mtx = "results/01_cellsnplite/{sample}/cellSNP.tag.OTH.mtx"
#     log:
#         'logs/{sample}/02_cellsnp_mode1a.log'
#     benchmark:
#         'benchmarks/{sample}/02_cellsnp_mode1a.txt'
#     resources:
#         mem_gb = 16,
#         cpus = 1
#     threads: 1
#     params:
#         wd = 'results/01_cellsnplite/{sample}',
#         vcf = config['REGION1'],
#         p = config['cellsnplite_p'],
#         minCOUNT = config['cellsnplite_minCOUNT']
#     shell:
#         '''
#         echo "start cellsnp-lite for <{wildcards.sample}>" >> {log}
#         cellsnp-lite -p {params.p} -s {input.bam} -O {params.wd} -R {params.vcf} -b {input.barcodes} --minCOUNT {params.minCOUNT}
#         echo "completed run cellsnp-lite for <{wildcards.sample}>" >> {log}
#         '''

rule cellsnp_mode1a:
    '''
    This is the rule to run cell snp lite in mode 1a (for droplet based dataset).
    '''
    input:
        barcodes = rules.unzip_barcodes.output.barcodes,
        bam = "../../data/cellranger61/out/{sample}/outs/possorted_genome_bam.bam"
    conda:
        'cellsnp_vireo'
    output:
        folder = directory("results/01_cellsnplite/{sample}")
    log:
        'logs/{sample}/01_cellsnp_mode1a.log'
    benchmark:
        'benchmarks/{sample}/01_cellsnp_mode1a.txt'
    resources:
        mem_mb = lambda wildcards, attempt: config['cellsnp_RAM'] * attempt,
        # mem_mb = 16000,
        cpus = config['cellsnp_CPU']
    threads: config['cellsnp_CPU']
    params:
        vcf = config['REGION'],
        p = config['cellsnp_CPU'],
        minCOUNT = config['cellsnplite_minCOUNT']
    shell:
        '''
        echo "create the output folder for <{wildcards.sample}>" > {log}
        mkdir -p {output.folder} 2>> {log}
        echo "folder created" >> {log}
        echo "start cellsnp-lite for <{wildcards.sample}>" >> {log}
        cellsnp-lite -p {params.p} -s {input.bam} -O {output.folder} -R {params.vcf} -b {input.barcodes} --minCOUNT {params.minCOUNT} --gzip
        echo "completed run cellsnp-lite for <{wildcards.sample}>" >> {log}
        '''

rule vireo_cellsnp_mode1a:
    '''
    This is the rule to run vireo on the ouput of cellsnp lite.
    '''
    input:
        folder_snplite = rules.cellsnp_mode1a.output.folder
    conda:
        'cellsnp_vireo'
    output:
        folder = directory("results/02_vireo/{sample}")
    log:
        'logs/{sample}/02_vireo.log'
    benchmark:
        'benchmarks/{sample}/02_vireo.txt'
    resources:
        mem_mb = lambda wildcards, attempt: config['vireo_RAM'] * attempt,
        # mem_mb = 16000,
        cpus = config['vireo_CPU']
    threads: config['vireo_CPU']
    params:
        donor_N = lambda w: config["SAMPLES"]["{}".format(w.sample)]['donor_N'],
        p = config['vireo_CPU'],
        random_seed = config['random_seed']
    shell:
        '''
        echo "create the output folder for <{wildcards.sample}>" > {log}
        mkdir -p {output.folder} 2>> {log}
        echo "folder created" >> {log}
        echo "start vireo for <{wildcards.sample}>" >> {log}
        vireo -p {params.p} -c {input.folder_snplite} -N {params.donor_N} -o {output.folder} --randSeed {params.random_seed}
        echo "completed run vireo for <{wildcards.sample}>" >> {log}
        '''