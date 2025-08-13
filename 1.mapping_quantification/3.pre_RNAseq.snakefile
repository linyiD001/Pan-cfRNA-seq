#config
tool_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/tools"

#reference
genome_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/reference/genome_index/star"
fasta_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/reference/fasta"
bed_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/reference/bed"
gtf_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/reference/gtf"

#input
data_dir="/dssg/home/acct-dahan/share/cfRNA/cleandata"
sample_ids=config["sample_id"]

#output
output_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/output"

#mapping and operate on bam file
thread_mapping=16
threads_decompress=4
map_steps=['spikein_long','univec','rRNA']
map_steps_sortbyName=['rRNA']
map_steps_rmdup=['rRNA']

#temp
temp_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/temp"

#pipeline
def get_all_inputs(wildcards):
        available_inputs = dict(
                bam=expand('{output_dir}/bam/{sample_id}/{map_step}.bam',
                        output_dir=output_dir, sample_id=sample_ids, map_step=map_steps),
                unmapped=expand('{output_dir}/unmapped/{sample_id}/{map_step}_{mate_index}.fastq.gz',
                        output_dir=output_dir, sample_id=sample_ids, map_step=map_steps, mate_index=[1, 2])
        )
        enabled_inputs = list(available_inputs.keys())
        inputs = []
        for key, l in available_inputs.items():
                if key in enabled_inputs:
                        inputs += l
        return inputs

rule all:
    input:
        get_all_inputs

map_command_pe = '''/dssg/home/acct-dahan/dahan-user3/.conda/envs/gy/bin/STAR --genomeDir {params.index} \
            --readFilesIn {input.fastq1} {input.fastq2} \
            --runThreadN {thread_mapping} \
            --outFileNamePrefix {params.output_prefix} \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx \
            --readFilesCommand gzip -d -c \
            --outSAMmultNmax 1 \
            --seedPerWindowNmax {params.seedPerWindowNmax} > {log} 2>&1
        mv {params.output_prefix}Aligned.out.bam {output.bam}
        {tool_dir}/bbmap/repair.sh in={params.output_prefix}Unmapped.out.mate1 in2={params.output_prefix}Unmapped.out.mate2 out={output.unmapped1} out2={output.unmapped2} >> {log} 2>&1
        rm -f {params.output_prefix}Unmapped.out.mate1 {params.output_prefix}Unmapped.out.mate2
        '''
rule map_spikein_long:
    input:
        fastq1 = data_dir + '/{sample_id}_1.fastq.gz',
        fastq2 = data_dir + '/{sample_id}_2.fastq.gz'
    output:
        bam = output_dir + '/bam/{sample_id}/spikein_long.bam',
        unmapped1 = output_dir + '/unmapped/{sample_id}/spikein_long_1.fastq.gz',
        unmapped2 = output_dir + '/unmapped/{sample_id}/spikein_long_2.fastq.gz'
    log:
        output_dir + '/log/{sample_id}/mapping_star/spikein_long/mapping.log'
    params:
        index = genome_dir + '/spikein_long',
        output_prefix = output_dir + '/log/{sample_id}/mapping_star/spikein_long/',
        seedPerWindowNmax = 20
    run:
        shell(map_command_pe)

rule map_univec:
    input:
        fastq1 = output_dir + '/unmapped/{sample_id}/spikein_long_1.fastq.gz',
        fastq2 = output_dir + '/unmapped/{sample_id}/spikein_long_2.fastq.gz'
    output:
        bam = output_dir + '/bam/{sample_id}/univec.bam',
        unmapped1 = output_dir + '/unmapped/{sample_id}/univec_1.fastq.gz',
        unmapped2 = output_dir + '/unmapped/{sample_id}/univec_2.fastq.gz'
    log:
        output_dir + '/log/{sample_id}/mapping_star/univec/mapping.log'
    params:
        index = genome_dir + '/univec',
        output_prefix = output_dir + '/log/{sample_id}/mapping_star/univec/',
        seedPerWindowNmax = 20
    run:
        shell(map_command_pe)

rule map_rRNA:
    input:
        fastq1 = output_dir + '/unmapped/{sample_id}/univec_1.fastq.gz',
        fastq2 = output_dir + '/unmapped/{sample_id}/univec_2.fastq.gz'
    output:
        bam = output_dir + '/bam/{sample_id}/rRNA.bam',
        unmapped1 = output_dir + '/unmapped/{sample_id}/rRNA_1.fastq.gz',
        unmapped2 = output_dir + '/unmapped/{sample_id}/rRNA_2.fastq.gz'
    log:
        output_dir + '/log/{sample_id}/mapping_star/rRNA/mapping.out'
    params:
        index = genome_dir + '/rRNA',
        output_prefix = output_dir + '/log/{sample_id}/mapping_star/rRNA/',
        seedPerWindowNmax = 20
    run:
        shell(map_command_pe)


