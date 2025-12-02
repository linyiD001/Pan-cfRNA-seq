#config
tool_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/tools"

#reference
genome_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/reference/genome_index/star"
fasta_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/reference/fasta"
bed_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/reference/bed"
gtf_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/reference/gtf"

#input
sample_ids=config["sample_id"]

#output
output_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/output"

#temp
temp_dir="/dssg/home/acct-dahan/dahan-user3/cfRNA/temp"


#mapping and operate on bam file
thread_mapping=16
threads_decompress=4
map_steps=['hg38']
map_steps_sortbyName=['hg38']
map_steps_rmdup=['hg38']
#count parameter
strandness="reverse" #{"no","forward","reverse"}
paired_end="True" #{"True","False"}
count_multimap_reads="True" #{"True","False"}
min_mapping_quality=0
count_overlapping_features="True" #{"True","False"}
count_levels=['hg38','hg38_rmdup']

######################################################################
######################################################################

#pipeline
def get_all_inputs(wildcards):
	available_inputs = dict(
		bam=expand('{output_dir}/bam/{sample_id}/{map_step}.bam',
			output_dir=output_dir, sample_id=sample_ids, map_step=map_steps),
		bam_sort_by_name=expand('{output_dir}/bam/{sample_id}/{map_step_sortbyName}_sortbyName.bam',
			output_dir=output_dir, sample_id=sample_ids, map_step_sortbyName=map_steps_sortbyName),
		bam_rmdup=expand('{output_dir}/bam/{sample_id}/{map_step_rmdup}_rmdup.bam',
			output_dir=output_dir, sample_id=sample_ids, map_step_rmdup=map_steps_rmdup),
		unmapped=expand('{output_dir}/unmapped/{sample_id}/{map_step}_{mate_index}.fastq.gz',
			output_dir=output_dir, sample_id=sample_ids, map_step=map_steps, mate_index=[1, 2]),
		count=expand('{output_dir}/counts/hg38/{sample_id}/featurecount',
			output_dir=output_dir, sample_id=sample_ids),
		count_rmdup=expand('{output_dir}/counts/hg38_rmdup/{sample_id}/featurecount',
			output_dir=output_dir, sample_id=sample_ids),
		samtools_stats=expand('{output_dir}/log/{sample_id}/samtool_stats/{map_step}.txt',
			output_dir=output_dir, sample_id=sample_ids, map_step=map_steps),
		samtools_stats_rmdup=expand('{output_dir}/log/{sample_id}/samtool_stats/{map_step_rmdup}_rmdup.txt',
			output_dir=output_dir, sample_id=sample_ids, map_step_rmdup=map_steps_rmdup),
		bam_read_pair=expand('{output_dir}/log/{sample_id}/read_pair/{map_step}.txt',
			output_dir=output_dir, sample_id=sample_ids, map_step=map_steps),
		bam_rmdup_read_pair=expand('{output_dir}/log/{sample_id}/read_pair/{map_step_rmdup}_rmdup.txt',
			output_dir=output_dir, sample_id=sample_ids, map_step_rmdup=map_steps_rmdup)

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

map_command_pe = '''/dssg/home/acct-dahan/share/Software/miniconda3/envs/cfRNA_processing/bin/STAR --genomeDir {params.index} \
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
map_command_circ = '''/dssg/home/acct-dahan/dahan-user3/.conda/envs/gy/bin/STAR --genomeDir {params.index} \
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

rule map_hg38:
	input:
		fastq1='{output_dir}/unmapped/{sample_id}/rRNA_1.fastq.gz',
		fastq2='{output_dir}/unmapped/{sample_id}/rRNA_2.fastq.gz'
	output:
		bam='{output_dir}/bam/{sample_id}/hg38.bam',
		unmapped1='{output_dir}/unmapped/{sample_id}/hg38_1.fastq.gz',
		unmapped2='{output_dir}/unmapped/{sample_id}/hg38_2.fastq.gz'
	log:
		'{output_dir}/log/{sample_id}/mapping_star/hg38/mapping.log'
	params:
		index=genome_dir+'/hg38',
		output_prefix='{output_dir}/log/{sample_id}/mapping_star/hg38/',
		seedPerWindowNmax=50
	run:
		shell(map_command_pe)

rule sort_bam_by_name:
	input:
		'{output_dir}/bam/{sample_id}/{map_step_sortbyName}.bam'
	output:
		'{output_dir}/bam/{sample_id}/{map_step_sortbyName}_sortbyName.bam'
	shell:
		'''
		/dssg/home/acct-dahan/share/Software/miniconda3/envs/cfRNA_processing/bin/samtools sort -n -T {temp_dir} -o {output} {input}
		'''

rule samtools_stats:
	input:
		'{output_dir}/bam/{sample_id}/{map_step}.bam'
	output:
		'{output_dir}/log/{sample_id}/samtool_stats/{map_step}.txt'
	wildcard_constraints:
		map_step='(hg38)|(circRNA)'
	shell:
		'''
		/dssg/home/acct-dahan/share/Software/miniconda3/envs/cfRNA_processing/bin/samtools stats {input} > {output}
		'''

rule samtools_stats_rmdup:
	input:
		'{output_dir}/bam/{sample_id}/{map_step_rmdup}_rmdup.bam'
	output:
		'{output_dir}/log/{sample_id}/samtool_stats/{map_step_rmdup}_rmdup.txt'
	shell:
		'''
		/dssg/home/acct-dahan/share/Software/miniconda3/envs/cfRNA_processing/bin/samtools stats {input} > {output}
		'''


rule remove_duplicates:
	input:
		'{output_dir}/bam/{sample_id}/{map_step_sortbyName}_sortbyName.bam'
	output:
		bam_rmdup='{output_dir}/bam/{sample_id}/{map_step_sortbyName}_rmdup.bam',
		metrics='{output_dir}/log/{sample_id}/remove_duplicate/{map_step_sortbyName}_rmdup'
	log:
		'{output_dir}/log/{sample_id}/picard_rmdup/{map_step_sortbyName}_remove_duplicates.log'
	shell:
		'''
		java -jar {tool_dir}/picard/picard.jar \
		MarkDuplicates REMOVE_DUPLICATES=true \
		ASSUME_SORT_ORDER=queryname \
		I={input} \
		O={output.bam_rmdup} \
		M={output.metrics} \
		READ_NAME_REGEX=null > {log} 2>&1
		'''

rule featurecounts:
	input:
		bam='{output_dir}/bam/{sample_id}/hg38_sortbyName.bam',
		gtf=gtf_dir+'/Homo_sapiens.GRCh38_ERCC.gtf'
	output:
		counts='{output_dir}/counts/hg38/{sample_id}/featurecount',
		summary='{output_dir}/counts/hg38/{sample_id}/featurecount.summary'
	params:
		strandness={'no': 0, 'forward': 1, 'reverse': 2}[strandness],
		paired_end={'True': '-p', 'False': ''}[paired_end],
		min_mapping_quality={min_mapping_quality},
		count_multimap_reads={'True': '-M','False': ''}[count_multimap_reads],
		count_overlapping_features={'True': '-O','False': ''}[count_overlapping_features]
	log:
		'{output_dir}/log/{sample_id}/featurecount/hg38.log'
	shell:
		'''
		/dssg/home/acct-dahan/share/Software/miniconda3/envs/cfRNA_processing/bin/featureCounts {params.count_overlapping_features} -t exon -g gene_id {params.count_multimap_reads} \
		-s {params.strandness} -Q {params.min_mapping_quality} \
		{params.paired_end} -a {input.gtf} -o {output.counts} {input.bam} > {log} 2>&1
		'''

rule featurecounts_rmdup:
	input:
		bam='{output_dir}/bam/{sample_id}/hg38_rmdup.bam',
		gtf=gtf_dir+'/Homo_sapiens.GRCh38_ERCC.gtf'
	output:
		counts='{output_dir}/counts/hg38_rmdup/{sample_id}/featurecount',
		summary='{output_dir}/counts/hg38_rmdup/{sample_id}/featurecount.summary'
	params:
		strandness={'no': 0, 'forward': 1, 'reverse': 2}[strandness],
		paired_end={'True': '-p', 'False': ''}[paired_end],
		min_mapping_quality={min_mapping_quality},
		count_multimap_reads={'True': '-M','False': ''}[count_multimap_reads],
		count_overlapping_features={'True': '-O','False': ''}[count_overlapping_features]
	log:
		'{output_dir}/log/{sample_id}/featurecount/hg38_rmdup.log'
	shell:
		'''
		/dssg/home/acct-dahan/share/Software/miniconda3/envs/cfRNA_processing/bin/featureCounts {params.count_overlapping_features} -t exon -g gene_id {params.count_multimap_reads} \
		-s {params.strandness} -Q {params.min_mapping_quality} \
		{params.paired_end} -a {input.gtf} -o {output.counts} {input.bam} > {log} 2>&1
		'''

rule featurecounts_circ_rmdup:
	input:
		bam='{output_dir}/bam/{sample_id}/circRNA_rmdup.bam',
		gtf=gtf_dir+'/lhsa_hg19_circRNA.txt'
	output:
		counts='{output_dir}/counts/circ_rmdup/{sample_id}/featurecount',
		summary='{output_dir}/counts/circ_rmdup/{sample_id}/featurecount.summary'
	params:
		strandness={'no': 0, 'forward': 1, 'reverse': 2}[strandness],
		paired_end={'True': '-p', 'False': ''}[paired_end],
		min_mapping_quality={min_mapping_quality},
		count_multimap_reads={'True': '-M','False': ''}[count_multimap_reads],
		count_overlapping_features={'True': '-O','False': ''}[count_overlapping_features]
	log:
		'{output_dir}/log/{sample_id}/featurecount/circRNA_rmdup.log'
	shell:
		'''
		featureCounts {params.count_overlapping_features} -t exon -g gene_id {params.count_multimap_reads} \
		-s {params.strandness} -Q {params.min_mapping_quality} \
		{params.paired_end} -a {input.gtf} -o {output.counts} {input.bam} > {log} 2>&1
		'''



rule count_bam_read_pairs:
	input:
		'{output_dir}/log/{sample_id}/samtool_stats/{map_step}.txt'
	output:
		'{output_dir}/log/{sample_id}/read_pair/{map_step}.txt'
	wildcard_constraints:
		map_step='(hg38)|(circRNA)'
	shell:
		'''
		awk 'BEGIN{{OFS="\t";FS="\t"}}/^SN/{{if($2 == "reads mapped and paired:") print int($3/2)}}' {input} > {output}
		'''

rule count_bam_rmdup_read_pairs:
	input:
		'{output_dir}/log/{sample_id}/samtool_stats/{map_step_rmdup}_rmdup.txt'
	output:
		'{output_dir}/log/{sample_id}/read_pair/{map_step_rmdup}_rmdup.txt'
	shell:
		'''
		awk 'BEGIN{{OFS="\t";FS="\t"}}/^SN/{{if($2 == "reads mapped and paired:") print int($3/2)}}' {input} > {output}
		'''

