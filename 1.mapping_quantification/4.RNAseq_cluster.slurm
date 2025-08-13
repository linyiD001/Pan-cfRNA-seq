#!/bin/bash

#SBATCH --job-name=test_cluster
#SBATCH --partition=debug64c512g
#SBATCH -N 1
#SBATCH --mail-type=end
#SBATCH --mail-user=bulubulu88@sjtu.edu.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --ntasks-per-node=8

for i in ~/cfRNA/code/0_pre/sample_ids/* ;do

sample_id=$(basename $i)

echo "$sample_id start at `date`"
snakemake --rerun-incomplete --nolock \
        -s pre_RNAseq.snakefile -j 2 --latency-wait 60 \
        --configfile $i \
        --cluster-config RNAseq_cluster.json \
        --cluster "sbatch -p {cluster.partition} -N {cluster.nodes} -n {cluster.ntasks} --job-name {cluster.jobname} --partition {cluster.partition} --ntasks-per-node=64 --output {cluster.output} --error {cluster.error}" > cluster_log/0.RNAseq_${sample_id}.log 2> cluster_log/0.RNAseq_${sample_id}.err
echo "$sample_id end at `date`"

done
