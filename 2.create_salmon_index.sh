##standard index for gencode v44 
# selective alignment [https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/]
# passes the entire genome as the decoy for transcript quantification adjustment

cd /dssg/home/acct-dahan/share/references/gencode

grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > gencode.v44.decoys.txt
sed -i.bak -e 's/>//g' gencode.v44.decoys.txt

cat gencode.v44.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gencode.v44.gentrome.fa.gz

salmon index -t gencode.v44.gentrome.fa.gz --keepDuplicates -d gencode.v44.decoys.txt -p 12 -i salmon_index --gencode 

#gencode tx names
gunzip -c gencode.v44.transcripts.fa.gz | grep '^>' | cut -d' ' -f1 | sed 's/^>//' > gencode.v44.transcript.names.txt
# gene names
cut -d'|' -f6 gencode.v44.transcript.names.txt > gencode.v44.gene.names.txt
# tx to gene
paste -d, gencode.v44.transcript.names.txt gencode.v44.gene.names.txt > gencode.v44.tx.to.gene.csv


grep ">" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d ">" -f 2 | cut -d " " -f 1 > GRCh38.primary_assembly.genome.chrnames.txt

###gencode + TE insertions index

salmon index \
-t <(cat gencode.v44.ucsc.rmsk.salmon.fa GRCh38.primary_assembly.genome.fa) \
-i ucsc.rm/gencode.v44.repeat.aware.index \
--keepDuplicates \
-p 16 \
-d GRCh38.primary_assembly.genome.chrnames.txt


