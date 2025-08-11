name=……
read1=/dssg/home/acct-dahan/share/cfRNA/unmapped/$name/rRNA_1.fastq.gz
read2=/dssg/home/acct-dahan/share/cfRNA/unmapped/$name/rRNA_2.fastq.gz
mkdir /dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome/$name
outdir=/dssg/home/acct-dahan/share/cfRNA/RESULT/MAPPING/salmon.repeataware/salmon.repeataware_stranded_allgenome/$name

cd $outdir

set +e 

/dssg/home/acct-dahan/share/Software/miniconda3/envs/salmon1/bin/salmon quant \
-i /dssg/home/acct-dahan/share/references/gencode/RE_allgenome/ucsc.rm/gencode.v44.repeat.aware.index \
--libType ISR \
-1 $read1 \
-2 $read2 \
-p 16 \
--validateMappings \
--gcBias \
--seqBias \
--recoverOrphans \
--rangeFactorizationBins 4 \
--output $outdir \
--writeUnmappedNames || { echo "salmon quant failed"; exit 1; }
