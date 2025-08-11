# ucsc genome browser
## hgdownload.soe.ucsc.edu:goldenPath/hg38/database
# path: /dssg/home/acct-dahan/share/references/ucsc.rm
## all the reference tables are stored here, details on access are in the `my_sql` directory scripts
## using MySQL queries to generate:
### 1. rmsk GTF
### 2. rmsk BED --> used to generate FASTA with `bedtools getfasta`
### 3. rmsk tx2gene
### 4. rmsk 'info'; contains the class and family details for each TE 'gene'
### 1-4 failed on HPC, unable to connect to ucsc server for unknown reason(try jcloud instead._succeed)
### 5. rmsk FA. 

##1
# wraps a MySQL query to generate the UCSC RMSK trancript (insertion) to gene (repName) CSV for `Salmon` gene level summarization

# generates analysis style tx to gene for good parsing of RNAseq results
# following column format (from GENCODE) of:
# enst, ensg, tx, gene, len, biotype

## host: genome-mysql.cse.ucsc.edu
## -N = --skip-column-names
## -A = --no-auto-rehash
## -D = --database

# sed converts to CSV; I think this approach allows for more flexible execution and control over STDOUT
####
mysql --user=genome --host=genome-mysql.soe.ucsc.edu -N -A -D hg38 -e \
	'select concat(repName, 
		       "_range=", genoName, ":", genoStart, "-", genoEnd, 
		       "_strand=", strand),
		concat(repName, 
		       "_range=", genoName, ":", genoStart, "-", genoEnd, 
		       "_strand=", strand),
		repName,
		ABS(genoEnd - genoStart) AS len,
		repClass
		from rmsk' | sed 's/\t/,/g' > rmsk.insert.csv


##2
# wraps a MySQL query to generate the UCSC RMSK bed required for bedtools getfasta and further analyses

## host: genome-mysql.cse.ucsc.edu
## -N = --skip-column-names
## -A = --no-auto-rehash
## -D = --database
## --batch = print as TSV, escape special characters

# awk:
## split the genoName column by '_' and store in object `a`
## print genoName column = if [ a[2] is empty ]; then print a[1]; else print a[2] .. print rest of columns as is
### this finds lines with chrX_IDv1_type naming and extracts the middle string
# needed to name the alt/patch/fix chromosomes according to gencode style

###
mysql --batch --user=genome --host=genome-mysql.soe.ucsc.edu -N -A -D hg38 -e \
	'select genoName,
		genoStart,
		genoEnd,
		concat(repName, 
		       "_range=", genoName, ":", genoStart, "-", genoEnd, 
		       "_strand=", strand),
		swScore,
		strand
		from rmsk' | awk -F '\t' -v OFS='\t' '{split($1, a, /_/); print (a[2] == "") ? a[1] : a[2], $2, $3, $4, $5, $6}' | sed 's/v1/.1/; s/v2/.2/1' > rmsk.insert.bed

###3
# wraps a MySQL query to generate the UCSC RMSK GTF required for tximport and further analyses
# tx ID naming matches reference fasta naming convention

## host: genome-mysql.cse.ucsc.edu
## -N = --skip-column-names
## -A = --no-auto-rehash
## -D = --database
## --batch = print as TSV, escape special characters

# awk:
## split the genoName column by '_' and store in object `a`
## print genoName column = if [ a[2] is empty ]; then print a[1]; else print a[2] .. print rest of columns as is
### this finds lines with chrX_IDv1_type naming and extracts the middle string
# needed to name the alt/patch/fix chromosomes according to gencode style


mysql --batch --user=genome --host=genome-mysql.cse.ucsc.edu -N -A -D hg38 -e \
	'select genoName,
		"hg38_rmsk",
		"exon",
		genoStart,
		genoEnd,
		swScore,
		strand,
		".",
		concat("gene_id", " ", "\"", repName, "\"", ";", " ", 
			"transcript_id", " ", "\"", repName, 
			"_range=", genoName, ":", genoStart, "-", genoEnd, 
			"_strand=", strand, "\"", ";")
		from rmsk' \
		| awk -F '\t' -v OFS='\t' '{split($1, a, /_/); print (a[2] == "") ? a[1] : a[2], $2, $3, $4, $5, $6, $7, $8, $9, $10}' | sed 's/v1/.1/; s/v2/.2/1'	>rmsk.insert.gtf
		
		
###4
# wraps a MySQL query to generate the UCSC RMSK 'info' table that describes TE name, class, and family

## host: genome-mysql.cse.ucsc.edu
## -N = --skip-column-names
## -A = --no-auto-rehash

## -D = --database
# `distinct` reduces the selection to distinct row entries
# this table should be used in a left join where x is the data table and y is this info table

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -A -D hg38 -e \
	'select distinct repName,
		repClass,
		repFamily
		from rmsk' >rmsk.insert.txt



####5
#!/bin/bash
# wraps a MySQL query to generate the UCSC RMSK trancript (insertion) to gene (repName) CSV for `Salmon` gene level summarization

## host: genome-mysql.cse.ucsc.edu
## -N = --skip-column-names
## -A = --no-auto-rehash
## -D = --database

# sed converts to CSV; I think this approach allows for more flexible execution and control over STDOUT
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -N -A -D hg38 -e \
	'select concat(repName, 
		       "_range=", genoName, ":", genoStart, "-", genoEnd, 
		       "_strand=", strand),
		repName
		from rmsk' | sed 's/\t/,/g' > rmsk.insert.short.csv
	
####5
##Usage:   bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf>
#Options: 
#        -fi     Input FASTA file(try GRCh38.primary_assembly.genome.fa.gz(some chromosome (KI GL etc) was not found in FASTA)))
#        -fo     Output file (opt., default is STDOUT
#        -bed    BED/GFF/VCF file of ranges to extract from -fi
#        -tab    Write output in TAB delimited format.
#-name	Use the name field and coordinates for the FASTA header 
#-nameOnly	Use the name field for the FASTA header


samtools faidx /dssg/home/acct-dahan/share/references/gencode/GRCh38.primary_assembly.genome.fa
#bedtools getfasta  -fi /dssg/home/acct-dahan/share/references/gencode/GRCh38.primary_assembly.genome.fa -bed ucsc.rm/rmsk.insert.bed -nameOnly -fo ucsc.rm/rmsk.insert.fa 

bedtools getfasta -fi /dssg/home/acct-dahan/share/references/gencode/GRCh38.primary_assembly.genome.fa -bed ucsc.rm/rmsk.insert.bed -s -nameOnly | sed 's/(+)*$//; s/(-)*$//g' >  ucsc.rm/rmsk.insert.fa


awk '{gsub(/\047/, "\"", $0); print}' \
/dssg/home/acct-dahan/share/references/gencode/RE_allgenome/ucsc.rm/rmsk.insert.gtf \
> /dssg/home/acct-dahan/share/references/gencode/RE_allgenome/ucsc.rm/rmsk.insert.corrected.pre.gtf

cat -t /dssg/home/acct-dahan/share/references/gencode/RE_allgenome/ucsc.rm/rmsk.insert.corrected.pre.gtf |head

awk -v OFS="\t" '{print $1, $2, $3, $4+1, $5, $6, $7, $8, substr($0, index($0,$9))}' \
/dssg/home/acct-dahan/share/references/gencode/RE_allgenome/ucsc.rm/rmsk.insert.corrected.pre.gtf \
> /dssg/home/acct-dahan/share/references/gencode/RE_allgenome/ucsc.rm/rmsk.insert.corrected.gtf

gunzip -c "/dssg/home/acct-dahan/share/references/gencode/gencode.v44.primary_assembly.annotation.gtf.gz" > "/dssg/home/acct-dahan/share/references/gencode/RE_allgenome/gencode.v44.ucsc.rmsk.salmon.gtf"
cat "/dssg/home/acct-dahan/share/references/gencode/RE_allgenome/ucsc.rm/rmsk.insert.corrected.gtf">> "/dssg/home/acct-dahan/share/references/gencode/RE_allgenome/gencode.v44.ucsc.rmsk.salmon.gtf"

cat <(gunzip -c /dssg/home/acct-dahan/share/references/gencode/gencode.v44.transcripts.fa.gz) RE_allgenome/ucsc.rm/rmsk.insert.fa > /dssg/home/acct-dahan/share/references/gencode/RE_allgenome/gencode.v44.ucsc.rmsk.salmon.fa

cat /dssg/home/acct-dahan/share/references/gencode/gencode.v44.tx.to.gene.csv /dssg/home/acct-dahan/share/references/gencode/RE_allgenome/ucsc.rm/rmsk.insert.short.csv > /dssg/home/acct-dahan/share/references/gencode/RE_allgenome/gencode.v44.ucsc.rmsk.tx.to.gene.csv


###
#cat /dssg/home/acct-dahan/share/references/gencode/gencode.v44.tx.to.gene.csv /dssg/home/acct-dahan/share/references/gencode/ucsc.rm/rmsk.insert.short.csv > /dssg/home/acct-dahan/share/references/gencode/gencode.v44.ucsc.rmsk.tx.to.gene.csv
