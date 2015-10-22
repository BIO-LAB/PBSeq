#README for PBSeq#
**Li Zhang (leo.zhang@nuaa.edu.cn)**

##Introduction#
PBSeq is a software for estimating gene and isoform expression levels from RNA-seq data. The PBSeq package provides an user-friendly interface and supports threads for parallel computation.
##Installation#
To compile PBSeq, simply run in the PBSeq folder.

`
$ bash setup.sh
`

Requirements:

* PBSeq software uses Python (v.2.7) to pre-process the RNA-seq data and C language to calculate the gene and isoform expression levels.
* In PBSeq, the Python codes use two special modules, [NumPy](http://www.numpy.org/) and [PP](http://www.parallelpython.com/) (parallel python).
* PBSeq uses Bowtie to align sequencing reads to transcriptreference sequences, so you must have [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) installed.

## Usage#

### Step 1. Aligning Sequencing Reads

Now, create index and align:

`
$ bowtie-build -f ensGene.fasta ensGene.ref_transcript.index
`
    
`
$ bowtie -t -f -p 4 -a -m 100 --suppress 2,6,7,8 ensGene.ref_transcript.index raw_data.fasta align_reads.output
`
	
	
	

Notice:

* PBSeq need use the transcript reference sequences, which can be downloaded from [UCSC](http://genome.ucsc.edu/cgi-bin/hgTables?command=start) and [Ensembl](http://asia.ensembl.org/index.html) websites. eg: ensGene.fasta, refGene.fasta, knownGene.fasta from [UCSC](http://genome.ucsc.edu/cgi-bin/hgTables?command=start) website and xxx.cdna.all.fa from [Ensembl](http://asia.ensembl.org/index.html) website.
* In Bowtie, '--suppress 2,6,7,8' must be set for the downstream analysis of PBSeq.
* If you want to get more usage information for Bowtie, please click [here](http://bowtie-bio.sourceforge.net/index.shtml) to visit the Bowtie website.
 

###Step 2. Pre-processing Annotation Inforamtion

PBSeq need to pre-process the annotation file and the corresponding reference sequences for the following steps.

`
$ python ./PBSeq/preprocessAnnotation.py -- Type ensGene --AnnotationFile ensGene.txt --SequenceFile ensGene.ref_transcript.fasta --OutputName ensGene
`
	

Options:

* -t/--Type: The type of annotation corresponds the transcript reference sequence. eg. refGene.fa-> refGene and xxx.cdna.all.fa -> Ensembl.
* -a/--AnnotationFile: the annotation includes the gene and isoform information corresponding the transcript reference sequence.
* -s/--SequenceFile: the corresponding transcript reference sequences.
* -o/--OutputName: The 'OutputName' is the header name of annotation files, eg: ensGene.Gene.Info. This file includes the names of gene and isoforms.

Notice:

* Currently, PBSeq supports four types of annotations, which are ensGene, refGene, knownGene and Ensembl.
* For the types of refGene, ensGene and knownGene, the reference sequences and annotation files(not gtf) are downloaded from UCSC website.
* For the type of Ensembl, the annotation files is xxx.gtf.
 

###Step 3. Pre-processing Alignment Files

Now PBSeq need to pre-process the alignment files for the following steps, which include pre-computing the probabilities for each alignment files, calculating the gene bias and extracting the count data of each genes.

`
$ python ./PBSeq/preprocessAlignment.py --Type ensGene --AligmentFiles align_reads1.output --Annotation ensGene
`

or

`
$ python ./PBSeq/preprocessAlignment.py --Type ensGene--AligmentFiles align_reads1.output --Annotation ensGene --InputGene target.GeneName
`

Options:

* -t/--Type: The type of annotation corresponds the transcript reference sequence. eg. refGene.fa-> refGene and xxx.cdna.all.fa -> Ensembl.
* -d/--AlignmentFiles: input the alignment files. eg: data1.output.
* -a/--Annotation: only need input the header name of annotation file, which is same as the '--OutputName' in the Step2.
* -i/--InputGene: optional parameters. If you only calculate a part of genes, you can choose this option. If you ignore this option, PBSeq will calculate all genes.
 

###Step 4. Calculating Expression Values

Now PBSeq starts calculating expression values.

`
$ python ./PBSeq/calculateExpression.py --Annotation ensGene --log
`

or

`
$ python ./PBSeq/calculateExpression.py --Annotation ensGene --InputGene target.GeneName --log
`

Options:

* -a/--Annotation: only need input the header name of annotation file, which is same as the '--OutputName' in the Step2.
* -i/--InputGene: optional parameters. If you only calculate a part of genes, you can choose this option. If you ignore this option, PGSeq will calculate all genes in the annotation file.
* -l/--log: optional parameters. If you want to detect the differential expressed genes or isoforms, you need choose this option. This option will calculate the logged values of expression for each genes and isoforms.

Output Files:

The PBSeq will produce four output files, which include gene.mean, gene.standard.deviation, isoform.mean and isoform.standard.deviation.

If you choose the '-l/--log', the additional four output files, which includes gene.mean.log, gene.standard.deviation.log, isoform.mean.log and isoform.standard.deviation.log.

Description of output files:

* ***gene/isoform.mean:*** these files contain the gene or isoforme expression. The first column is name of gene of isoform and the remaining columns are gene or isoform expression for each alignment files.
* ***gene/isoform.standard.deviation:*** these files contains the standard deviation of gene or isoform expression. The first column is name of gene of isoform and the remaining columns are the standard deviation of gene or isoform expression for each alignment files.
* ***gene/isoform.mean.log:*** same as gene/isoform.mean files.
* ***gene/isoform.standard.deviation.log:*** same as gene/isoform.standard.deviation files.
 

#Example

Here, we give two simple examples of PBSeq. 

##Example 1:

* Annotation file: refGene.chr1.txt
* Reference sequence file: refGene.chr1.fasta
* Sequencing reads: raw_data.fasta
When you use PBSeq, you need the above files at least. We suppose that the annotation file and the reference sequence file both are downloaded from UCSC website.

`
$ bowtie-build -f refGene.chr1.fasta refGene.chr1.ref_transcript
`

`
$ bowtie -t -f -p 4 -a -m 100 --suppress 2,6,7,8 refGene.ref_transcript raw_data.fasta ReadAligment.output
`

`
$ python ./PBSeq/preprocessAnnotation.py -t refGene -a refGene.chr1.txt -s refGene.chr1.fasta -o refGene.chr1
`

`
$ python ./PBSeq/preprocessAlignment.py -t refGene -a refGene.chr1 -d ReadAligment.output
`

`
$ python ./PBSeq/calculateExpression.py -a refGene.chr1 --log
`

##Example 2:

* Annotation file: Homo_sapiens.GRCH37.67.chr1.gtf
* Reference sequence file: Homo_sapiens.GRCH37.67.chr1.fasta
* Sequencing reads: raw_data.fasta
* Target Genes: test.Gene.Info (200 genes)
When you use PBSeq, you need the above files at least. We suppose that the annotation file and the reference sequence file both are downloaded from UCSC website.

`
$bowtie-build -f Homo_sapiens.GRCH37.67.chr1.fasta Ensembl.chr1.ref_transcript
`

`
$bowtie -t -f -p 4 -a -m 100 --suppress 2,6,7,8 Ensembl.chr1.ref_transcript raw_data.fasta ReadAligment.output
`

`
$pyhton ./PBSeq/preprocessAnnotation.py -t Ensembl -a Homo_sapiens.GRCH37.67.chr1.gtf -s Homo_sapiens.GRCH37.67.chr1.fasta -o Ensembl.chr1
`

`	
$python ./PBSeq/preprocessAlignment.py -t Ensembl -a Ensembl.chr1 -dReadAligment.output -i test.Gene.Info
`

`	
$python ./PBSeq/calculateExpression.py -a Ensembl.chr1 --log -i test.Gene.Info
`



#Differential Expression Analysis

Popular differentail expression analysis tools such as DESeq and edgeR do not consider the variance of gene/isoform expression and can not detect the differentially expressed isoforms.

We have proposed a BDSeq method to detect differentially expressed genes/isoforms. BDSeq method not only uses the gene/isoform expression level, but also considers the measurement uncertainty of expression level. 

PBSeq method is able to obtain the expression level and corresponding meausrement uncertainty, then use BDSeq method for differnetianly expression analysis.

 

#Authors

The PBSeq algorithm is developed by Li Zhang and Xuejun Liu. The PBSeq software is mainly implemented by Li Zhang.

 
