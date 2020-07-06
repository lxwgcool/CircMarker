**************
How to compile
**************
$ git clone https://github.com/lxwgcool/CircMarker.git
$ cd CircMarker/CircRnaDetectDraft/CircRNA_Detect_Draft/MakeFile/
$ make clean
$ make

Notice: 
Since CircMaker depends on zlib, please install zlib first (or load the zlib module first if you compile it in HPC).

The executable file is named "CircRnaDetectDraft"


*********************
How to use CircMarker
*********************
1: "CircRnaDetectDraft" and "libbamtools.so.2.4.0" need to be kept in the same folder.
2: Edit "config.ini" file to figure out the location of reads, reference, annotation file, Reads length and which chromosome you want to run. 
   The corresponding option name as below: 
   (1) Sequence Reads  : Reads1 and Reads2
   (2) Reference       : Reference
   (3) Annotation file : GTF
   (4) Reads Length    : ReadsLen
   (5) Which Chromosome: CurChromIndex (start from 1)
   Here is the example:
     
	[General]
	Reference=/home/lq/lxwg/WorkStudio/Prototype/CircRNA_Detect_Draft/Test_Data/Homo_sapiens.GRCh37.75.dna.chromosome.1.fa
	GTF=/home/lq/lxwg/WorkStudio/Prototype/CircRNA_Detect_Draft/Test_Data/chr1.gtf
	Reads1=/home/lq/lxwg/WorkStudio/Prototype/CircRNA_Detect_Draft/Test_Data/one_reads.fastq
	Reads2=
	[Parameters]
	MinSupportReads=1
	MaxSupportReads=999
	ReadsLen=101
	KmerRatio=30
	[Results]
	CurResult=
	SimulationResult=
	CIRIResult=
	CIRCexplorerResult=
	FindCircResult=
	CircBaseResult=
	[HitCompare]
	CIRIHit=
	CurHit=
	FindCircHit=
	CircExplorerHit=
	CurChromIndex=1
	[Mapping]
	BWA=

3: How to run circMarker:
   ./CircRnaDetectDraft ./config.ini

4: Result: 
  (1) ./Candidate_*.txt: Recorded the circular splicing junciton of each potential circular RNA
  (2) ./Brief_*.txt    : Breif format of "./Candidate_*.txt"
  (2) ./Brief_sum.txt  : The summary of all breif results from each chromosome.

5: Notice: 
  (1) The chromosome should be listed in gtf file sequencially, which means that the order should be 1,2,3,4,5 and etc. The chromosome should start at the first one.   

Enjoy
Best
Xin

=========
Update:
1: New feature: support circular RNA detection crouse the whole genome  --- July.24.2018
2: Delete the compiling option "-static" to fit the non root user --- July.25.2018

07/05/2020
1: Add new features for checking errors 
  (1) Check fastq file existence
  (2) Check fastq file existence
  (3) Check Annotation file (gtf) existence
2: Check sequence region names between *.fa and *.gtf files: some users might mistakenly use NCBI/EMSEMBL sequence level style (1,2, ..) for genome reference and UCSC style (chr1, chr2, ..) for gene annotation (Thanks for the suggestion of "TapscottLab"). 
3: Fix bug: out of index substr
  (1) CircMarker first sampling 8 kmers from reads to check if current reads may contribute a potential circular RNA. You may get this error if the reads is too short (e.g. after trimmomatic). This bug has been fixed in this version.
  
=========  
Q&A
Q1: Allow user to define the number of core: circMarker detects the number of available cores and use at least 25 cores. The parallelization is great. But when using shared cluster resources, core usage might be limited to some users
Comments: 
If you have enough computing resources in one computing node, you can get the best performance if the number of cores you applied equals to the number of chromosomes in your annotation (e.g. 24 cores for WES). Since circmarker uses multiple threads for parallel computing, all of the cores you requested should be in the same computing node. In another word, the cores across different HPC computing nodes cannot be used by circmarker at the same time. Actually, any multiple threads based software (such as BWA, bcl2fastq) can not use multiple cores across different nodes.

Q2: What's the format of the output and what they represent?
Comments: 
For "./Brief_sum.txt"
Col 1: Chromosme name
Col 2: CircRNA starting position
Col 3: CircRNA ending position
Col 4: Number of reads that support current circular RNA
Col 5: Strand
Col 6: Type of Circular RNA
       a) S: self-circle Case
       b) R: Regular-circle Case

  
