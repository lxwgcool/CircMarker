**************
How to compile
**************
$ git clone https://github.com/lxwgcool/CircMarker.git
$ cd CircMarker/CircRnaDetectDraft/CircRNA_Detect_Draft/MakeFile/
$ make clean
$ make

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
  (1) ./Candidate.txt: Recorded the circular splicing junciton of each potential circular RNA
  (2) ./Brief.txt    : Breif format of "./Candidate.txt"

5: Notice: 
  (1) The chromosome should be listed in gtf file sequencially, which means that the order should be 1,2,3,4,5 and etc. The chromosome should start at the first one.   

Enjoy
Best
Xin

