# anotomat
Annotates genes based on sequence coverage, takes start position as input and gives annotated files (dna, protein, gff).  
Use option "-h" or "--help" to get further instructions.  


  
## MANUAL  


	anotomat.py
	
	
		--genome <genome.fasta>
		--pos <start_positions.txt> *exact format see below*
		--cov <coverage> *file with '#chr bp coverage' from "samtools depth -a -Q 0 <file.bam>"*
		--mincov <int> minimal coverage of the START pos for a reannotation*default=0*
		--mincov_exon <int> minimal coverage for annotation of exon after intron*default=0*
		--gt use less stringency in defining introns *G._.G defines intron*(default=GT_AG)
		--cores <int> number of cores used *default=(cpu_count-1)*
		--name <gene name> *default='genes'*
		--out <name for output files> *default='genes'*

## EXAMPLE

	EXAMPLE for start_positions.txt:

	Bgt_chr-01	49207	~BgtAcSP-31373	#make a comment
	Bgt_chr-01	284818	~BgtE-20066	#tab separated???
	Bgt_chr-01	300256	~BgtE-20114	#special sign "~ "before gene-name!!!
                






