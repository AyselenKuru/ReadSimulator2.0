# ReadSimulator2.0

How to call ReadSimulator:
For a better overview, we list the individual parameters of the Read Simulator here again, including their default values and whether they are mandatory or optional.

Mandatory Parameters
-fasta: Fasta file from which the reads are simulated.
-fidx: Index file for the Fasta file.
-gtf: Genomic annotations.
-length: Desired length of the reads.
-gene_id: Gene ID of the desired gene.
-frlength: Average fragment length from which the reads are formed.
-SD: Standard deviation for the formation of the fragment length.
-od: Directory in which the fastaq files, the read mapping info, and the additional info are placed.
Optional Parameters
-readcounts: File containing the desired genes, their transcripts, and the number of reads, default: empty file.
-transcript_count_for_gene: Number of transcripts selected for the specified gene to simulate, default: 1 transcript.
-total_read_count: Number of reads for the specified gene, default: 100 reads.

-non_unique_frq: Probability that a read is completely or partially in an overlapping region, default: 50%.
-overlap_on_more_than_1_T: Allows to generate reads that are completely or partially in an overlapping region that exists on multiple transcripts, default: false.
-minoverlapsize: Allows to generate reads that are completely or partially in an overlapping region that have a certain size, which increases the probability that the read is completely in the overlapping region, default: 0.

-dot_mutation_rate: Probability of a point mutation for bases A,C,G,T, default: 0.0, 0.0, 0.0, 0.0.
-insertion_mutation_rate: Probability of an insertion for bases A,C,G,T, default: 0.0, 0.0, 0.0, 0.0.
-deletion_mutation_rate: Probability of a deletion for bases A,C,G,T, default: 0.0, 0.0, 0.0, 0.0.
-junction_mutation_rate: Probability of a point mutation at exon boundaries, default: 0.
-junction_size: Size of the region on which the junction mutation is formed, default: 3.
-sequence_error_rate: Probability of a sequencing error, default: 0.01%.

-biological: Selects fragment positions as biologically accurate as possible.
-bio_file: File path to the biological data.
(HUMAN GTF3.8)

For CCR9   	the gene_id	is	ENSG00000173585	total_transcript_count	4
For FYCO1  	the gene_id	is	ENSG00000163820 total_transcript_count	5
For HIVEP3 	the gene_id	is	ENSG00000127124	total_transcript_count	7
For LZTFL1	the gene_id	is	ENSG00000163818	total_transcript_count	15
For	SLAMF1	the gene_id	is	ENSG00000117090	total_transcript_count	3
For	SLAMF7	the	gene_id	is	ENSG00000026751	total_transcript_count	11
