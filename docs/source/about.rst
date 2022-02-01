About M3K
------------

The M3K software consists of several modules: initial quality control of a sample, splitting of FASTQ files into chunks, mapping of them to reference genome(s), removal of duplicate reads and their quantification, removal of contaminated cells and determining the number of viable cells, and finally, generation of the HTML report. The program supports rerunning of the analysis pipeline from each of the modules using an input parameter. More detailed descripting is provided in further text.

.. image:: https://github.com/cbiagii/M3K-docs/raw/main/workflow.png
   :width: 900px
   :align: center



Module 1: Quality control
'''''''''''''''''''''''''
This module serves as an initial sample quality check-up and consists of two parts, the FastQC software that generates basic sample statistics and quality check-ups [REF] and an additional program that calculates nucleotide frequencies of the first 30 nucleotide positions from Read 2 file. The latter program provides insights into potential sequencing problems. For example, let’s assume that the N nucleotide occurs at the last two nucleotide positions in UMI sequences at a high frequency for some reason. This would lead to a removal of the majority of reads by M3K software, since their UMI sequences would be marked as invalid, containing more than one mismatch. Thus, reads from faulty UMIs would be ignored, and the final output matrix would consist of only a fraction of reads of the original sample size.

Module 2: FASTQ file splitting
'''''''''''''''''''''''''
To improve the software efficacy, the FASTQ files are split into a number of equally-sized chunks. The number of chunks is given with an input parameter. 

Module 3: Barcode extraction and correction
'''''''''''''''''''''''''
The barcode (BC) and Unique Molecular Identifier (UMI) sequences are extracted from each row in Read 1 file; the length of a BC is 16 nucleotides, while UMI length can be either 10 or 12 nucleotides, depending on the chemistry version used in the experiment (v.2 and >v.3, respectively). Since some of these sequences will inevitably contain sequencing errors, and the correct BCs are identified by compared to a corresponding reference list of BCs, which are then added with their corresponding UMIs to the sequence names in the Read 2. For the remaining BCs, a sequence correction of those having a hamming distance of one to any of the BCs from the reference list is attempted, while the reads from the remaining BCs containing mismatches on multiple nucleotides are excluded from further analysis. 
The implementation of the BC correction is explained in further text. First, each BC from the reference list of BCs is split into four 4-mers, and is then encoded into an 8-bit integer number. Finally, these 4-mers are added into one of the corresponding hash tables :math:`M_i`, i=1,..,4, where indexes i denote the order of 4-mers in a BC sequence (i.e., :math:`M_1` contains first 4-mers of all reference BCs, etc.), and thus, lowering the search space to only 256 elements per 4-mer. 
Next, for each BC from Read 1 file, the same procedure is applied and a counter is utilized to keep track of how many of the sequence parts are matched. If a BC has exactly three 4-mers matched, then a bitwise comparison between a BC and the rest of the BCs is performed, to calculate the number of mismatching nucleotides. In case of a mismatch at a single nucleotide position, BC correction is attempted. A BC is considered to be corrected if and only if it matches one of the reference BCs uniquely. Otherwise, a corresponding sequence of the BC in question is discarded. 

Module 4: Mapping read sequences
'''''''''''''''''''''''''
Mapping of FASTQ reads from Read 2 to a reference genome is performed independently on each chunk (see FASTQ file splitting module) using a STAR aligner software (Dobin et al., Bioinformatics, 2013; v. 2.5.2b), retaining only uniquely mapped reads.

Module 5: BAM file merging
'''''''''''''''''''''''''
Merging of chunks containing mapped sequences as previously described is performed using Samtools (Li et al. Bioinformatics, 2009; v. 1.6). 

Module 6: Deduplication procedure and quantification of read counts
'''''''''''''''''''''''''
Duplicate reads are defined as sequences from the same cell (BC), having same UMI and same or slightly different starting positions (between 1 and 3 nucleotides), as described in (Sena et al., Sci. Reports, 2018), which are excluded from the analysis. 
Quantification procedure of the remaining reads is performed by overlapping each read with coding regions from a slightly altered transcriptome (overlapping coding regions belonging to same genes are previously merged together, and their union is used instead). 
Since UMI sequences are randomly generated, there is no reference list of UMIs that can be used for their correction. Therefore, we have developed an approach to address this issue. We reason that within a cell, it is highly unlikely that two (or more) different reads whose UMI sequences have a hamming distance of one between one another, will fall within the same coding region. Thus, we assume that these reads are duplicated reads due to a sequencing error at a single nucleotide position in one (or more) UMI sequences, with differing start positions, and conclude that the corresponding reads should be removed from the sample. For each exonic region overlapping at least two reads, the algorithm performs pairwise comparisons of hamming distances from the corresponding UMIs. Pairs having a hamming distance of 1 are marked as duplicate reads and one of them is consequently excluded from further analysis. The UMI sequence comparison utilizes an approach analogous to the described approach for the BC comparison (see Barcode extraction and correction module). An overview of the used approach to quantify reads overlapping the coding regions, is illustrated in Fig. XXX.

Module 7: Sample decontamination and cell calling
'''''''''''''''''''''''''
In the final module the number of viable cells is determined from a sample. In case a sample is derived from a PDX or a CDX mouse model (patient-derived xenografts or circulating tumour cells-derived xenografts, respectively) (Hidalgo et al., Cancer Discovery, 2014, Hodgkinson et al., Nat.Med. 2014), a prior step is necessary, i.e., the identification and removal of cells belonging to a host organism (typically a mouse). In this regard, a PDX/ CDX sample is first independently mapped to both a human (hg19) and a mouse genomes (mm10), reads are quantified to produce count matrices, which are then used to calculate a mapping score for each cell, 

.. math::
   \begin{align}
   F_i =&~ \frac{H_i}{H_i + M_i}, i=1,n
   \end{align}

where :math:`F_i` represents a fraction of reads mapped to a human genome, :math:`H_i` and :math:`M_i` denote total number of quantified reads, when cells are mapped to human and mouse genomes, respectively, and n is a total number of cells in a sample.
Next, :math:`F_i` scores are ranked in an ascending order and plotted as a smoothened curve in order calculate first derivatives in each point, after which the curve is again smoothened and second derivative is calculated (Fig. XXX). The local extreme point (a turning point of the tangent) of the curve is used as a cut-off, and all cells with a score lower than a cut-off are considered to be either host or dead cells, and are thus, removed from the analysis. Alternatively, a user can repeat the procedure and set a cut-off with an input parameter. 
Finally, the procedure for determining the total number of viable cells is similar to the afore-mentioned procedure for removing unwanted host cells. In brief, cells are ranked by their number of quantified reads in descending order and plotted. Next, the curve is smoothened and first derivative calculated, representing a cut-off for viable cells (see Fig. XXX). The cell calling option is also adjustable with a parameter. 

Results: 
'''''''''''''''''''''''''
XXXXX