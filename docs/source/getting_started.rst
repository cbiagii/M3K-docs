Getting Started
---------------

Here, you will be briefly guided through the basics of how to use M3K.

First of all, the input data for M3K are ...

Generating genome indexes
^^^^^^^^^^^^^^^^^^^^^^^^^
To generate the reference indexes it is necessary to have the **fasta sequence** and the **annotation GTF** file. We recommend this to be the first step because M3K depends on index to do the alignment step. 

M3K can be used in two main cases:
1) Data contains only *Homo sapiens* sequences;
2) The data contains a mixture of *Homo sapiens* and *Mus musculus* sequences.

So, it is necessary to create the indexes according to the user's needs.

For this, the **STAR** software is used to create the index as can be seen in the example below::

    STAR --runThreadN {threads.number} --runMode genomeGenerate --genomeDir {path.to.indexDir} 
    --genomeFastaFiles {fasta.reference} --sjdbGTFfile {gtf.reference} --sjdbOverhang 100 
    --outFileNamePrefix {index.prefix}

We provide in this `link <>`_ ready-made references for *Homo sapiens* and *Mus musculus*.


Modifying *config_paths.sh* file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
After cloning the M3K repository and accessing this folder, there is a file called *config_paths.sh*. We will need to edit this file according to the user's specifications. The file looks like this below::

    #!/bin/bash

    ### Paths (for execution on Cluster - to be setup by the user)
    export fastqc="/path/to/fastqc"
    export star="/path/to/STAR"
    export samtools="/path/to/samtools"
    export hg19_dir="/path/to/hsaIndex"
    export mm10_dir="/path/to/mmIndex"
    export gtf_hg19_mn="path/to/hsa_gtf_reference"
    export gtf_mm10_mn="path/to/mm_gtf_reference"

    ### Locally used file paths - no need to be changed by the user
    export bc_v2="barcodes-737K-august-2016.txt"
    export bc_v3="barcodes-3M-february-2018.txt"
    export version="0.341"
    export gzip="gzip.sh"

To start, the user needs to change the paths of the *fastqc*, *star* and *samtools* softwares (hint: type ``which fastq`` which will return the path).

Then it will be necessary to modify the reference paths for *Homo sapiens* and *Mus musculus*. The specified path will be the same used in the ``--genomeDir`` parameter plus the ``--outFileNamePrefix`` parameter (for example: if I used the ``--genomeDir`` parameter as **/path/to/index** and the ``--outFileNamePrefix`` parameter as **hg19**, I will replace it in the file *config_paths.sh* the variable *hg19_dir* by **/path/to/index/gh19**). The same logic mentioned above must be considered if it is necessary to create an index for *Mus musculus*.

Finally it will be necessary to modify the variables *gtf_hg19_mn* and *gtf_mm10_mn* according to the path of the GTF file used to create the index (file passed to the ``--sjdbGTFfile`` parameter).


Available modules on command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
There are 8 modules (0 to 8) available in M3K as shown below:
- **0**: Whole pipeline and QC reports;
- **1**: Split fastq files;
- **2**: Extract and correct barcodes;
- **3**: Map reads to reference genome(s);
- **4**: Merge BAM files;
- **5**: Deduplicate and quantify reads;
- **6**: Merge output files;
- **7**: Identify viable cells & remove host cells.


Running M3K for a single sample
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The parameters available to run M3K for a single sample are as follows:

- ``-f1`` or ``--file1``: Path to read 1 file
- ``-f2`` or ``--file2``: Path to read 2 file
- ``-o`` or ``--outpath``: Path to output folder
- ``-m`` or ``--module``: Module to run [default: 0]
- ``-u`` or ``--umi``: UMI length. 10 or 12 nts (depending on 10x chemistry). [default: 10]
- ``g`` or ``--genome``: Genome (hg, mm or pdx)
- ``-sp1`` or ``--smoothing1``: Smooting parameter [default: 0.5]
- ``-sp2`` or ``smoothing2``: Smooting parameter [default: 0.5]
- ``-sp3`` or ``smoothing3``: Smooting parameter [default: 0.5]
- ``-t`` or ``--threads``: Number of threads [default: 12]
- ``-c`` or ``--chunks``: Number of chunks to split fastq file [default: 24]
- ``ct`` or ``--cutoff``: [default: 10000]
- ``vct`` or ``--vcutoff``: Viable number of cells [default: 0]
- ``cct`` or ``--ccutoff``: Contaminated cells [default: 0]

Once this is done, just type::

    sh run_m3k.sh --file1 /path/to/sample_R1.fastq.gz --file2 /path/to/sample_R2.fastq.gz 
    --outpath /path/to/output --module [number of desired module] --umi [10 or 12 nts] 
    --genome hg --smoothing1 0.5 --smoothing2 0.5 --smoothing3 0.5 --threads [number of threads] 
    --chuncks [number of chuncks] --cutoff 10000 --vcutoff 0 ccutoff 0


Running M3K for multiple samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To run multiple samples just create the ``input.conf`` file in the following format::

    sample1_R1.fastq.gz sample1_R2.fastq.gz /path/to/output_sample1 0 10 hg 0.5 0.5 0.5 8 24 10000 2500
    sample2_R1.fastq.gz sample2_R2.fastq.gz /path/to/output_sample2 0 12 mm 0.5 0.5 0.5 8 36 20000 5000
    sample3_R1.fastq.gz sample3_R2.fastq.gz /path/to/output_sample3 0 12 pdx 0.5 0.5 0.5 8 36 20000 3000

Comments can be made using the ``#`` character. Do not leave any empty lines in the above file.

After that, just type the following command::
    
    run_many.sh /path/to/input.conf