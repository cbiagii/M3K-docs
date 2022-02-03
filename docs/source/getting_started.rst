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

- ``-f1``: Path to read 1 file
- ``-f2``: Path to read 2 file
- ``-o``: Path to output folder
- ``-m``: Module to run [default: 0]
- ``-u``: UMI length. Choose 10 if using v2 from 10x or 12 if using v3. [default: 10]
- ``-sp1``: Smooting parameter [default: 0.5]
- ``-sp2``: Smooting parameter [default: 0.5]
- ``-sp3``: Smooting parameter [default: 0.5]
- ``-t``: Number of threads [default: 12]
- ``-c``: Number of chunks to split fastq file [default: 24]

Once this is done, just type::
    sh run_m3k.sh -f1 /path/to/sample_R1.fastq.gz -f2 /path/to/sample_R2.fastq.gz 
    -o /path/to/output -m [number of desired module] -u [v2 uses 10 and v3 uses 12] 
    -sp1 0.5 -sp2 0.5 -sp3 0.5 -t [number of threads] -c [number of chuncks]


Running M3K for multiple samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
run_many.sh <input.conf>

- input.conf:
### DO NOT LEAVE ANY EMPTY ROWS ### (command example: S03334_T_FF_01_1.fastq.gz S03334_T_FF_01_2.fastq.gz /projects/cangen/milos/sc/Cleidson/S03334_T_FF_01_test 0 10 hg 0.5 0.5 0.5 8 24)
#
### S H A L L O W   S E Q. ### 
# Dec. 21
#PEC_JB_1.fastq.gz PEC_JB_2.fastq.gz /projects/cangen/milos/sc/shallow_seq/PEC_JB 0 12 hg 0.5 0.5 0.5 12 24 10000
#S03516_T_FF_41_1.fastq.gz S03516_T_FF_41_2.fastq.gz /projects/cangen/milos/sc/shallow_seq/S03516_T_FF_41 7 12 hg 0.5 0.5 0.5 12 24 10000 6000
#S03856_1.fastq.gz S03856_2.fastq.gz /projects/cangen/milos/sc/shallow_seq/S03856 2 12 hg 0.5 0.5 0.5 12 24 10000