Application in real data
------------------------

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