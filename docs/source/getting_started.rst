Getting Started
---------------

Here, you will be briefly guided through the basics of how to use M3K.
Once you are set, the following tutorials go straight into analysis of RNA velocity,
latent time, driver identification and many more.

First of all, the input data for scVelo are two count matrices of pre-mature (unspliced) and mature (spliced) abundances,
which can be obtained from standard sequencing protocols, using the `velocyto`_ or `loompy/kallisto`_
counting pipeline.



function help_msg() {
    echo -e "\nM3K Software v. $version (Peifer Lab)"
    echo -e "\nUsage:"
    echo -e "\nSingle samples: run_m3k.sh <input_parameters>"
    echo -e "Multiple samples: run_many.sh <input.conf>" 

    echo -e "\nInput parameters:" 
    echo -e "\t-f1\tFilename 1"
    echo -e "\t-f2\tFilename 2"
    echo -e "\t-o\tPath to output folder" 
    echo -e "\t-m\tModule [0]"
    echo -e "\t-u\tUMI length [10]"
    echo -e "\t-sp1\t Smooting parameter [0.5]"
    echo -e "\t-sp2\t Smooting parameter [0.5]"
    echo -e "\t-sp3\t Smooting parameter [0.5]"
    echo -e "\t-t\t Number of threads [12]"
    echo -e "\t-c\t Number of chunks [24]\n"
    
    echo -e "Available modules:"
    echo -e "\t0 - Whole pipeline and QC reports"
    echo -e "\t1 - Split fastq files"
    echo -e "\t2 - Extract and correct barcodes"
    echo -e "\t3 - Map reads to reference genome(s)"
    echo -e "\t4 - Merge BAM files"
    echo -e "\t5 - Deduplicate and quantify reads"
    echo -e "\t6 - Merge output files"
    echo -e "\t7 - Identify viable cells & remove host cells\n"
    
    echo -e "\tExample of the <input.conf> file is provided in the M3K online manual."
    echo -e "\n"
    exit 1
}


- input.conf:
### DO NOT LEAVE ANY EMPTY ROWS ### (command example: S03334_T_FF_01_1.fastq.gz S03334_T_FF_01_2.fastq.gz /projects/cangen/milos/sc/Cleidson/S03334_T_FF_01_test 0 10 hg 0.5 0.5 0.5 8 24)
#
### S H A L L O W   S E Q. ### 
# Dec. 21
#PEC_JB_1.fastq.gz PEC_JB_2.fastq.gz /projects/cangen/milos/sc/shallow_seq/PEC_JB 0 12 hg 0.5 0.5 0.5 12 24 10000
#S03516_T_FF_41_1.fastq.gz S03516_T_FF_41_2.fastq.gz /projects/cangen/milos/sc/shallow_seq/S03516_T_FF_41 7 12 hg 0.5 0.5 0.5 12 24 10000 6000
#S03856_1.fastq.gz S03856_2.fastq.gz /projects/cangen/milos/sc/shallow_seq/S03856 2 12 hg 0.5 0.5 0.5 12 24 10000


scVelo workflow at a glance
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Import scvelo as::

    import scvelo as scv

For beautified visualization you can change the matplotlib settings to our defaults with::

    scv.set_figure_params()

Read your data
''''''''''''''
Read your data file (loom, h5ad, csv, ...) using::

    adata = scv.read(filename, cache=True)

which stores the data matrix (``adata.X``),
annotation of cells / observations (``adata.obs``) and genes / variables (``adata.var``), unstructured annotation such
as graphs (``adata.uns``) and additional data layers where spliced and unspliced counts are stored (``adata.layers``) .

.. raw:: html

    <img src="http://falexwolf.de/img/scanpy/anndata.svg" style="width: 300px">

If you already have an existing preprocessed adata object you can simply merge the spliced/unspliced counts via::

    ldata = scv.read(filename.loom, cache=True)
    adata = scv.utils.merge(adata, ldata)

If you do not have a datasets yet, you can still play around using one of the in-built datasets, e.g.::

    adata = scv.datasets.pancreas()

The typical workflow consists of subsequent calls of preprocessing (``scv.pp.*``), analysis tools (``scv.tl.*``) and plotting (``scv.pl.*``).

Basic preprocessing
'''''''''''''''''''
After basic preprocessing (gene selection and normalization),
we compute the first- and second-order moments (means and uncentered variances) for velocity estimation::

    scv.pp.filter_and_normalize(adata, **params)
    scv.pp.moments(adata, **params)

Velocity Tools
''''''''''''''
The core of the software is the efficient and robust estimation of velocities, obtained with::

    scv.tl.velocity(adata, mode='stochastic', **params)

The velocities are vectors in gene expression space obtained by solving a stochastic model of transcriptional dynamics.
The solution to the deterministic model is obtained by setting ``mode='deterministic'``.

The solution to the dynamical model is obtained by setting ``mode='dynamical'``, which requires to run
``scv.tl.recover_dynamics(adata, **params)`` beforehand.

The velocities are stored in ``adata.layers`` just like the count matrices.

The velocities are projected into a lower-dimensional embedding by translating them into likely cell transitions.
That is, for each velocity vector we find the likely cell transitions that are in accordance with that direction.
The probabilities of one cell transitioning into another cell are computed using cosine correlation
(between the potential cell transition and the velocity vector) and are stored in a matrix denoted as velocity graph::

    scv.tl.velocity_graph(adata, **params)

Visualization
'''''''''''''

Finally, the velocities can be projected and visualized in any embedding (e.g. UMAP) on single cell level, as gridlines, or as streamlines::

    scv.pl.velocity_embedding(adata, basis='umap', **params)
    scv.pl.velocity_embedding_grid(adata, basis='umap', **params)
    scv.pl.velocity_embedding_stream(adata, basis='umap', **params)

For every tool module there is a plotting counterpart, which allows you to examine your results in detail, e.g.::

    scv.pl.velocity(adata, var_names=['gene_A', 'gene_B'], **params)
    scv.pl.velocity_graph(adata, **params)


.. _`velocyto`: http://velocyto.org/velocyto.py/tutorial/cli.html
.. _`loompy/kallisto`: https://linnarssonlab.org/loompy/kallisto/index.html