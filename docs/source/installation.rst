Installation
------------

M3K requires a **SLURM** environment.

GitHub
^^^^^^
To install M3K just clone the repository::

    git clone https://github.com/mnikolic/m3k


Dependencies
^^^^^^^^^^^^
We recommend to use Miniconda_.

- `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
- `STAR <https://github.com/alexdobin/STAR>`_
- `samtools <https://www.htslib.org>`_
- `R <https://www.r-project.org/>`_ - packages: *data.table*, *dplyr*, *DT*, *flexdashboard*, *ggplot2*, *htmltools*, *knitr*, *Matrix*, *R.utils*, *scales*.

It is possible to install the above programs directly from the source (see the manual of each tool) or it is also possible to install from the conda environment::

    conda install -c bioconda fastqc star samtools

There is also a `yml <https://raw.githubusercontent.com/cbiagii/M3K-docs/main/m3k.yml>`_ file with all the required dependencies above ready to be installed as a conda environment::
    
    conda env create -f m3k.yml

After installing just type `conda activate m3k` to load the environment and its dependencies.


If you run into issues, do not hesitate to approach us or raise a `GitHub issue`_.

.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _`Github issue`: https://github.com/mnikolic/m3k/issues/new/choose