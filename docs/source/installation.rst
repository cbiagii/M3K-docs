Installation
------------

M3K requires ... We recommend to use ...

PyPI
^^^^

#!/bin/bash

# TODO GTF files need to be re-generated to be on the safe side (Carlos) TODO

### Paths (for execution on Cluster - to be setup by the user)
export fastqc="/projects/cangen/milos/sc/testing/m3k/sw/FastQC/fastqc"
export star="/projects/cangen/milos/sc/testing/m3k/sw/STAR-2.5.2b/bin/Linux_x86_64/STAR"
export hg19_dir="/projects/cangen/milos/sc/testing/m3k/genome/HomoSapiensGRCh37_STAR"
export mm10_dir="/projects/cangen/milos/sc/testing/m3k/genome/mm10_STAR"
export samtools="/projects/cangen/milos/sc/testing/m3k/sw/samtools-1.2/samtools"

### Locally used file paths - no need to be changed by the user
export bc_v2="barcodes-737K-august-2016.txt"
export bc_v3="barcodes-3M-february-2018.txt"
export gtf_hg19_mn="gtf/hg37.75_exons_utrs_final_sorted2.gtf"
export gtf_mm10_mn="gtf/mm38.93_exons_uniq_merged_final.gtf"
export version="0.341"
export gzip="gzip.sh"

Install scVelo from PyPI_ using::

    pip install -U scvelo

``-U`` is short for ``--upgrade``.
If you get a ``Permission denied`` error, use ``pip install -U scvelo --user`` instead.


Development Version
^^^^^^^^^^^^^^^^^^^

To work with the latest development version, install from GitHub_ using::

    pip install git+https://github.com/theislab/scvelo@develop

or::

    git clone https://github.com/theislab/scvelo && cd scvelo
    git checkout --track origin/develop
    pip install -e .

``-e`` is short for ``--editable`` and links the package to the original cloned
location such that pulled changes are also reflected in the environment.

To contribute to scVelo, ``cd`` into the cloned directory and
install the latest packages required for development together with the pre-commit hooks::

    pip install -r requirements-dev.txt
    pre-commit install


Dependencies
^^^^^^^^^^^^

- `anndata <https://anndata.readthedocs.io/>`_ - annotated data object.
- `scanpy <https://scanpy.readthedocs.io/>`_ - toolkit for single-cell analysis.
- `numpy <https://docs.scipy.org/>`_, `scipy <https://docs.scipy.org/>`_, `pandas <https://pandas.pydata.org/>`_, `scikit-learn <https://scikit-learn.org/>`_, `matplotlib <https://matplotlib.org/>`_.


Parts of scVelo (directed PAGA and Louvain modularity) require (optional)::

    pip install python-igraph louvain


Using fast neighbor search via `hnswlib <https://github.com/nmslib/hnswlib>`_ further requires (optional)::

    pip install pybind11 hnswlib


If you run into issues, do not hesitate to approach us or raise a `GitHub issue`_.

.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _PyPI: https://pypi.org/project/scvelo
.. _Github: https://github.com/theislab/scvelo
.. _`Github issue`: https://github.com/theislab/scvelo/issues/new/choose