M3K: An efficient software for pre-processing datasets from single cell RNA-seq experiments
===========================================================================================

**M3K** is a single-cell RNA-seq preprocessing tool, based on **{REFERENCE}**.

RNA velocity enables the recovery of directed dynamic information by leveraging splicing kinetics.
scVelo generalizes the concept of RNA velocity
(`La Manno et al., Nature, 2018 <https://doi.org/10.1038/s41586-018-0414-6>`_)
by relaxing previously made assumptions with a stochastic and a dynamical model that solves the full
transcriptional dynamics. It thereby adapts RNA velocity to widely varying specifications such as non-stationary populations.

scVelo is compatible with scanpy_ and hosts efficient implementations of all RNA velocity models.

M3K key applications
^^^^^^^^^^^^^^^^^^^^^^^^^
- estimate RNA velocity to study cellular dynamics.
- identify putative driver genes and regimes of regulatory changes.
- infer a latent time to reconstruct the temporal sequence of transcriptomic events.
- estimate reaction rates of transcription, splicing and degradation.
- use statistical tests, e.g., to detect different kinetics regimes.

Latest news
^^^^^^^^^^^
- Aug/2021: `Perspectives paper out in MSB <https://doi.org/10.15252/msb.202110282>`_
- Feb/2021: scVelo goes multi-core
- Dec/2020: Cover of `Nature Biotechnology <https://www.nature.com/nbt/volumes/38>`_
- Nov/2020: Talk at `Single Cell Biology <https://coursesandconferences.wellcomegenomecampus.org/our-events/single-cell-biology-2020/>`_
- Oct/2020: `Helmholtz Best Paper Award <https://twitter.com/ICBmunich/status/1318611467722199041>`_
- Oct/2020: Map cell fates with `CellRank <https://cellrank.org>`_
- Sep/2020: Talk at `Single Cell Omics <https://twitter.com/fabian_theis/status/1305621028056465412>`_
- Aug/2020: `scVelo out in Nature Biotech <https://www.helmholtz-muenchen.de/en/aktuelles/latest-news/press-information-news/article/48658/index.html>`_

References
^^^^^^^^^^
Test

Support
^^^^^^^
Found a bug or would like to see a feature implemented? Feel free to submit an
`issue <https://github.com/cbiagii/M3K-docs/issues/new/choose>`_.
Have a question or would like to start a new discussion? Head over to
`GitHub discussions <https://github.com/cbiagii/M3K-docs/discussions>`_.
In either case, you can also always send us an `email <mailto:mail@scvelo.org>`_.
Your help to improve scVelo is highly appreciated.