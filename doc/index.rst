.. PyNAST documentation master file, created by
   sphinx-quickstart on Mon Jan 25 11:42:17 2010.

Downloading PyNAST: Latest stable release
=========================================
You can download the latest `stable release of PyNAST here <https://github.com/downloads/qiime/pynast/PyNAST-1.1.tgz>`_ and the `PyNAST OS X GUI (still PyNAST 1.0) here <https://github.com/downloads/qiime/pynast/PyNAST.app.zip>`_.

Downloading PyNAST: Development version
=======================================
If you want access to the latest-and-greatest features of PyNAST and can tolerate some instability we recommend that you check out the latest version from GitHub. You can do that with the following command: ::

    git clone git://github.com/qiime/pynast.git PyNAST

Installing PyNAST
=================
`Notes on installing and using the PyNAST command line application. <install.html>`_

`Notes on installing and using the PyNAST 1.0 Mac OS X GUI. <install_gui.html>`_

Stay up-to-date on PyNAST news
==============================
Subscribing to the PyNAST blog_ is the best way to keep up-to-date on news related to PyNAST. You can subscribe via RSS or e-mail on the front page of the blog. This is a very low traffic list, with currently around one e-mail per month or less.

The PyNAST blog is the primary means by which we will communicate information on bugs, new releases, and news to our users, so we highly recommend subscribing. We won't share subscriber information with anyone ever.

About PyNAST
============
PyNAST_ is a reimplementation of the NAST_ sequence aligner, which has become a popular tool for adding new 16s rDNA sequences to existing 16s rDNA alignments. This reimplementation is more flexible, faster, and easier to install and maintain than the original NAST implementation. PyNAST_ is built using the PyCogent Bioinformatics Toolkit.

The first versions of PyNAST (through PyNAST 1.0) were written to exactly match the results of the original NAST algorithm. Beginning with the post-PyNAST 1.0 development code, PyNAST no longer exactly matches the NAST output but is instead focused on getting better alignments. Users who wish to exactly match the results of NAST should download PyNAST 1.0.

Given a set of sequences and a template alignment, PyNAST_ will align the input sequences against the template alignment, and return a multiple sequence alignment which contains the same number of positions (or columns) as the template alignment. This facilitates the analysis of new sequences in the context of existing alignments, and additional data derived from existing alignments such as phylogenetic trees. Because any protein or nucleic acid sequences and template alignments can be provided, PyNAST_ is not limited to the analysis of 16s rDNA sequences.

PyNAST_ is presented in an open access `Bioinformatics Applications Note <http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btp636>`_.

Citing PyNAST
=============
If you make use of PyNAST_ in published work, please cite:

**PyNAST: a flexible tool for aligning sequences to a template alignment.** J. Gregory Caporaso, Kyle Bittinger, Frederic D. Bushman, Todd Z. DeSantis, Gary L. Andersen, and Rob Knight. January 15, 2010, DOI 10.1093/bioinformatics/btp636. Bioinformatics 26: 266-267.

Need help?
==========
For PyNAST_ support, you can contact `Greg Caporaso <gregcaporaso@gmail.com>`_.

.. _PyNAST: http://qiime.org/pynast
.. _blog: http://pynast.wordpress.com
.. _NAST: http://nar.oxfordjournals.org/cgi/content/full/34/suppl_2/W394
