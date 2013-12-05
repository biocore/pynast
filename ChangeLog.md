pynast ChangeLog
================

pynast 1.2.2
------------

Bug fix release, released on December 5th 2013

* Fixed bug where some temporary files were (still) being written to /tmp.
* PyNAST is now pip-installable.

PyNAST 1.2.1
------------

Bug fix release, released on November 14th 2013

* PyNAST is now distributed under the terms of the Modified BSD License (previously GPL).
* Improved handling of temporary files. Previously these were written to /tmp and the user didn't have any control over this. Temporary files are now written to tempfile.gettempdir() by default, which is determined by python on a platform-specific basis, and users can override this value by passing the temp_dir option to blast_align_unaligned_seqs, pynast_seqs, ipynast_seqs and the pynast command line interface.

PyNAST 1.2
----------

Feature release, released on November 9th 2012

* Required PyCogent version is now 1.5.3.
* If muscle is installed, required version is now 3.8.31.

PyNAST 1.1
----------

Feature release, released on March 31st 2010

* Switch from BLAST for database search to uclust for database search. BLAST is no longer available for database searching.
* Switch from BLAST as the default pairwise alignment to uclust as the default pairwise aligner.
* Addition of setup.py to facilitate installation.

PyNAST 1.0
----------

Initial release, released on January 25th 2010

