.. install_:

*************************************************************
Installing and using the PyNAST command line application
*************************************************************

Downloading PyNAST
==================
THIS DOCUMENTATION REFERS ONLY TO THE SVN VERSION OF PYNAST. SEE http://pynast.sourceforge.net FOR THE 1.0 DOCUMENTATION.

You can download the latest development branch of PyNAST here with the command:

	::
	 
		svn co https://pynast.svn.sourceforge.net/svnroot/pynast PyNAST

Required software
=================
PyNAST_ is built on the PyCogent_ package, and uses NCBI's BLAST_ software. You must have PyCogent_ 1.4.0 and BLAST_ 2.2.22 installed to run PyNAST_.

Optional software
=================
If you'd like to perform pairwise alignments using MUSCLE_, MAFFT_, or ClustalW_, you must have those programs installed on your machine and in your system path. Currently tested versions are MUSCLE_ v3.6, MAFFT v6.602b, and ClustalW 1.81.

Installation steps
==================
#. Download PyCogent_ 1.4.0 and its dependencies, Python_ 2.5.1 or greater (but less than Python 3.0) and NumPy 1.3.0 or greater.

#. Install BLAST_. Versions 2.2.16 through 2.2.21 have been tested extensively with PyNAST_, but other versions should work. (Note: we've had trouble with BLAST_ installed via package managers, so it may be best to download directly from NCBI and install per their instructions.)

#. From your command terminal on an OS X or Linux system, change to the directory where you wish to install PyNAST_. You can either download `PyNAST 1.1 from here <https://sourceforge.net/projects/pynast/files/PyNAST%20releases/PyNAST-1.0.tar.gz/download>`_, or if you want the latest development version you can checkout the latest version of PyNAST_ from the SVN repository with the command:
	::
      
		svn co https://pynast.svn.sourceforge.net/svnroot/pynast PyNAST
		
If you downloaded from svn, you will have a new folder in the current working directory called ``PyNAST``. If you downloaded PyNAST-1.1, after untar/unzipping ``PyNAST-1.1.tar.gz`` will have a new directory named ``PyNAST-1.1``. **For consistency, all instructions below will refer to this directory as** ``PyNAST``. You may choose to rename ``PyNAST-1.1`` as ``PyNAST``.

#. Run setup.py. You may need to do this as root (see :ref:`customizing_your_installation` below if this is not an option, or if you'd like to install the PyNAST library code and/or scripts in non-default locations):
	::

		cd PyNAST
		python setup.py install

#. Change to the PyNAST/tests directory:
	::

		cd PyNAST/tests

#. Run the following commands. All tests should pass, unless you don't have MUSCLE_, MAFFT_, and/or ClustalW_ installed. These are optional external software packages, and you will get one test failure per missing software package. You can ignore test failures which indicate that these programs cannot be found.
	::

		python test_logger.py
		python test_util.py

#. If all tests pass, you can get the usage information for the command line version of PyNAST_ with the following command anywhere on your system:
	::
	
		pynast -h

.. _customizing_your_installation:
		
Customizing your installation
=============================

PyNAST consists of library code and a script. By default the script will be installed in ``/usr/local/bin``. This can be customized with the ``-install_scripts`` option like:

::
	
	python setup.py install --install-scripts=/Users/caporaso/bin/
	
You can similarly install the library code in an alternate location using the ``--install-purelib`` option:

::
	
	python setup.py install --install-purelib=/Users/caporaso/temp/


A combination of these options is also possible:	
::
	
	python setup.py install --install-scripts=/Users/caporaso/bin/ --install-purelib=/Users/caporaso/temp/

For a complete discussion of customizations related to the setup.py script, `see this page <http://docs.python.org/install/index.html#alternate-installation-the-home-scheme>`_.

Using the PyNAST command line application
=========================================

After installing the PyNAST_ software as described above, you should download the sample candidate sequences and template alignment. You can then apply the PyNAST_ command line tool as follows:
::
	
	pynast -i candidate_seqs_sample.fasta -t template_sample.fasta

This will result in three files being written to the current working directory: candidate_seqs_sample_pynast_aligned.fasta, candidate_seqs_sample_pynast_log.txt, and candidate_seqs_sample_pynast_fail.fasta, which correspond to the alignment, the run log, and the list of sequences which failed to align, respectively.

To get usage information for the PyNAST_ command line application run:
::
	
	pynast -h
	
	
.. _PyCogent: http://pycogent.sourceforge.net
.. _Python: http://www.python.org
.. _NumPy: http://numpy.scipy.org/
.. _MUSCLE: http://www.drive5.com/muscle/
.. _PyNAST: http://pynast.sourceforge.net
.. _ClustalW: http://www.ebi.ac.uk/Tools/clustalw2/index.html
.. _BLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
.. _MAFFT: http://align.bmr.kyushu-u.ac.jp/mafft/online/server/