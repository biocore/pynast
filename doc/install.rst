.. install_:

*************************************************************
Installing and using the PyNAST command line application
*************************************************************

Downloading PyNAST
==================
You can download the latest development branch of PyNAST here with the command:

	::
	 
		svn co https://pynast.svn.sourceforge.net/svnroot/pynast PyNAST

Required software
=================
PyNAST_ is built on the PyCogent_ package, and uses uclust_. You must have PyCogent_ 1.4.1 and uclust `v1.1.579 <http://www.drive5.com/uclust/downloads1_1_579.html>`_ installed to run PyNAST_. You should first obtain these software packages, and install them according to the instructions provided by their authors.

Optional software
=================
If you'd like to perform pairwise alignments using BLAST_, MUSCLE_, MAFFT_, or ClustalW_, you must have those programs installed on your machine and in your system path. Currently tested versions are BLAST_ 2.2.22, MUSCLE_ v3.6, MAFFT v6.602b, and ClustalW 1.81. Note that PyNAST makes use of the legacy BLAST software, not BLAST+.

Installation steps
==================
#. Download PyCogent_ 1.4.1 and its dependencies, Python_ 2.6 or greater (but less than Python 3.0) and NumPy 1.3.0 or greater.

#. Download and install uclust_. Binaries are available, or you can install from source.

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

		cd tests

#. Run the test suite with the following command. All tests should pass, unless you don't have BLAST_, MUSCLE_, MAFFT_, and/or ClustalW_ installed. These are optional external software packages, and you will get one test failure per missing software package. You can ignore test failures which indicate that these programs cannot be found.
	::

		python all_tests.py

#. If all tests pass, you can get the usage information for the command line version of PyNAST_ with the following command anywhere on your system:
	::
		
		cd
		pynast -h

.. _customizing_your_installation:
		
Customizing your installation
=============================

PyNAST consists of library code and a script. By default the script will be installed in ``/usr/local/bin``. This can be customized with the ``-install_scripts`` option like::
	
	python setup.py install --install-scripts=/home/pynast_user/bin/
	
You can similarly install the library code in an alternate location using the ``--install-purelib`` option::
	
	python setup.py install --install-purelib=/home/pynast_user/lib/


A combination of these options is also possible::
	
	python setup.py install --install-scripts=/home/pynast_user/bin/ --install-purelib=/home/pynast_user/lib/

For a complete discussion of customizations related to the setup.py script, `see this page <http://docs.python.org/install/index.html#alternate-installation-the-home-scheme>`_.

If you specify an alternate directory for ``--install-purelib``, you'll need to ensure that python knows where to look for the pynast module. Following the example above, you would do this with the following commands::

	echo "export PYTHONPATH=/home/pynast_user/lib/:$PYTHONPATH" >> /home/pynast_user/.bashrc
	source /home/pynast_user/.bashrc
	
Similarly, if you specify an alternate directory for ``--install-scripts``, you'll need to ensure that the shell knows where to look for executable files. Following the example above, you would do this with the following commands::

	echo "export PATH=/home/pynast_user/bin/:$PATH" >> /home/pynast_user/.bashrc
	source /home/pynast_user/.bashrc
	
	

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
.. _uclust: http://www.drive5.com/uclust/