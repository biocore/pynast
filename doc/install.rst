.. install_:

********************************************************
Installing and using the PyNAST command line application
********************************************************

PyNAST is installable with ``pip``, but also has non-python dependencies so installation is a little (but not much) more complicated than ``pip install pynast``.

Required software
=================
PyNAST_ is built on the PyCogent_ package, and uses uclust_. 

You must have uclust `v1.2.22q <http://www.drive5.com/uclust/downloads1_2_22q.html>`_ installed to run PyNAST_. You should first obtain uclust, and install it according to the instructions provided by their authors.

PyNAST_ also depends on PyCogent_ and NumPy_. If you install PyNAST using the ``pip`` instructions below, these will be installed for you. Alternatively, you can download and install each according to the instructions on the project websites.

Optional software
=================
If you'd like to perform pairwise alignments using BLAST_, MUSCLE_, MAFFT_, or ClustalW_, you must have those programs installed on your machine and in your system path. Currently tested versions are BLAST_ 2.2.22, MUSCLE_ v3.8.31, MAFFT v6.602b (**MAFFT v6.925b is known to NOT work with PyNAST**), and ClustalW 1.81 or 1.83. Note that PyNAST makes use of the legacy BLAST software, not BLAST+.

pip installation of the latest stable PyNAST release: the easy way
==================================================================

#. Download and install uclust_. Binaries are available (`uclust v1.2.22q binaries <http://www.drive5.com/uclust/downloads1_2_22q.html>`_).

#. Run ``pip install numpy``

#. Run ``pip install pynast``

That's it! You should now have a working ``PyNAST`` installation. (Note that you must run ``pip` in two steps, due to an `issue with PyCogent <https://github.com/pycogent/pycogent/issues/59>`_.)

pip installation of the latest development version of PyNAST: the easy way
==========================================================================

#. Download and install uclust_. Binaries are available (`uclust v1.2.22q binaries <http://www.drive5.com/uclust/downloads1_2_22q.html>`_).

#. Run ``pip install numpy``

#. Run ``pip install git://github.com/qiime/pynast.git``

That's it! You should now have a working ``PyNAST`` installation. (Note that you must run ``pip` in two steps, due to an `issue with PyCogent <https://github.com/pycogent/pycogent/issues/59>`_.)

Manual installation: the harder way
===================================
#. Download PyCogent_ 1.5.3 (`src <http://sourceforge.net/projects/pycogent/files/PyCogent/1.5.3/PyCogent-1.5.3.tgz/download>`_) and its dependencies, Python_ 2.6 or greater (but less than Python 3.0) and NumPy 1.3.0 or greater. PyNAST was tested with Python 2.7.1 and 2.7.2 and NumPy 1.5.1, though other versions may work as well.

#. Download and install uclust_. Binaries are available (`uclust v1.2.22q binaries <http://www.drive5.com/uclust/downloads1_2_22q.html>`_).

#. From your command terminal on an OS X or Linux system, change to the directory where you wish to install PyNAST_. You can either download `PyNAST from PyPI <https://pypi.python.org/pypi/pynast>`_, or if you want the latest development version you can checkout the latest version of PyNAST_ from the GitHub repository with the command: ::

    git clone git://github.com/qiime/pynast.git pynast

#. If you downloaded from GitHub, you will have a new folder in the current working directory called ``pynast``. If you downloaded PyPI, after untar/unzipping the ``tar.gz`` file, you will have a new directory named ``pynast-<version>``, where ``<version>`` is the PyNAST version number. Change to whichever of these directories is relevant for your install procedure, for example::

    cd pynast

#. Run ``setup.py``. You may need to do this as root (see :ref:`customizing_your_installation` below if this is not an option, or if you'd like to install the PyNAST library code and/or scripts in non-default locations)::

    python setup.py install

#. To test your installation, you should run the test suite with the following command. All tests should pass, unless you don't have MUSCLE_, MAFFT_, and/or ClustalW_ installed. These are optional external software packages, and you will get one test failure per missing software package. You can ignore test failures which indicate that these programs cannot be found. ::

    python tests/all_tests.py

#. If all tests pass, you are ready to use PyNAST. You can get the usage information for the command line version of PyNAST_ with the following command anywhere on your system: ::

    cd
    pynast -h

.. _customizing_your_installation:

Customizing your installation
=============================
PyNAST consists of library code and a script. By default the script will be installed in ``/usr/local/bin``. This can be customized with the ``--install_scripts`` option: ::

    python setup.py install --install-scripts=$HOME/bin/

You can similarly install the library code in an alternate location using the ``--install-purelib`` option: ::

    python setup.py install --install-purelib=$HOME/lib/

A combination of these options is also possible: ::

    python setup.py install --install-scripts=$HOME/bin/ --install-purelib=$HOME/lib/

For a complete discussion of customizations related to the setup.py script, `see this page <http://docs.python.org/install/index.html#alternate-installation-the-home-scheme>`_.

If you specify an alternate directory for ``--install-purelib``, you'll need to ensure that python knows where to look for the pynast module. Following the example above, you would do this with the following commands: ::

    echo "export PYTHONPATH=$HOME/lib/:$PYTHONPATH" >> $HOME/.bashrc
    source $HOME/.bashrc

Similarly, if you specify an alternate directory for ``--install-scripts``, you'll need to ensure that the shell knows where to look for executable files. Following the example above, you would do this with the following commands: ::

    echo "export PATH=$HOME/bin/:$PATH" >> $HOME/.bashrc
    source $HOME/.bashrc

Using the PyNAST command line application
=========================================
After installing the PyNAST_ software as described above, you should download the sample candidate sequences and template alignment. You can then apply the PyNAST_ command line tool as follows: ::

    pynast -i candidate_seqs_sample.fasta -t template_sample.fasta

This will result in three files being written to the current working directory: :file:`candidate_seqs_sample_pynast_aligned.fasta`, :file:`candidate_seqs_sample_pynast_log.txt`, and :file:`candidate_seqs_sample_pynast_fail.fasta`, which correspond to the alignment, the run log, and the list of sequences which failed to align, respectively.

To get usage information for the PyNAST_ command line application run: ::

    pynast -h

.. _PyCogent: http://pycogent.sourceforge.net
.. _Python: http://www.python.org
.. _NumPy: http://numpy.scipy.org/
.. _MUSCLE: http://www.drive5.com/muscle/
.. _PyNAST: http://qiime.org/pynast
.. _ClustalW: http://www.ebi.ac.uk/Tools/clustalw2/index.html
.. _BLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.22/
.. _MAFFT: http://align.bmr.kyushu-u.ac.jp/mafft/online/server/
.. _uclust: http://www.drive5.com/uclust/
