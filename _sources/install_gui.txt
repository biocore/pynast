.. Install GUI

*************************************
Installing and using the Mac OS X GUI
*************************************
Download the draft version of the `PyNAST OS X GUI here <https://github.com/downloads/qiime/pynast/PyNAST.app.zip>`_. Unzip the downloaded file to extract the PyNAST_ application. Depending on your system settings, the PyNAST_ application will either be called PyNAST or PyNAST.app. Ensure that your system meets the requirements listed below. If all requirements are met, double-click on the PyNAST_ application to launch PyNAST_. Note that YOU DO NOT NEED PyCogent_ or the PyNAST_ API/command line interface installed to use the PyNAST_ GUI.

Requirements for the PyNAST GUI
===============================
    * An Intel Mac running OS X 10.5 (Leopard).
    * Python_ 2.5 or greater (but less than Python 3.0) and NumPy_ 1.3.0 or greater.
    * ``blastall``, ``formatdb``, and ``bl2seq`` installed in ``/usr/bin/``, ``/usr/local/bin/``, or ``$HOME/bin``. These are all part of NCBI's 'legacy' BLAST_ package, *NOT the BLAST+ package*. Versions 2.2.16 through 2.2.21 have been tested extensively with the PyNAST_ GUI, but other versions should work. (Due to current limitations of the PyNAST_ GUI you need to have the required external software installed in one of these specific locations on your system.)

Optional for the PyNAST GUI
===========================
    * MUSCLE_ installed in ``/usr/bin/``, ``/usr/local/bin/``, or ``$HOME/bin`` if you want to use that for pairwise aligning.
    * ClustalW_ installed in ``/usr/bin/``, ``/usr/local/bin/``, or ``$HOME/bin`` if you want to use that for pairwise aligning.

Limitations in the draft release of the PyNAST GUI
==================================================
    * Not all pairwise aligners are available. Missing options are pair_hmm and MAFFT.
    * Users must place external executables in specific locations for PyNAST_ to find them, rather than PyNAST_ looking in user-defined locations. This will be addressed by adding a preferences box where users can define where these executables are stored.
    * No help text within the application.

.. _PyCogent: http://pycogent.sourceforge.net
.. _Python: http://www.python.org
.. _NumPy: http://numpy.scipy.org/
.. _MUSCLE: http://www.drive5.com/muscle/
.. _PyNAST: http://qiime.org/pynast
.. _ClustalW: http://www.ebi.ac.uk/Tools/clustalw2/index.html
.. _BLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.22/
