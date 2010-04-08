#!/usr/bin/env python

from __future__ import division
import sys
from cogent import LoadSeqs, DNA
from cogent.util.misc import remove_files
from cogent.core.alignment import DenseAlignment
from cogent.app.util import get_tmp_filename
from cogent.app.muscle import align_unaligned_seqs as muscle_align_unaligned_seqs
from cogent.app.mafft import align_unaligned_seqs as mafft_align_unaligned_seqs
from cogent.app.clustalw import align_unaligned_seqs as clustal_align_unaligned_seqs
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.unit_test import TestCase, main
from pynast.util import (align_two_seqs, reintroduce_template_spacing, 
 adjust_alignment, nearest_gap, pynast_seq,
 introduce_terminal_gaps, UnalignableSequenceError, pynast_seqs,
 pair_hmm_align_unaligned_seqs, blast_align_unaligned_seqs, ipynast_seqs,
 remove_template_terminal_gaps)
from pynast.logger import NastLogger

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The PyNAST Project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

class PyNastTests(TestCase):
    """ Tests of the PyNAST functionality
    """
    
    def setUp(self):
        """ """
        
        self.full_length_test1_input_seqs =\
         LoadSeqs(data=input_seqs1_fasta,moltype=DNA,aligned=False)
        self.full_length_test1_input_seqs_lines = input_seqs1_fasta.split('\n')
        self.full_length_test1_template_aln = \
         LoadSeqs(data=pynast_test_template_fasta1,moltype=DNA,aligned=DenseAlignment)
        self.full_length_test1_expected_aln = \
         LoadSeqs(data=input_seqs1_aligned_fasta,moltype=DNA,aligned=DenseAlignment)
        self.full_length_test1_expected_fail = \
         LoadSeqs(data=input_seqs1_fail_fasta,moltype=DNA,aligned=False)
         
        self.full_length_test2_input_seqs =\
         LoadSeqs(data=input_seqs2_fasta,moltype=DNA,aligned=False)
        self.full_length_test2_template_aln = \
         LoadSeqs(data=pynast_test_template_fasta2,moltype=DNA,aligned=DenseAlignment)
        
        self.input_seqs_gaps = input_seqs_gaps.split('\n')
        
        self.files_to_remove = []
        
        self.log_filename = \
            get_tmp_filename(prefix='PyNastTest', suffix='.log')
        self.files_to_remove.append(self.log_filename)
        # touch the log file, so we don't get an error trying to remove it
        # if a test doesn't create it
        open(self.log_filename,'w').close()

    def tearDown(self):
        """ Clean up temporary files created by the tests
        """
        remove_files(self.files_to_remove)

    def test_pynast_logging(self):
        """pynast_seqs() should write log file with correct contents
        """
        logger = NastLogger(self.log_filename)
        seqs = [('1','ACGTACGTTAATACCCTGGTAGT'),
                ('2','AA')]
        # testing for side effect - do not collect return value
        pynast_seqs(seqs, db_aln2, min_len=5, logger=logger)

        log_file = open(self.log_filename, 'r')
        header = log_file.readline()
        contents = log_file.read()
        log_file.close()

        self.assertEqual(contents, expected_logfile_contents)

    def test_pynast_logging_for_stringent_user_requirements(self):
        """pynast_seqs() should record info if best hit does not meet min requirements
        """
        logger = NastLogger(self.log_filename)
        seqs = [('1','ACGTACGTTAATACCCTGGTAGT')]
        # testing for side effect - do not collect return value
        pynast_seqs(seqs, db_aln2, min_len=500, logger=logger)

        log_file = open(self.log_filename, 'r')
        header = log_file.readline()
        contents = log_file.read()
        log_file.close()

        self.assertEqual(contents, expected_stringent_logfile_contents)
        

    def test_pynast_seqs_fail(self):
        """ pynast_seqs: returns expected fail list for sample data
        """
        actual = pynast_seqs(\
         MinimalFastaParser(self.full_length_test1_input_seqs_lines),\
         self.full_length_test1_template_aln,\
         min_len=1000,min_pct=75.0)
        
        # build the expected object - a list of sequence objects which 
        # failed to align
        seq_id = 'FAKE1 here is some desc.73602 tag1;tag2, tag3:tag4'
        expected = [\
         DNA.makeSequence(self.full_length_test1_expected_fail.getSeq(seq_id),\
         Name=seq_id)]
        
        self.assertEqual(actual[1],expected)
      
    def test_pynast_seqs_exact_matches(self):
        """ pynast_seqs: perfectly aligns several exact template matches
        """
        template_aln = self.full_length_test1_template_aln
        
        # Build the expected result object, which is a list of 
        # dna sequence objects where names include the aligned span
        expected_seqs = []
        for n in template_aln.Names:
            expected_seqs.append(\
             DNA.makeSequence(\
              str(template_aln.getGappedSeq(n)),\
              Name='%s 1..%d' % (n,len(template_aln.getSeq(n).degap()))))
          
        expected_aln = LoadSeqs(data=expected_seqs,\
            moltype=DNA,aligned=DenseAlignment) 
        input_seqs = self.full_length_test1_template_aln.degap()
        
        # run pynast_seqs on the input sequences
        actual = pynast_seqs(input_seqs.todict().items(),\
         template_aln,\
         min_len=1000,min_pct=75.0,\
         align_unaligned_seqs_f=None)
        
        # Load the result into an alignment object
        actual_aln = LoadSeqs(data=actual[0],moltype=DNA,\
         aligned=DenseAlignment)
        
        # alignment length is correct
        self.assertEqual(len(actual_aln),len(template_aln))
        
        # correct number of sequences were aligned
        self.assertEqual(actual_aln.getNumSeqs(),expected_aln.getNumSeqs())
        
        # same collection of seq ids is returned
        actual_names = actual_aln.Names
        actual_names.sort()
        expected_names = expected_aln.Names
        expected_names.sort()
        self.assertEqual(actual_names,expected_names)
        
        # all sequence lengths match expected sequence lengths (ie, no
        # missing bases)
        for seq_id in actual_aln.Names:
            self.assertEqual(\
             len(actual_aln.getSeq(seq_id)),\
             len(expected_aln.getSeq(seq_id)))      

        # resulting list of dna sequence objects is as expected
        # (this would take care of some of the above tests, but testing
        # aspects individually makes it easier to diagnose failures)
        actual[0].sort()
        expected_seqs.sort()
        self.assertEqual(actual[0],expected_seqs)
        
        # fail list is empty
        self.assertEqual(actual[1],[])

    def test_pynast_seqs_aligned_full_length(self):
        """ pynast_seqs: pynast results at least 95% identical to NAST results
        
            A note on this test: In the initial versions of PyNAST, I
            wanted the alignments to be exactly like those resulting from
            NAST (e.g., in PyNAST 1.0). I've since abandoned that, in favor
            of getting improved alignments. This test was modified after 
            PyNAST 1.0, and I'm now only testing that the alignments
            are similar to those derived from NAST. This test may be
            of little use, but it is a nice test of the code on 
            full-length sequences, so I hesitate to delete it. 
            -Greg (24 Mar 2010)
        
        """
        template_aln = self.full_length_test1_template_aln
        expected_aln = self.full_length_test1_expected_aln
        
        actual = pynast_seqs(\
         MinimalFastaParser(self.full_length_test1_input_seqs_lines),\
         template_aln,\
         align_unaligned_seqs_f=None) 

        # Build the expected result object, which is a list of 
        # dna sequence objects where names include the aligned span
        expected_seqs = []
        for n in expected_aln.Names:
            expected_seqs.append(\
             DNA.makeSequence(str(expected_aln.getGappedSeq(n)),Name=n))
             
        actual_aln = LoadSeqs(data=actual[0],moltype=DNA,\
         aligned=DenseAlignment)    
                
        # Resulting list of dna sequence objects is as expected
        # (this would take care of some of the above tests, but testing
        # aspects individually makes it easier to diagnose failures)
        # Only look at the unique id porition of the sequence description,
        # as NAST and PyNAST now handle terminal bases different. NAST
        # does local alignments, so sometimes loses terminal bases. PyNAST
        # does global alignments, so the candidate only lose terminal bases
        # if they introduce terminal gaps in the template alignments.
        a_list = [(a.Name.split()[0], a) for a in actual[0]]
        e_list = [(e.Name.split()[0], e) for e in expected_seqs]
        a_list.sort()
        e_list.sort()
        
        for a,e in zip(a_list,e_list):
            # first component of names are equal
            self.assertEqual(a[0],e[0])
            a_seq = a[1]
            e_seq = e[1]
            count_same = 0
            for i in range(len(a_seq)):
                if a_seq[i] == e_seq[i]: count_same += 1
            percent_same = count_same/len(a_seq)
            self.assertTrue(percent_same >= 0.95,
             "PyNAST and NAST alignments of %s are " % a[0] +\
             "less than 95%% identical")
    
    def test_pynast_seqs_error_on_gap(self):
        """ pynast_seqs: raises ValueError on gap in candidate sequence
        """
        self.assertRaises(ValueError,pynast_seqs,
         MinimalFastaParser(self.input_seqs_gaps),\
         self.full_length_test1_template_aln,\
         min_len=1000,min_pct=75.0)
        
    def test_pynast_seqs_simple(self):
        """pynast_seqs: fns with simple test data
        """
        candidate_seqs = [\
         ('1','ACGTACGTTAATACCCTGGTAGT'),\
         ('2','ACGTACGTTAATACCCTGGTAGT'),\
         ('3','AA')]
         
        expected_aln = [\
         DNA.makeSequence('ACGTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---',Name='1'),\
         DNA.makeSequence('ACGTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---',Name='2')]
        expected_fail = [DNA.makeSequence('AA',Name='3')]
        
        actual = pynast_seqs(candidate_seqs,db_aln2,min_len=5,min_pct=75.0)
        self.assertEqual(actual,(expected_aln,expected_fail))
        
         
        # all fail when min_len restricts matches
        expected_aln = []
        expected_fail = [\
         DNA.makeSequence('ACGTACGTTAATACCCTGGTAGT',Name='1'),\
         DNA.makeSequence('ACGTACGTTAATACCCTGGTAGT',Name='2'),\
         DNA.makeSequence('AA',Name='3')]
        
        actual = pynast_seqs(candidate_seqs,db_aln2,min_len=5000,min_pct=75.0)
        
        self.assertEqual(actual,(expected_aln,expected_fail))
    
    def test_pynast_seqs_simple_alt_pairwise(self):
        """pynast_seqs: fns with alt pairwise aligner
        """
        # tests that the order of the returned sequences is correct
        # as this is easy to screw up
        candidate_seqs = [('1','AGCCCCTTTT')]
        template_aln = LoadSeqs(data=dict([
            ('2','ACCC-----CCTTTT')]),\
            moltype=DNA,aligned=DenseAlignment)
        expected_aln = [DNA.makeSequence('AGCC-----CCTTTT',Name='1')]
        expected_fail = []
        
        actual = pynast_seqs(candidate_seqs,template_aln,
                     min_len=5,min_pct=75.0,\
                     align_unaligned_seqs_f=pair_hmm_align_unaligned_seqs)
        self.assertEqual(actual,(expected_aln,expected_fail))
        
        
        # tests that the aligner was actually applied, as it's
        # nearly impossible to get different alignments with
        # different aligners on these short test sequences --
        # therefore test with a fake aligner that alters the sequence
        def fake_aligner(seqs,moltype,params={}):
            return LoadSeqs(data=[('candidate','AGGGGGTTTT'),
                                   ('template', 'ACCCCCTTTT')],moltype=DNA)

        candidate_seqs = [('1','ACCCCCTTTT')]
        template_aln = LoadSeqs(data=dict([
             ('2','ACCC-----CCTTTT')]),\
             moltype=DNA,aligned=DenseAlignment)
        expected_aln = [DNA.makeSequence('AGGG-----GGTTTT',Name='1')]
        expected_fail = []
        actual = pynast_seqs(candidate_seqs,template_aln,
                              min_len=5,min_pct=75.0,\
                              align_unaligned_seqs_f=fake_aligner)
        self.assertEqual(actual,(expected_aln,expected_fail))

       
    def test_ipynast_seqs_simple(self):
        """ipynast_seqs: fns with simple test data
        """
        candidate_seqs = [\
         ('1','ACGAACGTTAATACCCTGGAAGT'),\
         ('2','ACGTACGTTAATACCCTGGTAGT'),\
         ('3','AA')]
         
        expected = [\
         (DNA.makeSequence(\
          'ACGAACGT-TA--ATA-C-----CC-T-G-GAA-G-T---',Name='1'),0),\
         (DNA.makeSequence(\
          'ACGTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---',Name='2'),0),\
         (DNA.makeSequence('AA',Name='3'),1)]
        
        actual = list(ipynast_seqs(\
         candidate_seqs,db_aln2,min_len=5,min_pct=75.0))
        
        self.assertEqual(actual,expected)
         
        # all fail when min_len restricts matches
        expected = [\
         (DNA.makeSequence('ACGAACGTTAATACCCTGGAAGT',Name='1'),2),\
         (DNA.makeSequence('ACGTACGTTAATACCCTGGTAGT',Name='2'),2),\
         (DNA.makeSequence('AA',Name='3'),1)]
        
        actual = list(ipynast_seqs(\
         candidate_seqs,db_aln2,min_len=5000,min_pct=75.0))
        
        self.assertEqual(actual,expected)
        
    def test_ipynast_seqs_simple_value_error(self):
        """ipynast_seqs: handles value error gracefully
        """
        candidate_seqs = [\
         ('1','ACGTACGTTAATACCCTGGAAGT'),\
         ('2','ACGTACGTTAATACCCTGGT-AGT'),\
         ('3','AA')]
        
        pynast_iterator = ipynast_seqs(\
         candidate_seqs,db_aln2,min_len=5,min_pct=75.0)
        
        self.assertRaises(ValueError,list,pynast_iterator)
        
    def test_ipynast_seqs_real_data(self):
        """ipynast_seqs_real_data: sanity check with real data
        """
        actual = list(ipynast_seqs(\
         self.full_length_test2_input_seqs.items(),\
         self.full_length_test2_template_aln,\
         min_len=5,min_pct=75.0))
        # correct number of results returned
        self.assertEqual(len(actual),1)
         
        actual = list(ipynast_seqs(\
         self.full_length_test1_input_seqs.items(),\
         self.full_length_test1_template_aln,\
         min_len=5,min_pct=75.0))
        # correct number of results returned
        self.assertEqual(len(actual),6)

        
    def test_pynast_seqs_simple_status_callback(self):
        """pynast_seqs: status callback functions as expected
        """
        candidate_seqs = [\
         ('1','ACGTACGTTAATACCCTGGTAGT'),\
         ('2','ACGTACGTTAATACCCTGGTAGT'),\
         ('3','AA')]
         
        expected_aln = [\
         DNA.makeSequence('ACGTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---',Name='1'),\
         DNA.makeSequence('ACGTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---',Name='2')]
        expected_fail = [DNA.makeSequence('AA',Name='3')]
        
        class StatusTracker(object):
            completed_seqs_count = 0
            def update_completed_seqs_count(self,x):
                self.completed_seqs_count += 1
        
        st = StatusTracker()
        self.assertEqual(st.completed_seqs_count,0)
        results = pynast_seqs(candidate_seqs,db_aln2,min_len=5,min_pct=75.0,\
         status_callback_f=st.update_completed_seqs_count)
        
        self.assertEqual(st.completed_seqs_count,3)

    def test_pynast_seq_simple(self):
        """pynast_seq: fns as exp with simple example
        """
        candidate_sequence =\
         DNA.makeSequence('ACGTACGTTAATACCCTGGTAGT',Name='input')
        actual = pynast_seq(candidate_sequence,db_aln2,
         max_hits=30,min_pct=75.0,
         min_len=5,align_unaligned_seqs_f=None)
        
        # check individual components of result object
        expected_template_hit = '5'
        expected_aligned_seq = 'ACGTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---'
        expected_aligned_seq_id = 'input 1..23'
        
        self.assertEqual(actual[0],expected_template_hit)
        self.assertEqual(str(actual[1]),expected_aligned_seq)
        self.assertEqual(actual[1].Name,expected_aligned_seq_id)
        
        # check full result object
        expected = ('5',\
         DNA.makeSequence('ACGTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---',\
         Name='input 1..23')) 
        self.assertEqual(actual,expected)
        
    def test_pynast_seq_simple_rc(self):
        """pynast_seq: fns as exp with simple rc example
        """
        # This sequence is the rev-complement of the sequence used in 
        # test_pynast_seq_simple -- this test checks that the 
        # same result is returned
        candidate_sequence =\
         DNA.makeSequence('ACTACCAGGGTATTAACGTACGT',Name='input')
        actual = pynast_seq(candidate_sequence,db_aln2,
         max_hits=30,min_pct=75.0,
         min_len=5,align_unaligned_seqs_f=None)
        
        # check individual components of result object
        expected_template_hit = '5'
        expected_aligned_seq = 'ACGTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---'
        expected_aligned_seq_id = 'input RC:1..23'
        
        self.assertEqual(actual[0],expected_template_hit)
        self.assertEqual(str(actual[1]),expected_aligned_seq)
        self.assertEqual(actual[1].Name,expected_aligned_seq_id)
        
        # check full result object
        expected = ('5',\
         DNA.makeSequence('ACGTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---',\
         Name='input RC:1..23')) 
        self.assertEqual(actual,expected)        
    
    def test_pynast_seq_10116(self):
        """pynast_seq: real seq that introduces 5' gaps in pw aligned template
        
            The pairwise alignment of this sequence to the template alignment
             results in five prime gaps in the pairwise aligned template. This
             caused a bug in early versions of PyNAST because too many terminal
             gaps were being reintroduced. Therefore keeping this as a real
             test case, essentially of the introduce_terminal_gaps 
             functionality.
        
        """
        candidate_sequence =\
         LoadSeqs(data=input_seq_10116.split('\n'),moltype=DNA).\
         getSeq('10116')
        template_aln = self.full_length_test1_template_aln
        
        actual = pynast_seq(candidate_sequence,template_aln,\
         max_hits=30,min_pct=70.0,min_len=150,\
         align_unaligned_seqs_f=None)
         
        self.assertEqual(len(actual[1]),len(template_aln))
    
        
    def test_pynast_seq_14990(self):
        """pynast_seq: aligning handles input seq longer than best template seq
        """
        template_aln =\
         LoadSeqs(data=template_14990_trimmed.split('\n'),\
         moltype=DNA,aligned=DenseAlignment) 
        candidate_sequence =\
         LoadSeqs(data=input_seq_14990.split('\n'),moltype=DNA).\
         getSeq('14990')
        expected = ('14990_5_and_3_prime_lost_four_bases_each',\
         template_aln.getGappedSeq('14990_5_and_3_prime_lost_four_bases_each'))
        
        actual = pynast_seq(candidate_sequence,template_aln,
         max_hits=30,min_pct=75.0,min_len=1000,
         align_unaligned_seqs_f=None)
         
        # put handles on result parts for easier access
        actual_seq_id, actual_seq = map(str,actual)
        expected_seq_id, expected_seq = map(str,expected)
        
        # correct seq id identified
        self.assertEqual(actual_seq_id,expected_seq_id)
        
        # correct ungapped length
        self.assertEqual(len(actual_seq.replace('-','')),\
                         len(expected_seq.replace('-','')))

        # correct gapped length
        self.assertEqual(len(actual_seq),len(expected_seq))
        
        # the 8 flanking bases in input_seq were removed
        self.assertEqual(len(actual_seq.replace('-','')),\
                         len(candidate_sequence)-8)
        
        # aligned seqs are equal
        self.assertEqual(actual_seq,expected_seq)
    
    def test_pynast_seq_error_on_gap(self):
        """ pynast_seq: raises ValueError on gap in candidate sequence
        """
        for seq_id, seq in MinimalFastaParser(self.input_seqs_gaps):
            # error when gap(s) in seq
            cs = DNA.makeSequence(seq,Name=seq_id)
            self.assertRaises(ValueError,pynast_seq,cs,db_aln2,\
             max_hits=1,min_pct=75.0,min_len=5,align_unaligned_seqs_f=None)
             
            seq = seq.replace('-','').replace('.','')
            # no error when no gaps in seq
            cs = DNA.makeSequence(seq,Name=seq_id)
            r = pynast_seq(cs,db_aln2,\
             max_hits=1,min_pct=70.0,min_len=5,align_unaligned_seqs_f=None)
             
        
    def test_align_two_seqs_with_muscle(self):
        """ align_two_seqs: fns for simple alignments with muscle
        """
        # Only a few trivial cases are tested as it is not the place to
        # test how the aligners functions
        f = muscle_align_unaligned_seqs

        # perfect alignment
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(s1,s2))
    
        # gap added to s2
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATCCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACAT-CCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
    
        # gap added to s1
        s1 = DNA.makeSequence('ACGTACGTACATCCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACAT-CCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
    
        # single mismatch
        s1 = DNA.makeSequence('ACGTACGTACATTCCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(s1,s2))
    
        # truncated sequence (3')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACATACCCT------')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
    
        # truncated sequence (5')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('CGTACATACCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('-----CGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))   
    
        # truncated sequence (5' and 3')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('CGTACATACCCTGGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('-----CGTACATACCCTGGT---')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))   
        
    def test_align_two_seqs_with_pair_hmm(self):
        """ align_two_seqs: fns for simple alignments with pair_hmm alignment
        """
        # Only a few trivial cases are tested as it is not the place to
        # test how the aligners functions
        f = pair_hmm_align_unaligned_seqs

        # perfect alignment
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(s1,s2))
    
        # gap added to s2
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATCCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACAT-CCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
    
        # gap added to s1
        s1 = DNA.makeSequence('ACGTACGTACATCCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACAT-CCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
    
        # single mismatch
        s1 = DNA.makeSequence('ACGTACGTACATTCCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(s1,s2))
    
        # truncated sequence (3')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACATACCCT------')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
    
        # truncated sequence (5')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('CGTACATACCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('-----CGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2)) 
    
        # truncated sequence (5' and 3')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('CGTACATACCCTGGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('-----CGTACATACCCTGGT---')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))   
        
        
    def test_align_two_seqs_with_blast(self):
        """ align_two_seqs: fns for simple alignments with blast (bl2seq)
        """
        # Only a few trivial cases are tested as it is not the place to
        # test how the aligners functions
        f = blast_align_unaligned_seqs
    
        # perfect alignment
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(s1,s2))
            
        # gap added to s2
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATCCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACAT-CCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
            
        # gap added to s1
        s1 = DNA.makeSequence('ACGTACGTACATCCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACAT-CCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
            
        # single mismatch
        s1 = DNA.makeSequence('ACGTACGTACATTCCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(s1,s2))
    
        # truncated sequence (3')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACATACCCT------')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
        # reversed order works as well (ie., extended sequence 3')
        self.assertEqual(align_two_seqs(s2,s1,f),(exp2,exp1))
    
        # truncated sequence (5')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('CGTACATACCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('-----CGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))   
        # reversed order works as well (ie., extended sequence 5')
        self.assertEqual(align_two_seqs(s2,s1,f),(exp2,exp1))   
    
        # truncated sequence (5' and 3')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('CGTACATACCCTGGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('-----CGTACATACCCTGGT---')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))   
        # reversed order works as well (ie., extended sequence 5' and 3')
        self.assertEqual(align_two_seqs(s2,s1,f),(exp2,exp1)) 
        
        # staggered ends
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGT')
        s2 = DNA.makeSequence(     'CGTACATACCCTGGTAGTTT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGT-----')
        exp2 = DNA.makeSequence('-----CGTACATACCCTGGTAGTTT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
        # reversed order works as well 
        self.assertEqual(align_two_seqs(s2,s1,f),(exp2,exp1))
        
    def test_align_two_seqs_with_clustal(self):
        """ align_two_seqs: fns for simple alignments with clustal
        """
        # Only a few trivial cases are tested as it is not the place to
        # test how the aligners function
        f = clustal_align_unaligned_seqs

        # perfect alignment
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(s1,s2))
    
        # gap added to s2
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATCCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACAT-CCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
    
        # gap added to s1
        s1 = DNA.makeSequence('ACGTACGTACATCCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACAT-CCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
    
        # single mismatch
        s1 = DNA.makeSequence('ACGTACGTACATTCCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(s1,s2))
    
        # truncated sequence (3')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACATACCCT------')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
    
        # truncated sequence (5')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('CGTACATACCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('-----CGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))   
    
        # truncated sequence (5' and 3')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('CGTACATACCCTGGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('-----CGTACATACCCTGGT---')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
        
    def test_align_two_seqs_with_mafft(self):
        """ align_two_seqs: fns for simple alignments with mafft
        """
        # Only a few trivial cases are tested as it is not the place to
        # test how the aligners functions
        f = mafft_align_unaligned_seqs

        # perfect alignment
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(s1,s2))
    
        # gap added to s2
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATCCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACAT-CCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
    
        # gap added to s1
        s1 = DNA.makeSequence('ACGTACGTACATCCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACAT-CCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
    
        # single mismatch
        s1 = DNA.makeSequence('ACGTACGTACATTCCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(s1,s2))
    
        # truncated sequence (3')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('ACGTACGTACATACCC------T')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))
    
        # truncated sequence (5')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('CGTACATACCCTGGTAGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('-----CGTACATACCCTGGTAGT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))   
    
        # truncated sequence (5' and 3')
        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('CGTACATACCCTGGT')
        exp1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp2 = DNA.makeSequence('-----CGTACATACCCTG---GT')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))   
  
    def test_align_two_seqs_with_fake_aligner(self):
        """ align_two_seqs: fns for simple alignments with fake_aligner
        """
        # Test a fake aligner function which uses the params dict
        def f(seqs,moltype,params={}):
            try:
                res = params['res']
            except KeyError:
                res = 'AAAAAAAAAA'
            seqs = [('template',str(res)), ('candidate',str(res))]
            seqs = LoadSeqs(data=seqs,moltype=moltype,aligned=DenseAlignment)
            return seqs

        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp1 = DNA.makeSequence('AAAAAAAAAA')
        exp2 = DNA.makeSequence('AAAAAAAAAA')
        self.assertEqual(align_two_seqs(s1,s2,f),(exp1,exp2))   

        s1 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        s2 = DNA.makeSequence('ACGTACGTACATACCCTGGTAGT')
        exp1 = DNA.makeSequence('BBB')
        exp2 = DNA.makeSequence('BBB')
        self.assertEqual(align_two_seqs(s1,s2,f,params={'res':'BBB'}),\
         (exp1,exp2))

    def test_reintroduce_template_spacing_template(self):
        """ reintroduce_template_spacing: template example from DeSantis2004
        """
        template = DNA.makeSequence('ATAC-----GTA-AC----GTA---C---G-T-AC-GG')
        pw_aligned_template = DNA.makeSequence('ATACGT-A-ACGTACGTAC--GG')
        pw_aligned_candidate= DNA.makeSequence('C-ACGTTAAACGT-CGTACCCGG')
        template_expected = \
         DNA.makeSequence('ATAC-----GT-A-AC----GTA---C---G-T-AC--GG')
        
        actual = reintroduce_template_spacing(\
         template,pw_aligned_template,pw_aligned_candidate)
        self.assertEqual(actual[0],template_expected)
        
    def test_reintroduce_template_spacing_candidate(self):
        """ reintroduce_template_spacing: candidate example from DeSantis2006
        """
        template = DNA.makeSequence('ATAC-----GTA-AC----GTA---C---G-T-AC-GG')
        pw_aligned_template = DNA.makeSequence('ATACGT-A-ACGTACGTAC--GG')
        pw_aligned_candidate= DNA.makeSequence('C-ACGTTAAACGT-CGTACCCGG')
        candidate_expected = \
         DNA.makeSequence('C-AC-----GTTAAAC----GT----C---G-T-ACCCGG')
        
        actual = reintroduce_template_spacing(\
         template,pw_aligned_template,pw_aligned_candidate)
        self.assertEqual(actual[1],candidate_expected)
        
    def test_reintroduce_template_spacing_new_gaps(self):
        """ reintroduce_template_spacing: new gaps example from DeSantis2006
        """
        template = DNA.makeSequence('ATAC-----GTA-AC----GTA---C---G-T-AC-GG')
        pw_aligned_template = DNA.makeSequence('ATACGT-A-ACGTACGTAC--GG')
        pw_aligned_candidate= DNA.makeSequence('C-ACGTTAAACGT-CGTACCCGG')
        new_gaps_expected = [11,36]
        
        actual = reintroduce_template_spacing(\
         template,pw_aligned_template,pw_aligned_candidate)
        self.assertEqual(actual[2],new_gaps_expected)
        
    def test_reintroduce_template_spacing(self):
        """ reintroduce_template_spacing: example from DeSantis2006
        """
        template = DNA.makeSequence('ATAC-----GTA-AC----GTA---C---G-T-AC-GG')
        pw_aligned_template = DNA.makeSequence('ATACGT-A-ACGTACGTAC--GG')
        pw_aligned_candidate= DNA.makeSequence('C-ACGTTAAACGT-CGTACCCGG')
        template_expected = \
         DNA.makeSequence('ATAC-----GT-A-AC----GTA---C---G-T-AC--GG')
        candidate_expected = \
         DNA.makeSequence('C-AC-----GTTAAAC----GT----C---G-T-ACCCGG')
        new_gaps_expected = [11,36]
        
        actual = reintroduce_template_spacing(\
         template,pw_aligned_template,pw_aligned_candidate)
        self.assertEqual(actual,\
         (template_expected,candidate_expected,new_gaps_expected))
        
    def test_reintroduce_template_spacing_no_change(self):
        """ reintroduce_template_spacing: no changes
        """
        template = DNA.makeSequence('AT-CG')
        actual = reintroduce_template_spacing(\
         template,template,template)
        self.assertEqual(actual,(template,template,[]))        
        
        # different seqs but pw alignment matches template pattern
        template = DNA.makeSequence('ATC-G')
        pw_aligned_template = DNA.makeSequence ('ATC-G')
        pw_aligned_candidate = DNA.makeSequence('ATCCG')
        template_expected = DNA.makeSequence ('ATC-G')
        candidate_expected = DNA.makeSequence('ATCCG')
        
        actual = reintroduce_template_spacing(\
         template,pw_aligned_template,pw_aligned_candidate)
        self.assertEqual(actual,(template_expected,candidate_expected,[]))  
        
    def test_reintroduce_template_spacing_middle(self):
        """ reintroduce_template_spacing: change to non-terminal character
        """
        template = DNA.makeSequence('GTA---C')
        pw_aligned_template = DNA.makeSequence( 'GTAC')
        pw_aligned_candidate = DNA.makeSequence('GT-C')
        template_expected = DNA.makeSequence( 'GTA---C')
        candidate_expected = DNA.makeSequence('GT----C')
        new_gaps_expected = []
        
        actual = reintroduce_template_spacing(\
         template,pw_aligned_template,pw_aligned_candidate)
        self.assertEqual(actual,\
         (template_expected,candidate_expected,new_gaps_expected))
        
        template = DNA.makeSequence('ATAC-----GTA-AC')
        pw_aligned_template = DNA.makeSequence( 'ATACGT-A-AC')
        pw_aligned_candidate = DNA.makeSequence('C-ACGTTAAAC')
        template_expected = DNA.makeSequence( 'ATAC-----GT-A-AC')
        candidate_expected = DNA.makeSequence('C-AC-----GTTAAAC')
        new_gaps_expected = [11]
        
        actual = reintroduce_template_spacing(\
         template,pw_aligned_template,pw_aligned_candidate)
        self.assertEqual(actual,\
         (template_expected,candidate_expected,new_gaps_expected)) 
        
        # single gap in new spot
        template = DNA.makeSequence('GTA-AC')
        pw_aligned_template = DNA.makeSequence( 'GT-A-AC')
        pw_aligned_candidate = DNA.makeSequence('GTTAAAC')
        template_expected = DNA.makeSequence( 'GT-A-AC')
        candidate_expected = DNA.makeSequence('GTTAAAC')
        new_gaps_expected = [2]
        
        actual = reintroduce_template_spacing(\
         template,pw_aligned_template,pw_aligned_candidate)
        self.assertEqual(actual,\
         (template_expected,candidate_expected,new_gaps_expected))   
        
        # existing gap extended
        template = DNA.makeSequence('AC-GG')
        pw_aligned_template = DNA.makeSequence( 'AC--GG')
        pw_aligned_candidate = DNA.makeSequence('ACCCGG')
        template_expected = DNA.makeSequence( 'AC--GG')
        candidate_expected = DNA.makeSequence('ACCCGG')
        new_gaps_expected = [2]
        
        actual = reintroduce_template_spacing(\
         template,pw_aligned_template,pw_aligned_candidate)
        self.assertEqual(actual,\
         (template_expected,candidate_expected,new_gaps_expected)) 

    def test_reintroduce_template_spacing_leading_trailing_gaps_ignored(self):
        """ reintroduce_template_spacing: lead/trailing template gaps ignored
        """       
        # leading gaps
        template = DNA.makeSequence('----AC-GG')
        pw_aligned_template = DNA.makeSequence( 'AC--GG')
        pw_aligned_candidate = DNA.makeSequence('ACCCGG')
        template_expected = DNA.makeSequence( 'AC--GG')
        candidate_expected = DNA.makeSequence('ACCCGG')
        new_gaps_expected = [2]
        
        actual = reintroduce_template_spacing(\
         template,pw_aligned_template,pw_aligned_candidate)
        self.assertEqual(actual,\
         (template_expected,candidate_expected,new_gaps_expected))    
         
        # trailing gaps
        template = DNA.makeSequence('AC-GG---')
        pw_aligned_template = DNA.makeSequence( 'AC--GG')
        pw_aligned_candidate = DNA.makeSequence('ACCCGG')
        template_expected = DNA.makeSequence( 'AC--GG')
        candidate_expected = DNA.makeSequence('ACCCGG')
        new_gaps_expected = [2]
        
        actual = reintroduce_template_spacing(\
         template,pw_aligned_template,pw_aligned_candidate)
        self.assertEqual(actual,\
         (template_expected,candidate_expected,new_gaps_expected))   
           
        # leading/trailing gaps
        template = DNA.makeSequence('-AC-GG---')
        pw_aligned_template = DNA.makeSequence( 'AC--GG')
        pw_aligned_candidate = DNA.makeSequence('ACCCGG')
        template_expected = DNA.makeSequence( 'AC--GG')
        candidate_expected = DNA.makeSequence('ACCCGG')
        new_gaps_expected = [2]
        
        actual = reintroduce_template_spacing(\
         template,pw_aligned_template,pw_aligned_candidate)
        self.assertEqual(actual,\
         (template_expected,candidate_expected,new_gaps_expected))   
        
    def test_adjust_alignment_paper_example(self):
        """ adjust_alignment: example from DeSantis2006
        """
        template = \
          DNA.makeSequence('ATAC-----GT-A-AC----GTA---C---G-T-AC--GG')
        candidate = \
          DNA.makeSequence('C-AC-----GTTAAAC----GT----C---G-T-ACCCGG')
        new_gaps = [11,36]
        # IS THERE A TYPO IN THEIR EXAMPLE? THEY CHANGE GT-A-AC TO 
        # GT-AAC, BUT THAT DOESN'T REALLY MAKE SENSE GIVEN THAT THE
        # TEMPLATE ALIGNMENT IS GTA-AC...
        template_expected = \
         DNA.makeSequence('ATAC-----GTA-AC----GTA---C---G-T-AC-GG')
        candidate_expected = \
         DNA.makeSequence('C-AC----GTTAAAC----GT----C---G-TACCCGG')
      
        actual = adjust_alignment(template,candidate,new_gaps)
        self.assertEqual(actual,(template_expected,candidate_expected))
        
    def test_adjust_alignment(self):
        """ adjust_alignmnet: simple adjustments handled as expected
        """
        # remove a 3' gap
        t = DNA.makeSequence('AA-GGC---ATTAA')
        c = DNA.makeSequence('AATCCTT--AAAAA')
        new_gaps = [2]
        t_expected = DNA.makeSequence('AAGGC---ATTAA')
        c_expected = DNA.makeSequence('AATCCTT-AAAAA')
        self.assertEqual(adjust_alignment(t,c,new_gaps),\
         (t_expected,c_expected))
         
        # remove a 5' gap
        t = DNA.makeSequence('AA-GGC----TTAA')
        c = DNA.makeSequence('AATCCTT--AAAAA')
        new_gaps = [9]
        t_expected = DNA.makeSequence('AA-GGC---TTAA')
        c_expected = DNA.makeSequence('AATCCTT-AAAAA')
        self.assertEqual(adjust_alignment(t,c,new_gaps),\
         (t_expected,c_expected))
         
         # multiple gaps to remove
        t = DNA.makeSequence('AA-GGC----TTAA')
        c = DNA.makeSequence('AATCCTT--AAAAA')
        new_gaps = [2,9]
        t_expected = DNA.makeSequence('AAGGC---TTAA')
        c_expected = DNA.makeSequence('AATCCTTAAAAA')
        self.assertEqual(adjust_alignment(t,c,new_gaps),\
         (t_expected,c_expected))
        
         
    def test_adjust_alignment_multiple_adjancent_new_gaps(self):
        """ adjust_alignmnet: multiple adjacent new gaps handled as expected
        """
        t = DNA.makeSequence('AA--GC---ATTAA')
        c = DNA.makeSequence('AATCCTT--AAAAA')
        new_gaps = [2,3]
        t_expected = DNA.makeSequence('AAGC---ATTAA')
        c_expected = DNA.makeSequence('AATCCTTAAAAA')
        actual = adjust_alignment(t,c,new_gaps)
        # print ''
        # print actual[0]
        # print t_expected
        self.assertEqual(actual,(t_expected,c_expected))
        
        t = DNA.makeSequence('AATTGCG---CAT')
        c = DNA.makeSequence('AA---CTTTTAAA')
        new_gaps = [7,8,9]
        t_expected = DNA.makeSequence('AATTGCGCAT')
        c_expected = DNA.makeSequence('AACTTTTAAA')
        actual = adjust_alignment(t,c,new_gaps)
        # print ''
        # print actual[0]
        # print t_expected
        self.assertEqual(actual,(t_expected,c_expected))
        
        t = DNA.makeSequence('AATTGCG---CAT')
        c = DNA.makeSequence('AA-CTTTTTA-A-')
        new_gaps = [7,8,9]
        t_expected = DNA.makeSequence('AATTGCGCAT')
        c_expected = DNA.makeSequence('AACTTTTTAA')
        actual = adjust_alignment(t,c,new_gaps)
        # print ''
        # print actual[0]
        # print t_expected
        self.assertEqual(actual,(t_expected,c_expected))
        
    def test_nearest_gap(self):
        """nearest_gap: functions with single gap in seq
        """
        seq = 'AAA-AAAA'
        for pos in range(len(seq)):
            self.assertEqual(nearest_gap(seq,pos),3)
            
        seq = '-ACGTACGT'
        for pos in range(len(seq)):
            self.assertEqual(nearest_gap(seq,pos),0)
            
        seq = 'ACGTACGT-'
        for pos in range(len(seq)):
            self.assertEqual(nearest_gap(seq,pos),8)
            
    def test_nearest_gap_mutliple_gaps(self):
        """nearest_gap: handles multiple gaps in same sequence
        """
        seq = 'ACG-TT-AACC--TAAT'
        self.assertEqual(nearest_gap(seq,0),3)
        self.assertEqual(nearest_gap(seq,1),3)
        self.assertEqual(nearest_gap(seq,2),3)
        self.assertEqual(nearest_gap(seq,3),3)
        self.assertEqual(nearest_gap(seq,4),3)
        self.assertEqual(nearest_gap(seq,5),6)
        self.assertEqual(nearest_gap(seq,6),6)
        self.assertEqual(nearest_gap(seq,7),6)
        self.assertEqual(nearest_gap(seq,8),6)
        self.assertEqual(nearest_gap(seq,9),11)
        self.assertEqual(nearest_gap(seq,10),11)
        self.assertEqual(nearest_gap(seq,11),11)
        self.assertEqual(nearest_gap(seq,12),12)
        self.assertEqual(nearest_gap(seq,13),12)
        self.assertEqual(nearest_gap(seq,14),12)
        self.assertEqual(nearest_gap(seq,15),12)
        self.assertEqual(nearest_gap(seq,16),12)
        
    def test_nearest_gap_ambiguous(self):
        """nearest_gap: handles ambiguous cases by chosing the 5' position
        
            Not certain that this is how this should be handled... Maybe 
            revisit by seeing which way gives the better alignment?
        """
        seq = '-A-A-A-'
        self.assertEqual(nearest_gap(seq,1),0)
        self.assertEqual(nearest_gap(seq,3),2)
        self.assertEqual(nearest_gap(seq,5),4)
        
        
    def test_nearest_gap_handles_error(self):
        """nearest_gap: errors are handled correctly
        """
        seq = 'AA-AAA'
        self.assertRaises(IndexError,nearest_gap,seq,22)
        self.assertRaises(IndexError,nearest_gap,seq,-1)
        
        seq = 'AAA'
        self.assertRaises(UnalignableSequenceError,nearest_gap,seq,1)
        
    def test_introduce_terminal_gaps_simple(self):
        """introduce_terminal_gaps: functions as expected
        """
        # no terminal gaps
        template = DNA.makeSequence('AAA',Name='t')
        aligned_candidate = DNA.makeSequence('AAA',Name='ac')
        aligned_template = DNA.makeSequence('AAA',Name='at')
        actual = introduce_terminal_gaps(\
            template,aligned_template,aligned_candidate)
        expected = DNA.makeSequence('AAA',Name='ac')
        self.assertEqual(actual,expected)
        
        # 5' terminal gaps only
        template = DNA.makeSequence('-AAA',Name='t')
        aligned_candidate = DNA.makeSequence('AAA',Name='ac')
        aligned_template = DNA.makeSequence('AAA',Name='at')
        actual = introduce_terminal_gaps(\
            template,aligned_template,aligned_candidate)
        expected = DNA.makeSequence('-AAA',Name='ac')
        self.assertEqual(actual,expected)
        
        template = DNA.makeSequence('-----AAA',Name='t')
        aligned_candidate = DNA.makeSequence('AAA',Name='ac')
        aligned_template = DNA.makeSequence('AAA',Name='at')
        actual = introduce_terminal_gaps(\
            template,aligned_template,aligned_candidate)
        expected = DNA.makeSequence('-----AAA',Name='ac')
        self.assertEqual(actual,expected)
        
        # 3' terminal gaps only
        template = DNA.makeSequence('ACG--',Name='t')
        aligned_candidate = DNA.makeSequence('ACG',Name='ac')
        aligned_template = DNA.makeSequence('AAA',Name='at')
        actual = introduce_terminal_gaps(\
            template,aligned_template,aligned_candidate)
        expected = DNA.makeSequence('ACG--',Name='ac')
        self.assertEqual(actual,expected)
        
        template = DNA.makeSequence('ACCTG----',Name='t')
        aligned_candidate = DNA.makeSequence('ACGGG',Name='ac')
        aligned_template = DNA.makeSequence('ACCTG',Name='at')
        actual = introduce_terminal_gaps(\
            template,aligned_template,aligned_candidate)
        expected = DNA.makeSequence('ACGGG----',Name='ac')
        self.assertEqual(actual,expected)
        
        # 5' and 3' terminal gaps
        template = DNA.makeSequence('---AC--CTG----',Name='t')
        aligned_candidate = DNA.makeSequence('ACTTGGG',Name='ac')
        aligned_template = DNA.makeSequence( 'AC--CTG',Name='at')
        actual = introduce_terminal_gaps(\
            template,aligned_template,aligned_candidate)
        expected = DNA.makeSequence('---ACTTGGG----',Name='ac')
        self.assertEqual(actual,expected)
        
    def test_introduce_terminal_gaps_existing_terminal_template_gaps(self):
        """introduce_terminal_gaps: aligned template already has terminal gaps
        """
        
        # one 5' gap in aligned_template
        template = DNA.makeSequence('---AAA',Name='t')
        aligned_candidate = DNA.makeSequence('AAAA',Name='ac')
        aligned_template = DNA.makeSequence('-AAA',Name='at')
        actual = introduce_terminal_gaps(\
            template,aligned_template,aligned_candidate)
        expected = DNA.makeSequence('--AAAA',Name='ac')
        self.assertEqual(actual,expected)
        
        # multiple 5' gaps in aligned_template
        template = DNA.makeSequence('---AAA',Name='t')
        aligned_candidate = DNA.makeSequence('AAAAAA',Name='ac')
        aligned_template = DNA.makeSequence( '---AAA',Name='at')
        actual = introduce_terminal_gaps(\
            template,aligned_template,aligned_candidate)
        expected = DNA.makeSequence('AAAAAA',Name='ac')
        self.assertEqual(actual,expected)
        
        # one 3' gap in aligned_template
        template = DNA.makeSequence('AAA---',Name='t')
        aligned_candidate = DNA.makeSequence('AAAA',Name='ac')
        aligned_template = DNA.makeSequence( 'AAA-',Name='at')
        actual = introduce_terminal_gaps(\
            template,aligned_template,aligned_candidate)
        expected = DNA.makeSequence('AAAA--',Name='ac')
        self.assertEqual(actual,expected)
        
        # multiple 3' gaps in aligned_template
        template = DNA.makeSequence('AAA---',Name='t')
        aligned_candidate = DNA.makeSequence('AAAAAA',Name='ac')
        aligned_template = DNA.makeSequence( 'AAA---',Name='at')
        actual = introduce_terminal_gaps(\
            template,aligned_template,aligned_candidate)
        expected = DNA.makeSequence('AAAAAA',Name='ac')
        self.assertEqual(actual,expected)
        
        # 5 prime, 3 prime gaps in aligned_template
        template = DNA.makeSequence('--CAA---',Name='t')
        aligned_candidate = DNA.makeSequence('GCAAT',Name='ac')
        aligned_template = DNA.makeSequence( '-CAA-',Name='at')
        actual = introduce_terminal_gaps(\
            template,aligned_template,aligned_candidate)
        expected = DNA.makeSequence('-GCAAT--',Name='ac')
        self.assertEqual(actual,expected)
        
        # internal, 5', 3' gaps
        template = DNA.makeSequence('--CATA---',Name='t')
        aligned_candidate = DNA.makeSequence('GCA-AT',Name='ac')
        aligned_template = DNA.makeSequence( '-CATA-',Name='at')
        actual = introduce_terminal_gaps(\
            template,aligned_template,aligned_candidate)
        expected = DNA.makeSequence('-GCA-AT--',Name='ac')
        self.assertEqual(actual,expected)
        
    def test_remove_template_terminal_gaps(self):
        """ removing terminal gaps functions as expected """
        # no template terminal gaps
        candidate = DNA.makeSequence('--CGTTGG-',Name='c')
        template  = DNA.makeSequence('ACCGT-GGA',Name='t')
        actual = remove_template_terminal_gaps(candidate,template)
        expected = (DNA.makeSequence('--CGTTGG-',Name='c 1..6'),template)
        self.assertEqual(actual[0].Name,expected[0].Name)
        self.assertEqual(actual[1].Name,expected[1].Name)
        self.assertEqual(actual,expected)
        
        candidate = DNA.makeSequence('',Name='c')
        template  = DNA.makeSequence('',Name='t')
        actual = remove_template_terminal_gaps(candidate,template)
        expected = (candidate,template)
        self.assertEqual(actual[0].Name,expected[0].Name)
        self.assertEqual(actual[1].Name,expected[1].Name)
        self.assertEqual(actual,expected)
        
        # 5' template terminal gaps
        candidate = DNA.makeSequence('ACCGTTGGA',Name='c')
        template  = DNA.makeSequence('--CGT-GGA',Name='t')
        actual = remove_template_terminal_gaps(candidate,template)
        expected = (DNA.makeSequence('CGTTGGA',Name='c 3..9'),
                    DNA.makeSequence('CGT-GGA',Name='t'))
        self.assertEqual(actual[0].Name,expected[0].Name)
        self.assertEqual(actual[1].Name,expected[1].Name)
        self.assertEqual(actual,expected)
        
        candidate = DNA.makeSequence('ACCGTTGGA',Name='c')
        template  = DNA.makeSequence('-CCGT-GGA',Name='t')
        actual = remove_template_terminal_gaps(candidate,template)
        expected = (DNA.makeSequence('CCGTTGGA',Name='c 2..9'),
                    DNA.makeSequence('CCGT-GGA',Name='t'))
        self.assertEqual(actual[0].Name,expected[0].Name)
        self.assertEqual(actual[1].Name,expected[1].Name)
        self.assertEqual(actual,expected)
        
        # 3' template terminal gaps
        candidate = DNA.makeSequence('ACCGTTGGA',Name='c')
        template  = DNA.makeSequence('ACCGT-GG-',Name='t')
        actual = remove_template_terminal_gaps(candidate,template)
        expected = (DNA.makeSequence('ACCGTTGG',Name='c 1..8'),
                    DNA.makeSequence('ACCGT-GG',Name='t'))
        self.assertEqual(actual[0].Name,expected[0].Name)
        self.assertEqual(actual[1].Name,expected[1].Name)
        self.assertEqual(actual,expected)
        
        candidate = DNA.makeSequence('ACCGTTGGA',Name='c')
        template  = DNA.makeSequence('ACCGT-G--',Name='t')
        actual = remove_template_terminal_gaps(candidate,template)
        expected = (DNA.makeSequence('ACCGTTG',Name='c 1..7'),
                    DNA.makeSequence('ACCGT-G',Name='t'))
        self.assertEqual(actual[0].Name,expected[0].Name)
        self.assertEqual(actual[1].Name,expected[1].Name)
        self.assertEqual(actual,expected)
        
        # 5' and 3' template terminal gaps
        candidate = DNA.makeSequence('ACCGTTGGA',Name='c')
        template  = DNA.makeSequence('--CGT-GG-',Name='t')
        actual = remove_template_terminal_gaps(candidate,template)
        expected = (DNA.makeSequence('CGTTGG',Name='c 3..8'),
                    DNA.makeSequence('CGT-GG',Name='t'))
        self.assertEqual(actual[0].Name,expected[0].Name)
        self.assertEqual(actual[1].Name,expected[1].Name)
        self.assertEqual(actual,expected)
        
        # name constructed correctly when contains RC
        candidate = DNA.makeSequence('ACCGTTGGA',Name='c RC')
        template  = DNA.makeSequence('--CGT-GG-',Name='t')
        actual = remove_template_terminal_gaps(candidate,template)
        expected = (DNA.makeSequence('CGTTGG',Name='c RC:3..8'),
                    DNA.makeSequence('CGT-GG',Name='t'))
        self.assertEqual(actual[0].Name,expected[0].Name)
        self.assertEqual(actual[1].Name,expected[1].Name)
        self.assertEqual(actual,expected)
        
        # ValueError on unaligned seqs
        candidate = DNA.makeSequence('ACCGTTGGA',Name='c')
        template  = DNA.makeSequence('-CGT-GG-',Name='ct')
        self.assertRaises(ValueError,\
         remove_template_terminal_gaps,candidate,template)

    def test_pynast_seq_3037(self):
        """ uclust as pairwise aligner fixes problematic bl2seq alignment
          
             Strange alignment issues were found with this sequence in 
              PyNAST 1.0. This tests that a good alignment is achieved
              with this seqeunce in later versions.
        """
        template_alignment = LoadSeqs(data=template_128453.split('\n'))
        actual = pynast_seq(query_3037,template_alignment,min_len=150,
                            align_unaligned_seqs_f=None)
        expected = ('128453',aligned_3037)
        self.assertEqual(actual,expected)
        
        
query_3037 = DNA.makeSequence("CTGGGCCGTGTCTCAGTCCCAGTGTGGCTGATCATCCTCTCAGACCAGCTAAGGATCGTCGCCTTGGTGCGCCTTTACCACACCAACTAGCTAAAGGCGATAAATCTTTGATCTCGCGATATCATCCGGTATTAGCAGCAATTTCTCGCTGTTATTCCGAACCTGAGGGCAGATTCCCACGCGTTACGCACCCGTGCGCCACTAAGGCCG",Name=">v15D30.1.08_100583")

aligned_3037 = DNA.makeSequence("----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------C-G----------------------------------------------------------------------------------GC-------------------------------CT--T--AG-T-GG-C-GC-A--C-------------GGG-TGCGT-A--AC-GC-G-T-G-GG---A-A--T-CT-G--C-C-CTC--AG-G------------------------------------------------------------------T-TC----GGA-AT-AA-CAG-------------------------C-G-A-----------------------GAA-A---TTG-CTG-CTAA-TA---CC-G--G-AT-G----------A--------------------T-------------------------------------AT-C-----------------------------------------------------------------------------------------------------------------------G-CG-A--------------------------------------------------------------------------------------------------------------------------------------G-A-T---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CAAA--G-A----------------------------------------------------------------------------------------------------------------------------------------TTT-A----------------------------------------------------------------------------------------------------------------------------------T---C-G--------------C----C-T--------------------------------------------------TT--A--G-CT-A----G---TTGG-T-G-TG-G-T----AAA-GG-C-G-C-ACCA--A-GG-C-G--A-CG-A------------TCC-T-T------AG-CT-G-G-TCT-G-AG----A--GG-AT--G-AT-C-AG-CCAC-A-CTGGG--A-C-TG-A-GA-C-AC-G-G-CCCAG------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------",Name="v15D30.1.08_100583 1..210")

template_128453 = """>128453
------------------------------------------------------------------------------------------------------AACTTGAGAGTTT-GA--T-TC-T-G-GCTC-AG-AA-CGAA-C-GC--TGG-C--G-GC-A-TG--C----T-T--AACACA-T-GC-A-AGT-CGA-A-CGA---------A-G------------------------------------------GC----------------------------------------------------TTC-G----------------------------------------------------------------------------------GC-------------------------------CT--T--AG-T-GG-C-GC-A--C-------------GGG-TGCGT-A--AC-GC-G-T-G-GG---A-A--T-CT-G--C-C-TTC--AG-G------------------------------------------------------------------T-AC----GGA-AT-AA-CTA-------------------------G-G-G-----------------------GAA-A---CTC-GAG-CTAA-TA---CC-G--T-AT-G----------A--------------------T-------------------------------------AT-C-----------------------------------------------------------------------------------------------------------------------G-AG-A--------------------------------------------------------------------------------------------------------------------------------------G-A-T---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CAAA--G-A----------------------------------------------------------------------------------------------------------------------------------------TTT-A----------------------------------------------------------------------------------------------------------------------------------T---C-G--------------C----C-T---G-AA-G---AT---G-A-----G-CCC-GCG--T-TGG--A------TT--A--G-CT-A----G---TTGG-T-A-GG-G-T----AAA-GG-C-T-T-ACCA--A-GG-C-G--A-CG-A------------TCC-A-T------AG-CT-G-G-TCT-G-AG----A--GG-AT--G-AT-C-AG-CCAC-A-CTGGG--A-C-TG-A-GA-C-AC-G-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TG---GG-G-A-ATA-TTGGA-C-AA-T-GG--GG-GA-A----A-C-CC-T-GA-TC-CA-GCAA-TGCC-G-CG-T---G-A-G--T--GA-A-G--A--A-G-G-CC-----TT-AG---------G-G-T-T-G-T--A---AA-G-CTC--------TT-TT-A-C--C-CGG----GA-T--G---A-----------------------T--AA------------------------------T-GA-CA-GT-A-C-CG-G-GA-G---------AA-----------TAAGC-TCC-GG-C-TAA---C--T-CCGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GG-AG-GGA-GCT-A-G-CG-TTGT-T-CGG-AA-TT-A--C-T--GGGC-GTA----AA-GCGT-AC--G-TA-G-G-C-G------------G--T-TT-A-A-T-AA----G-T-C-A---G-GGG-TG-A-AA-GC--CC-AGA-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CT-C-T-GG-AA-C----T-G-C-C-T-T--------T--GA-G-A-C-T-G-TTA--G-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-A-G-A-A-C-A-----T-AG--AA-G-A------------G-GT-A-AG-T----GG--AATT-CCG-A-GT--GT-A-GAG-GTGAAA-TT-CGT-AGAT-A-TT-C-GGA--AG-A-AC-A-CC-AG--T--G--GC-GAA-G--G-C---G----A--C-T-TACTG------G-TC-TA--------------------------------------------------------------TA-G-T-T--GA--CG-----CT-GA-GG--T-A-CGA--AA-G-C--------------G-TGGG-TAG-C-A-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCG-T-AAA--C-GATG-AT--AA-CT---------A-GC--T--G-T-CC-G-GG-T--A--------------------------------------------------------------------------------------CAT-GG--------------------------------------------------------------------------------------------------------------------------------------------------T-A-T-CT--G-G-G-T-GG-C------GG--A----GC-TAA--CG-C-A-T--T--AA-GT--T----A-TCC-GCC-T-G-GG-GAG-TA---CGG-----T-C--G-C-A-A-GAT-T--AAA-ACTC-AAA---------GAAA-TTG-ACGGG-G-G-CCTG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-AAG-CAAC-G-CG-C-AG-A-A-CC-TT-A-CC-AGCGT-TT-G-AC-A-T-C-------------CTGA-T-C-------------G-CG-G-AAA--GT--G-GA-G-A-C--A-C-A-TT-C-T-T--T-C-----AG-------------------------------------T--TC-GG-----------------------------------------CT----G--------GA-TCA-G-A--GA---------------------------------------------------C-A-G-G-T-GCTG-CA-TGG-CT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-CA--CC--T-CTAG--T-T-G-C-C---AT-C-A--T----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TAAG----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T----T-G------------G----G---C-A--CT---------------T-T-A-G-A-GG-A--AC-T-G-CCG--G-T------------------------------------G-A---TAA----------------------------------G-C-C-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAGT-C---CTC-A-T-G-G-C-C-CTT----AC-G--CG-C-T-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-CGGT-G-A--C-AGA-GG-GC--------------------------------------------------------------------------------------------------C-G-C-A-A--G-CCTG-C--A---------------------------------------A-AG-G-T-----------T--A-G-CT---A----------A--TCT-C--------A-AAAAG-CC-G-T-C-T-CAG-TTC--------GGA-T-TGTTC-TC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-AGAGC-A-T-G-AA-G-GC-GGAAT-CG-C-TA--G-TA-AT-C-G-C----GGA-TC-A-G-C-------AT--GCC-GC-G-GT-G-AAT-ACGT-T-CCCAGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-TG-GG-A--G---TTG-G-AT-TC-ACC--C-GAA------G--G-CGC-TG-C-G-C-T-AA-C-C-C-----------------------------------------------------------G-CA-A---------------------------------------------------------------------------------------------------G--GG-A--GG-C--A---GG-CGA--CC--ACG-G----T-GGG-TT-TAG------------------------CG--ACT-GGGG-TG-AAG-TCGTAACAA-GGTAG-CCGT-AGGGGAA-CCTG-CGGC-TGGATCACCTCCTTTCTAAGGA---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"""

db_aln2 = LoadSeqs(data=dict([
('1','ACGT--ACGTAC-ATA-C-----CC-T-G-GTA-G-T---'),
('2','AGGTTTACGTAG-ATA-C-----CC-T-G-GTA-G-T---'),\
('3','AGGTACT-CCAC-ATA-C-----CC-T-G-GTA-G-T---'),
('4','TCGTTCGT-----ATA-C-----CC-T-G-GTA-G-T---'),
('5','ACGTACGT-TA--ATA-C-----CC-T-G-GTA-G-T---')]),\
moltype=DNA,aligned=DenseAlignment)

template_14990_trimmed = """>14990_5_and_3_prime_lost_four_bases_each
--------------------------------------------------------------------------------------------------------------------------------------AG-GA-CGAA-C-GC--TGG-C--G-GC-G-TG--C----C-T--AATACA-T-GC-A-AGT-CGA-G-CGG---------A-A---ATTTTA--------------------------TTGG---TG----------------------------------------------------CTT-G----------------------------------------------------------------------------------CAC-CTT-------------------TAAAAT-TT--T--AG-C-GGCG-G--A--C-------------GGG-TGAGT-A--AC-AC-G-T-G-GG---TAA--C-CTAC--C-T--TA--TA-G------------------------------------------------------------------A-TT----GGG-AT-AA-CTC-------------------------C-G-G-----------------------GAA-A---CCG-GGG-CTAATAC---CG-A----AT-A---------------------------------A-TA-C-T--T--T----------------TTA---AC-------------------------------------------------------------------------------------------------------------------------A-CA-T--------------------------------------------------------------------------------------------------------------------------------------G-T-T--TGA---------------A--A---G-T-T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA--G-A-C-GG-----T-----T-----------------------------------------------------------------------------------------------------------------------TCG--------------------------------------------------------------------------------------------------------------------------G--C--TG--T---C-A--------------C----T-A---T-AA-G---AT---G-G-----A-CCC-GCG--G-CGC--A------TT--A--G-CT-A----G---TTGG-T-G-AG-G-T----AAC-GG-C-T-C-ACCA--A-GG-C-A--A-CG-A------------TGC-G-T------AG-CC-G-A-CCT-G-AG----A--GG-GT--G-AT-C-GG-CCAC-A-CTGGG--A-C-TG-A-GA-C-AC-G-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TA---GG-G-A-ATC-TTCCA-C-AA-T-GG--AC-GA-A----A-G-TC-T-GA-TG-GA-GCAA-CGCC-G-CG-T---G-A-G--T--GA-A-G--A--A-G-G-AT-----TT-CG---------G-T-T-C-G-T--A---AA-A-CTC--------TG-TT-G-C--A-AGG----GA-A--G---AACAAGT---AGCG-TA----G--T--AA-C---T----G-----G--C-GCT-ACC-TT-GA-CG-GT-A-C-CT-T-GT-T---------AG-----------AAAGC-CAC-GG-C-TAA---C--T-ACGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GT-AG-GTG-GCA-A-G-CG-TTGT-C-CGG-AA-TT-A--T-T--GGGC-GTA----AA-GCGC-GC--G-CA-G-G-T-G------------G--T-TC-C-T-T-AA----G-T-C-T---G-ATG-TG-A-AA-GC--CC-CCG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CC-G-G-GG-AG------G-GTC-A-T-T--------G--GA-A-A-C-T-G-GGG--A-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-T-G-----C-AG--AA-G-A------------G-GA-T-AG-T----GG--AATT-CCA-A-GT--GT-A-GCG-GTGAAA-TG-CGT-AGAG-A-TT-T-GGA--GG-A-AC-A-CC-AG--T--G--GC-GAA-G--G-C---G----A--C-T-GTCTG------G-TC-TG--------------------------------------------------------------TA-A-C-T--GA--CA-----CT-GA-GG--C-G-CGA--AA-G-C--------------G-TGGG-GAG-C-A-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCG-T-AAA--C-GATG-AG--TG-CT---------A-AG--T--G-T-TG-G-GG-G--G--T------------------------------------------------------------------------------------TT-CC----------------------------------------------------------------------------------------------------------------------------------------------G---C-C-C-CT--C-A-G-T-GC-T------GC--A----GC-TAA--CG-C-A-T--T--AA-GC--A----C-TCC-GCC-T-G-GG-GAG-TA---CGG-----T-C--G-C-A-A-GAC-T--GAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-AAG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-AGGTC-TT-G-AC-A-TCC--------------CGG-T-G-------------A-CC-A-C-T--AT--G-GA-G-A-C--A-T-A--G-T-T-T--C-C-----CC-------------------------------------T--TC-G------------------------------------------GG----G----G--CAA-CGG---T--GA---------------------------------------------------C-A-G-G-T-GGTG-CA-TGG-TT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TA--TT--C-TTAG--T-T-G-C-C---AT-C-A--T----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TCAG----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T----T-G------------G----G---C-A--CT---------------C-T-A-A-G-GA-G--AC-T-G-CCG--G-T------------------------------------G-A---TAA----------------------------------A-C-C-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAAT-C---ATC-A-T-G-C-C-C-CTT----AT-G--AC-C-T-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-ACGG-T-A--C-AAA-CG-GT--------------------------------------------------------------------------------------------------T-G-C-C-A--A-CCCG-C--G---------------------------------------A-GG-G-G-----------G--A-G-CT---A----------A--TCC-G------A-T-AAAAC-CG-T-T-C-T-CAG-TTC--------GGA-T-TGTAG-GC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-CCTAC-A-T-G-AA-G-CC-GGAAT-CG-C-TA--G-TA-AT-C-G-C----GGA-TC-A-G-C-------AT--GCC-GC-G-GT-G-AAT-ACGT-T-CCCGGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-CG-AG-A--G---TTT-G-TA-AC-ACC--C-GAA------G--T-CGG-TG-A-G-G-T-AA-C-C-T-----------------------------------------------------------T-TA-----------------------------------------------------------------------------------------------------T--GG-A-C-C-C--A---CC-CGC--CG--AAG-G----T-GGG-AT-AAA------------------------TA--ATT-GGGG-TG-AAT-TCTTAACAA-GGTAC-CCGT-ATCGGAA-GGTG-CGGC-TGG------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

input_seq_14990 = """>14990
GCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGGAAATTTTATTGGTGCTTGCACCTTTAAAATTTTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTACCTTATAGATTGGGATAACTCCGGGAAACCGGGGCTAATACCGAATAATACTTTTTAACACATGTTTGAAAGTTGAAAGACGGTTTCGGCTGTCACTATAAGATGGACCCGCGGCGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCAACGATGCGTAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGAAAGTCTGATGGAGCAACGCCGCGTGAGTGAAGAAGGATTTCGGTTCGTAAAACTCTGTTGCAAGGGAAGAACAAGTAGCGTAGTAACTGGCGCTACCTTGACGGTACCTTGTTAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGTGGTTCCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGATAGTGGAATTCCAAGTGTAGCGGTGAAATGCGTAGAGATTTGGAGGAACACCAGTGGCGAAGGCGACTGTCTGGTCTGTAACTGACACTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTGGGGGGTTTCCGCCCCTCAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCGGTGACCACTATGGAGACATAGTTTCCCCTTCGGGGGCAACGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATTCTTAGTTGCCATCATTCAGTTGGGCACTCTAAGGAGACTGCCGGTGATAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGACGGTACAAACGGTTGCCAACCCGCGAGGGGGAGCTAATCCGATAAAACCGTTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTCGGTGAGGTAACCTTTATGGACCCACCCGCCGAAGGTGGGATAAATAATTGGGGTGAATTCTTAACAAGGTACCCGTATCGGAAGGTGCGGCTGGATCA"""

input_seq_10116 = """>10116
CTGGTCCGTGTCTCAGTACCAGTGTGGGGGACCTTCCTCTCAGAACCCCTACGCATCGTCGCCTTGGTGGGCCGTTACCCCACCAACTATCTAATCAGACGCGAGCCCATCTCTGAGCGAATTTCTTTGATATTCAAATCATGCGATTTAAATATGTTATGAGGTATTACCATCCGTTTCCAGAAGCTATCCCTCTCTCAGAGGCAGGTTGCTCACGTGTTACTCACCCGTTCGCCACTCAACTCTTCATCGGTGAGTGCAAGCACTCGGTGATGAAGAAGTTTCGTTCGACTTGCATGTATTAGGCACGCCGCCAGCGTTCATCCTGAGCCAGGATCAAACTCTG"""

expected_fail1 = [('FAKE1 here is some desc.73602 tag1;tag2, tag3:tag4',\
'AGGCGGCTACCTGGACCAACACTGACACTGAGGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGCGAACTGGATGTTGGGTGCAATTTGGCACGCAGTATCGAAGCTAACGCGTTAAGTTCGCCGCCTGGGGAGTACGGTCGCAAGACTTAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGTATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAACTTTCCATAGATGGATTGGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCCTTAGTTGCCAGCACGTAATGGTGGGAACTCTAAGGAGACCGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTAGGGGACCAGGGCTACACACGTACTACAATGGTAGGGACAGAGGGCTGCAAACCCGCGAGGGCAAGCCAATCCCAGAAACCCTATCTCAGTCCGGATTGGAGTTTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGCAGATCAGCATTGCTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTTGTTGCACCAGAAGCAGGTAGCTTAACCTTCGGGAGGGCGCTCACGGTGTGGCCGATGACTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTTTGAGCATGACGTCATCGTCCTGTCGGGCGTCCTCACAAATTACCTGCATTCAGAGATGCGTATCGGCACAGGCCGGTATGCGAAAGTCCCATCATGGGGCCTTAGCTCAGCTGGGAGAGCACCTGCTTTGCAAGCAGGGGGTCGTCGGTTCGATCCCGACAGGCTCCACCATTTGAGTGAAACGACTTTGGGTCTGTAGCTCAGGTGGTTAGAGCGCACCCCTGATAAGGGTGAGGTCGGTGGTTCGAGTCCTCCCAGACCCACCACTCTGAATGTAGTGCACACTTAAGAATTTATATGGCTCAGCGTTGAGGCTGAGACATGTTCTTTTATAACTTGTGACGTAGCGAGCGTTTGAGATATCTATCTAAACGTGTCGTTGAGGCTAAGGCGGGGACTTCGAGTCCCTAAATAATTGAGTCGTATGTTCGCGTTGGTGGCTTTGTACCCCACACAACACGGCGTATGGCCCCGAGGCAACTTGGGGTTATATGGTCAAGCGAATAAGCGCACACGGTGGATGCCTAGGCGGTCAGAGGCGATGAAGGACGTGGTAGCCTGCGAAAAGTGTCGGGGAGCTGGCAACAAGCTTTGATCCGGCAATATCCGAATGGGGAAACCCGG')]


input_seqs2 = """>AKIW1129_fasta.screen.Contig1 description field
GAGTTTGATCATGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAATGACAGAGGAGCTTGCTCCTCTCGATTTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTTATAGCTTGGGATAACTCCGGGAAACCGGGGCTAATACCGAATAATACTTTTGGACACATGTTCGAAAGTTGAAAGATGGTTCTGCTATCACTATAAGATGGACCCGCGCTGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCCACGATGCATAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGAAAGTCTGATGGAGCAACGCCGCGTGAGTGAAGAAGGATTTCGGTTCGTAAAACTCTGTTGTAAGGGAAGAACAAGTACAGTAGTAACTGGCTGTACCTTGACGGTACCTTATTAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGTGGTCCTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGGGACTTGAGTGCAGAAGAGGATAGTGGAATTCCAAGTGTAGCGGTGAAATGCGTAGAGATTTGGAGGAACACCAGTGGCGAAGGCGACTGTCTGGTCTGTAACTGACACTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTGGGGGGTTTCCGCCCCTCAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCATTGACCACTGTAGAGATACAGTTTTCCCTTCGGGGACAACGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTTCTTAGTTGCCATCATTTAGTTGGGCACTCTAAGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGACGGTACAAACGGTTGCCAACCCGCGAGGGGGAGCTAATCCGATAAAACCGTTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTCGGTGAGGTAACCTTTTGGAGCCAGCCGCCGAAGGTGGGATAGATGATTGGGGTGAAGTCGTAACAAGGT"""

input_seqs_gaps = """>FAKE1 here is some desc.73602 tag1;tag2, tag3:tag4
AGGCGGCTACCTGGACCAACACTGACACTGAGGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGCGAACTGGATGTTGGGTGCAATTTGGCACGCAGTATCGAAGCTAACGCGTTAAGTTCGCCGCCTGGGGA
GTACGGTCGCAAGACTTAAA-----CTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGTATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAACTTTCCATAGATGGATTGGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCC
TTGTCCTTAGTTGCCAGCACGTAAT---------GGTGGGAACTCTAAGGAGACCGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTAGGGGACCAGGGCTACACACGTACTACAATGGTA-GGGACAGAGGGCTGCAAACCCGCGAGGGCAAGCCAATCCCAGAAACCCTATCTCAGTCCGGATTGGAGTTTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGCAGATCAGCATTGCTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTTGTTGCACCAGAA
GCAGGTAGCTTAACCTTCGGGAGGGCGCTCACGGTGTGGCCGATGACTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTTTGAGCATGACGTCATCGTCCTGTCGGGCGTCCTCACAAATTACCTGCATTCAGAGATGCGTATCGGCACAGGCCGGTATGCGAAAGTCCCATCATGGGGCCTTAGCTCAGCTGGGAGAGCACCTGCTTTGCAAGCAGGGGGTCGTCGGTTCGATCCCGACAGGCTCCACCATTTGAGTGAAACGACTTTGGGTCTGTAGCTCAGGTGGTTAGAGCGCACCCCTGATAAGGGTGAGGTCGGTGGTTCGAGTCCTCCC-----------------AGACCCACCACTCTGAATGTAGTGCACACTTAAGAATTTATATGGCTCAGCGTTGAGGCTGAGACATGTTCTTTTATAACTTGTGACGTAGCGAGCGTTTGAGATATCTATCTAAACGTGTCGTTGAGGCTAAGGCGGGGACTTCGAGTCCCTAAATAATTGAGTCGTATGTTCGCGTTGGTGGCTTTGTACCCCACACAACACGGCGTATGGCCCCG--AGGCAACTTGGGGT
TATATGGTCAAGCGAATAAGCGCACACGGTGGATGCCTAGGCGGTCAGAGGCGA----TGAAGGACGTGGTAGCCTGCGAAAAGTGTCGGGGAGCTGGCAACAAGCTTTGATCCGGCAATATCCGAATGGGGAAACCCGG
>AKIW1129_fasta.screen.Contig1 description field
GAGTTTGATCATGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAG-------------CGAATGACAGAGGAGCTTGCTCCTCTCGATTTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTTATAGCTTGGGATAACTCCGGGAAACCGGGGCTAATAC-CGAATAATACTTTTGGACACATGTTCGAAAGTTGAAAGATGGTTCTGCTATCACTATAAGATGGACCCGCGCTGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGG-CCACGATGCATAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGAAAGT------CTGATGGAGCAACGCCGCGTGAGTGAAGAAGGATTTCGGTTCGTAAAACTCTGTTGTAAGGGAAGAACAAGTACAGTAGTAACTGGCTGTACCTTGACGGTACCTTATTAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGTGGTCCTT--------TAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGGGACTTGAGTGCAGAAGAGGATAGTGGAATTCCAAGTGTAGCGGTGAAATGCGTAGAGATTTGGAGGAACACCAGTGGCGAAGGCGACTGTCTGGTCTGTAACTGACACTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTGGGGGGTTTCCGCCCCTCAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCAT-------------TGACCACTGTAGAGATACAGTTTTCCCTTCGGGGACAACGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGC----AACCCTTGTTCTTAGTTGCCATCATTTAGTTGGGCACTCTAAGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCT-GGGCTACACACGTGCTACAATGGACGGTACAAACGGTTGCCAACCCGCGAGGGGGAGCTAATCCGATAAAACCGTTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTAC-------ATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCG---------AAGTCGGTGAGGTAACCTTTTGGAGCCAGCCGCCGAAGGTGGGATAGATGATTGGGGTGAAGTCGTAACAAGGT
"""

expected_logfile_contents = \
"""1\t23\t\t5\t100.00\t23
2\t2\tNo search results.
"""

expected_stringent_logfile_contents = \
"""1\t23\tNo search results.
"""

pynast_test_template_fasta1 = """>128618
----------------------------------------------------------------------------------------------------------GGAGAGTTT-GA--T-CC-T-G-GCTC-AG-GA-CGAA-C-GC--TGG-C--G-GC-G-TG--C----C-T--AATACA-T-GC-A-AGT-CGA-G-CGG---------A-C---CG-A----------------------------CGGG---AG----------------------------------------------------CTT-G----------------------------------------------------------------------------------CTC-TCT-------------------TA--G--GT--C--AG-C-GGCG-G--A--C-------------GGG-TGAGT-A--AC-AC-G-T-G-GG---TAA--C-CTGC--C-T--GT--AA-G------------------------------------------------------------------A-CT----GGG-AT-AA-CTC-------------------------C-G-G-----------------------GAA-A---CCG-GGG-CTAATAC---CG-G----AT-G---------------------------------C-TT-G-A--T--T----------------GAA---CC-------------------------------------------------------------------------------------------------------------------------G-CA-T--------------------------------------------------------------------------------------------------------------------------------------G-G-T--TCC---------------A--A--TC-A-T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------AAAA--G-G-T-GG-----C-----T----------------------------------------------------------------------------------------------------------------------TTCA--------------------------------------------------------------------------------------------------------------------------G--C--TA--C---C-A--------------C----T-T---A-CA-G---AT---G-G-----A-CCC-GCG--G-CGC--A------TT--A--G-CT-A----G---TTGG-T-G-AG-G-T----AAC-GG-C-T-C-ACCA--A-GG-C-G--A-CG-A------------TGC-G-T------AG-CC-G-A-CCT-G-AG----A--GG-GT--G-AT-C-GG-CCAC-A-CTGGG--A-C-TG-A-GA-C-AC-G-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TA---GG-G-A-ATC-TTCCG-C-AA-T-GG--AC-GA-A----A-G-TC-T-GA-CG-GA-GCAA-CGCC-G-CG-T---G-A-G--T--GA-T-G--A--A-G-G-TT-----TT-CG---------G-A-T-C-G-T--A---AA-A-CTC--------TG-TT-G-T--T-AGG----GA-A--G---AACAAGT---ACCG-TT----C--G--AA-T---A----G-----GG-C-GGT-ACC-TT-GA-CG-GT-A-C-CT-A-AC-C---------AG-----------AAAGC-CAC-GG-C-TAA---C--T-ACGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GT-AG-GTG-GCA-A-G-CG-TTGT-C-CGG-AA-TT-A--T-T--GGGC-GTA----AA-GCGC-GC--G-CA-G-G-C-G------------G--T-TT-C-T-T-AA----G-T-C-T---G-ATG-TG-A-AA-GC--CC-CCG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CC-G-G-GG-AG------G-GTC-A-T-T--------G--GA-A-A-C-T-G-GGG--A-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-T-G-----C-AG--AA-G-A------------G-GA-G-AG-T----GG--AATT-CCA-C-GT--GT-A-GCG-GTGAAA-TG-CGT-AGAG-A-TG-T-GGA--GG-A-AC-A-CC-AG--T--G--GC-GAA-G--G-C---G----A--C-T-CTCTG------G-TC-TG--------------------------------------------------------------TA-A-C-T--GA--CG-----CT-GA-GG--C-G-CGA--AA-G-C--------------G-TGGG-GAG-C-G-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCG-T-AAA--C-GATG-AG--TG-CT---------A-AG--T--G-T-TA-G-AG-G--G--T------------------------------------------------------------------------------------TT-CC----------------------------------------------------------------------------------------------------------------------------------------------G---C-C-C-TT--T-A-G-T-GC-T------GC--A----GC-AAA--CG-C-A-T--T--AA-GC--A----C-TCC-GCC-T-G-GG-GAG-TA---CGG-----T-C--G-C-A-A-GAC-T--GAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-AAG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-AGGTC-TT-G-AC-A-T-C--------------CTC-T-G-------------A-CA-A-C-C--CT--A-GA-G-A-T--A-G-G--G-C-T-T--C-C-----CC-------------------------------------T--TC-G------------------------------------------GG----G----G---CA-GAG---T--GA---------------------------------------------------C-A-G-G-T-GGTG-CA-TGG-TT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TG--AT--C-TTAG--T-T-G-C-C---AG-C-A--T----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TCAG----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T----T-G------------G----G---C-A--CT---------------C-T-A-A-G-GT-G--AC-T-G-CCG--G-T------------------------------------G-A---CAA----------------------------------A-C-C-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAAT-C---ATC-A-T-G-C-C-C-CTT----AT-G--AC-C-T-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-GCAG-A-A--C-AAA-GG-GC--------------------------------------------------------------------------------------------------A-G-C-G-A--A-GCCG-C--G---------------------------------------A-GG-C-T-----------A--A-G-CC---A----------A--TCC-C------A-C-AAATC-TG-T-T-C-T-CAG-TTC--------GGA-T-CGCAG-TC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-ACTGC-G-T-G-AA-G-CT-GGAAT-CG-C-TA--G-TA-AT-C-G-C----GGA-TC-A-G-C-------AT--GCC-GC-G-GT-G-AAT-ACGT-T-CCCGGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-CG-AG-A--G---TTT-G-TA-AC-ACC--C-GAA------G--T-CGG-TG-A-G-G-T-AA-C-C-T-----------------------------------------------------------T-TT--------------------------------------------------------------------------------------------------------GG-A-G-C-C--A---GC-CGC--CG--AAG-G----T-GGG-AC-AGA------------------------TG--ATT-GGGG-TG-AAG-TCGTAACAA-GGTAG-CCGT-ATCGGAA-GGTG-CGGC-TGGATCACCTCCTTTCT--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
>81187
---------------------------------------------------------------------------------------------------------------AGAGTTTGAT-CC-T-G-GCTC-AG-AG-TGAA-C-GC--TGG-C--G-GC-A-TG--C----C-T--AACACA-T-GC-A-AGT-CGA-A-CG----------G-TAA-CA-G------------------------------GC-C-CG----------------------------------------------------CAA-G----------------------------------------------------------------------------------GG---T------------------G-CT--G--AC--G--AG-T-GG-C-GG-A--C-------------GGG-TGAGG-A--AC-AC-A-T-C-GG---A-A--T-TT-G--C-C-CAG--AC-G------------------------------------------------------------------T-GG----GGG-AT-AA-CGT-------------------------A-G-G-----------------------GAA-A---CTT-ACG-CTAA-TA---CC-G--C-AT-A----------C--------------------G-------------------------------------TC-C-----------------------------------------------------------------------------------------------------------------------T-AC-G--------------------------------------------------------------------------------------------------------------------------------------G-G-A---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA--G-C-G-GG-----G--GA-T--C--------------------------------------------------------------------------------------------------------------------GCA-A----------------------------------------------------------------------------------------------------------------------A----CC-TC--G---C-G--------------C----G-G---T-TG-G---AT---G-A-----G-CCG-ATG--T-CGG--A------TT--A--G-CT-A----G---TTGG-C-G-GG-G-T----AAG-AG-C-C-C-ACCA--A-GG-C-G--A-CG-A------------TCC-G-T------AG-CT-G-G-TCT-G-AG----A--GG-AT--G-AT-C-AG-CCAC-A-TTGGG--A-C-TG-A-GA-C-AC-G-G-CCCAAA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TG---GG-G-A-ATA-TTGGA-C-AA-T-GG--GG-GC-A----A-C-CC-T-GA-TC-CA-GCAA-TGCC-G-CG-T---G-T-G--T--GA-A-G--A--A-G-G-CC-----TT-CG---------G-G-T-T-G-T--A---AA-G-CAC--------TT-TT-A-T--C-AGG----AA-C--G---AA-ACGC---GCTT-GG----T--G--AA-T---A----G-----CA-G-GTG-AAC--T-GA-CG-GT-A-C-CT-G-AG-G---------AA-----------TAAGC-ACC-GG-C-TAA---C--T-TCGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GA-AG-GGT-GCA-A-G-CG-TTAC-T-CGG-AA-TT-A--C-T--GGGC-GTA----AA-GGGT-GC--G-TA-G-G-T-G------------G--T-TG-T-T-T-AA----G-T-C-T---G-CTG-TG-A-AA-GC--CC-CGG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CC-T-G-GG-AA-T----G-G-C-A-G-T--------G--GA-T-A-C-T-G-GGC--A-G-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-A-G-A-A-T-G-----C-GG--TA-G-A------------G-GG-T-AG-T----GG--AATT-CCC-G-GT--GT-A-GCA-GTGAAA-TG-CGT-AGAG-A-TC-G-GGA--GG-A-AC-A-CC-AG--T--G--GC-GAA-G----C---G----G--C-T-ACCTG------G-AC-CA--------------------------------------------------------------GC-A-T-T--GA--CA-----CT-CA-AG--C-A-CGA--AA-G-C--------------G-TGGG-GAG-C-A-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCC-T-AAA--C-GATG-TC--TA-CT---------A-GT--T--G-T-CG-G-GT-C--T---------------------------------------------------------------------------------------TA-AT--------------------------------------------------------------------------------------------------------------------------------------------------T-G-A-CT--T-G-G-T-AA-C------GC--A----GC-TAA--CG-C-G-T--G--AA-GT--A----G-ACC-GCC-T-G-GG-GAG-TA---CGG-----T-C--G-C-A-A-GAT-T--AAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-A-CCCG----C-A--C-A-A-GCG-GT-G--G--AT-GA-T--GT-GGA-TT-AATT-C-G-ATG-CAAC-G-CG-A-AA-A-A-CC-TT-A-CC-TACC--TT-G-AC-A-T-G--------------GCT-G-G-------------A-AT-C-C-C--GG--A-GA-G-A-T--T-T-G--G-G-A-G--T-GC----TC-------------------------------------G--AA-A------------------------------------------GA---GA----A---CC-AGT---A--CA---------------------------------------------------C-A-G-G-T-GCTG-CA-TGG-CT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TG--TC--A-TTAG--T-T-G-C-T---A--C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------G------------G----G---C-A--CT---------------C-T-A-A-T-GA-G--AC-T-G-CCG--G-T------------------------------------G-A---CAA----------------------------------A-C-C-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAGT-C---CTC-A-T-G-G-C-C-CTT----AT-G--GG-T-A-GG-GC-TT-CAC-ACGTC-A--TA--CAATG---G-TACA-T-A--C-AGA--C-GC--------------------------------------------------------------------------------------------------C-G-C-C-A--A-CCCG-C--G---------------------------------------A-GG-G-G-----------G--A-G-CT---A----------A--TCG-C------A-G-AAAGT-GT-A-T-C-G-TAG-TCC--------GGA-T-TGTAG-TC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-ACTGC-A-T-G-AA-G-TT-GGAAT-CG-C-TA--G-TA-AT-C-G-C----GGA-TC-A-G-C-------AT--GTC-GC-G-GT-G-AAT-ACGT-T-CCCGGGTCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-TG-GG-A--G---CGG-G-TT-TT-ACC--A-GAA------G--T-AGG-TA-G-C-T-T-AA-C-C-------------------------------------------------------------G-CA-A------------------------------------------------------------------------------------------------------GG-A--GG-G--C---GC-TTA--CC--ACG-G----T-AGG-AT-TCG------------------------TG--ACT-GGGGTGAAGTCGTAACAAGGTAAC----C-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
>58677
--------------------------------------------------------------------------------------------------------------------------C--T-G-GCTC-AG-GA-CGAA-C-GC--TGG-C--G-GC-G-TG--C----C-T--AATACA-T-GC-A-AGT-CGA-G-CGG---------A-C---CA-A-------------------------------A-T-CG------------------------------------------------GAGCTTGCT----------------------------------------------------------------------------------CTGG--------------------T-TT--G--GT--C--AG-C-GG-C-GG-A--C-------------GGG-TGAGT-A--AC-AC-G-T-G-GG---CAA--C-CT-G--C-C-CGC--AA-G------------------------------------------------------------------A-CC----GGG-AT-AA-CTC-------------------------C-G-G-----------------------GAA-A---CCG-GAG-CTAA-TA---CC-G--G-AT-A----------A--------------------C-A--C-C-G--A--A-----------------GA---CC-G-----------------------------------------------------------------------------------------------------------------------C-AT-G--------------------------------------------------------------------------------------------------------------------------------------G---T--C-T---------------T--T-G-G-T-T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA--G-G-C-GG-----C-CTTTG-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GC-TG--T---C-A--------------C----T-T---G-CG-G---AT---G-G-----G-CCC-GCG--G-CGC--A------TT--A--G-CT-A----G---TTGG-T-G-AG-G-T----AAC-GG-C-T-C-ACCA--A-GG-C-G--A-CG-A------------TGC-G-T------AG-CC-G-G-CCT-G-AG----A--GG-GT--G-AC-C-GG-CCAC-A-CTGGG--A-C-TG-A-GA-C-AC-G-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TA---GG-G-A-ATC-TTCCG-C-AA-T-GG--GC-GA-A----A-G-CC-T-GA-CG-GA-GCGA-CGCC-G-CG-T---G-A-G--C--GA-A-G--A--A-G-G-CC-----TT-CG---------G-G-T-C-G-T--A---AA-G-CTC--------TG-TT-G-T--G-AGG----GA-C--G---AAGGAGC---GCCG-TT----C--G--AA-G---A----G-----GG-C-GGC-GCG-GT-GA-CG-GT-A-C-CT-C-AC-G---------AG-----------AAAGC-CCC-GG-C-TAA---C--T-ACGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GT-AG-GGG-GCG-A-G-CG-TTGT-C-CGG-AA-TT-A--T-T--GGGC-GTA----AA-GCGC-GC--G-CA-G-G-C-G------------G--T-CC-C-T-T-AA----G-T-C-T---G-ATG-TG-A-AA-GC--CC-ACG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CC-G-T-GG-AG-G----G-T-C-A-T-T--------G--GA-A-A-C-T-G-GGG--G-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-T-G-----C-AG--GA-G-A------------G-GA-G-AG-C----GG--AATT-CCA-C-GT--GT-A-GCG-GTGAAA-TG-CGT-AGAG-A-TG-T-GGA--GG-A-AC-A-CC-AG--T--G--GC-GAA-G--G-C---G----G--C-T-CTCTG------G-CC-TG--------------------------------------------------------------CA-A-C-T--GA--CG-----CT-GA-GG--C-G-CGA--AA-G-C--------------G-TGGG-GAG-C-A-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCG-T-AAA--C-GATG-AG--TG-CT---------A-AG--T--G-T-TA-G-AG-G----------------------------------------------------------------------------------------GGTC-ACAC--------------------------------------------------------------------------------------------------------------------------------------------------C-C-TT--T-A-G-T-GC-T------GC--A----GC-TAA--CG-C-G-A--T--AA-GC--A----C-TCC-GCC-T-G-GG-GAG-TA---CGG-----C-C--G-C-A-A-GGC-T--GAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-AAG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-AGGTC-TT-G-AC-A-T-C--------------CCC-T-G-------------A----C-A-A--CC--CAAG-A-G-A--T-T-G--G-G-C-G--T-TC----CC-----------------------------------CCTT-CG-G------------------------------------------GG---GG----A---CA-GGG---T--GA---------------------------------------------------C-A-G-G-T-GGTG-CA-TGG-TT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-CG--CC--T-CTAG--T-T-G-C-C---AG-C-A--T----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TCAG----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T----T-G------------G----G---C-A--CT---------------C-T-A-G-A-GG-G--AC-T-G-CCG--G-C------------------------------------G-A---CAA----------------------------------G-T-C-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAAT-C---ATC-A-T-G-C-C-C-CTT----AT-G--AC-C-T-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-GCGG-T-A--C-AAA-GG-GC--------------------------------------------------------------------------------------------------T-G-C-G-A--A-CCCG-C--G---------------------------------------A-GG-G-G-----------G--A-G-CG---A----------A--TCC-C------A-A-AAAGC-CG-C-T-C-T-CAG-TTC--------GGA-T-TGCAG-GC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-CCTGC-A-T-G-AA-G-CC-GGAAT-CG-C-TA--G-TA-AT-C-G-C----GGA-TC-A-G-C-------AT--GCC-GC-G-GT-G-AA-TACGT-T-CCCGGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-CG-AG-A--G---CTT-G-CA-AC-ACC--C-GAA------G--T-CGG-TG-A-G-G-C-AA-C-C-C-----------------------------------------------------------G-CA-A---------------------------------------------------------------------------------------------------G--GG-A--GC-C--A---GC-CGC--CG--AAG-G----T-GGG-GC-AAG------------------------TG--ATT-GGGG-TG-AAG-TCGTAACAA-GGTAG-CCGT-ACCGGAA--GTG-CGGCTGGATCACCCTCCTT-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
>14308
-------------------------------------------------------------------------------------------------------------------------------------------------------TGG-C--G-GC-G-TG--C----C-T--AACACA-T-GC-A-AGT-CGC-G-CGA---------G-A---AA-G----------------------------CTGC-T--C----------------------------------------------------TTT-G----------------------------------------------------------------------------------AG--CAGT----------------T--A--G--TA--A--AG-C-GG-C-GG-A--C-------------GGG-TGAGT-A--AC-GC-G-T-G-AG---TAA--T-CT-A--C-C-TTT--AA-G------------------------------------------------------------------T-CT----GAT-AT-AA-CTT-------------------------C-T-C-----------------------GAA-A---GGG-AAG-CTAA-TT---TC-G--G-AT-A---------------------------------T-TA-T-G--C--T----------------GCC---TG-G-----------------------------------------------------------------------------------------------------------------------A-TA-A--------------------------------------------------------------------------------------------------------------------------------------C-C-A--G-G---------------C--T-G-C-A-T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CAAA--G-G-C-GG-----C-----T-----------------------------------------------------------------------------------------------------------------------TTT--------------------------------------------------------------------------------------------------------------------------T--GC-CT--C---C-G--------------C----T-T---T-TA-G---AT---G-T-----G-CTC-GCG--T-CCC--A------TT--A--G-CT-T----G---TTGG-T-G-AG-A-T----AAC-AG-C-T-C-ACCA--A-GG-C-T--G-CG-A------------TGG-G-T------AG-CC-G-A-CCT-G-AG----A--GG-GT--G-AT-C-GG-CCAC-A-CTGGG--A-C-TG-A-GA-C-AC-G-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-T-GC-A-G-TG---GG-G-A-ATC-TTTCG-C-AA-T-GA--GC-GC-A----A-G-CT-T-GA-CG-AA-GCGA-CGCC-G-CG-T---G-A-G--T--GA-T-G--A--A-G-G-CC-----TT-CG---------G-G-T-C-G-T--A---AA-G-CTC--------TG-TC-C-T--C-AGG----GA-A--G---AACATCT---TAGT-AG----T--G--AA-T--------A-----AC-T-GCT-AGGCTT-GA-CG-GT-A-C-CT-G-AG-A---------AG-----------AAAGC-TCC-GG-C-TAA---C--T-ACGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GT-AG-GGG-GCA-A-G-CG-TTGT-C-CGG-AA-TC-A--T-T--GGGC-GTA----AA-GGGT-GC--G-CA-G-G-C-G------------G--T-CT-G-G-C-AA----G-T-C-A---A-GTG-TG-A-AA-TG--TA-TCG-G--------------------------------------------------------------------CT-T-AA-------------------------------------------------------------------------CT-G-A-TA-CA------C-TGC-G-C-T--------T--GA-A-A-C-T-G-TCA--G-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-G-G-----C-AA--GA-G-A------------A-GA-G-AG-C----GG--AATT-CCT-A-GT--GT-A-GCG-GTGAAA-TG-CGT-AGAT-A-TT-A-GGA--AG-A-AC-A-CC-AG--T--G--GC-GAA-A--G-C---G----G--C-T-CTCTG------G-CT-TG--------------------------------------------------------------AC-C-C-T--GA--CG-----CT-GA-GG--C-A-CGA--AA-G-C--------------T-AGGG-GAG-C-G-AACG--GG-ATTA-G-ATA-C-----CC-C-G-GTA-G-T----C-CT--G-G-CTG-T-AAA--C-GCTG-GA--TA-CT---------A-GG--T--G-T-TG-G--G-G--G--T------------------------------------------------------------------------------------TC-AA----------------------------------------------------------------------------------------------------------------------------------------------C---T-C-C-CT--C-A-G-T-GC-T------GC--A----GT-TAA--CG-C-G-T--T--AA-GT--A----T-CCC-GCC-T-G-GG-GAT-TA---CGA-----C-C--G-C-A-A-GGT-T--GAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-GCCT-G----C-A--C-A-A-GCG-GC-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-AAG-CAAC-G-CG-C-AG-A-A-CC-TT-A-CC-AGGGC-TT-G-AC-A-T-C------------CCGTGAC-T-------------A-TC-T-G-T--CA--A-CA-G-C-A--G-A-A--T-T-T-G---------GTCC------------------------------------T--TT-G------------------------------------------GA----T----C---AC-ACG-G-T--GA---------------------------------------------------C-A-G-G-T-GGTG-CA-TGG-CT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-C-TA--TC--C-TTAG--T-T-G-C-C---AG-C-A--T----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TAAG----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T----T-G------------G----G---G-A--CT---------------C-T-A-G-G-GA-G--AC-T-G-CCA--G-T------------------------------------C-A---AAA----------------------------------A-C-T-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAGT-C---ATC-A-T-G-C-C-C-CTT----AT-G--CT-C-T-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-CCTG-T-A--C-AGA-GG-GC--------------------------------------------------------------------------------------------------T-G-C-T-A--T-ACCG-C--A---------------------------------------A-GG-T-T-----------T--A-G-CC---A----------A--T-C-C------T-C-AAAAC-AG-G-T-C-C-CAG-TTC--------GGA-T-TGCTG-GC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-CCTGC-A-T-G-AA-G-CT-GGAGT-CG-C-TA--G-TA-AT-C-G-C----GGA-TC-A-G-A-------AT--GCC-GC-G-GT-G-AAT-CCGT-T-CCCAGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-CC-CG-A--G---TTG-G-AT-GC-ACC--A-GAA------G--T-CG---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
>100011
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T--AACACA-T-GC-A-AGT-CGA-A-C-G-----------A---TA-A----------------------------CCTGG--AG----------------------------------------------------CT--G----------------------------------------------------------------------------------CTC-T-A-------------------GG-GA--AT--T--AG-T-GG-C-GA-A--C-------------GGA-GTGAG-T--AC-AC-G-T-G-AG---TAA--C-CT-G--C-C-CTT--GA-C------------------------------------------------------------------T-CT----GGG-AT-AA-CCT-------------------------C-C-G-----------------------GAA-A---CGG-AAG-CTAA------CC-G--G-AT-A---------------------------------T-GA-C-G--C--------------------AC---GGAG-----------------------------------------------------------------------------------------------------------------------G-CA-T-------------------------------------------------------------------------------------------------------------------------------------CT-C----CTG---------------T--G-C-G-T-G-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA--G-A----------------------------------------------------------------------------------------------------------------------------------------ACT------------------------------------------------------------------------------------------------------------------------------------T---C-G--------------G----T-C---A-AG-G---AT---G-G-----A-CTC-GCG--G-CCT--A------TC--A--G-GT-A----G---TTGG-T-G-AG-G-T----AAC-GGCC-C---ACCA--A-GC-C----TACG-A------------CGG-G-T------A--CC-G-G-CCT-G-AG----A--GG-GT--G-AC-C-GG-CCAC-A-CTGGG--A-C-TG-A-TA-C-AC-G-G-CC-AGA-CTCC-TAC-G--G-G---G-GC-A-GC-ACGGTG---GG-G-A-ATA-TTGCA-C-AA-T-GG--GC-GA-A----A-G-CC-T-GA-TG-CA-GCA--CGCC-G-CG-T---A-G-G--G----------A--C-G-G-CC-----TT-CG---------G-G-T-T-G--------AA-C-CT---------TT-TT-A-T--T-AGG----GA-A--G---AAGC------------------------A-A---------------------------GT-GA-CG-GT-A-C-CT-G-TA------------A-----------AAAGC-ACC-GG-C-TAA---C--T-ACGT--GCCA--G-C---A--GCGG-----GG--TA-AT--AC---GT-AG-GGT-GCG-A-G-CG-TTGT-C-CGG-AA-TT-A--T-T--GGGC-GTA----AA-GAGC-TC--G-TA-G-G-C-G------------G--T-CT-G-T-C-GC----G-T-C-T---G-C-G-TGAG-AA-A---AC-CAG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CC-T-C-GG-GC-T----T-G-C-A-G-T--------G--GA-T-A-C-G-G-GCA--G-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-A-G-A-G-T-------C-GG--TA-G-G------------G-GA-G-AA-T----GG--AATT-C---G-GT--GT---GCG-GTGGAA-TG-CGC-AGAT-A-TC-A-GGA--GG-A-CC---CC-GA--T--G--GC-GAA-T--G-C---A----G--T-T-CTCTG------G-C--CG-------------------------------------------------------------TA--A-C-T--GA--CA-----CT-GA-G---A-T-CGA--AA-G-C--------------G-TGGG-A---C-G-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCG-T-AA---C-GTTG-CG--CT--T---------A-GA--T--G-T-GG-G-GA-C--C-------------------------------------------------------------------------------------ATTC-CACG------------------------------------------------------------------------------------------------------------------------------------------------G-T-T--T--C-C-G-T-GT-C------G---A----GC-TAA--CG-C-A-T--T--AA-TG--C----G-CCC-GCC-T-G-GG-GAG-TA---CGG-----C----G-C-A-A-GGC-T--AAA--CTC-AAG------------A-TTG-ACGGG-G-G-CCCG----C-A--C-A-C-GCG-AG-----------A-T--GC-GGA-TT-AATT-G-A-TCG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-AAGGC-TT-G-AC-A-T-A------------C-ACG-A-G-------------A-TA---C-G-GGCCAGAAA-T-G-G----T----C-A-A-C----------TC---------------------------------------TTTGG------------------------------------------AC----------AC-TC-AGT---G--AA---------------------------------------------------C-A-G-G-T-GGTG-CA-TGG-TT--GTC-GTC-A-GC-TC---G-TG-TC-A--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-C-TG--TG--G-TTAG--T-T-G-C-C---AG-C-A--C--G-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TAA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TG---G----T-G------------G----G---A-A--CT---------------C-A-T-A-G-GA-G--AC-T-G-CC---G-G------------------------------------G-T---CAA----------------------------------C-T---G--G-A-GG----AGG-T--GGGG-A-TGAC-GTC--AAAT-A---ATC-A-T-G--CC-C-CTT----AT-G--TC-T-T-GG-GC-TT-CAC-GTATG-C--TA--CAATG---C-CGGT-A-A--C-AAA-GG-GC--------------------------------------------------------------------------------------------------T-G-C-A-A--T-ACCG-T--A---------------------------------------A-GG-T-G-----------G--A---CG---A----------A--TCC-C------A-A-AAA-C-CG-G-T-C-T-CAG-TTC--------GGA-T-TGAGG-TC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-ACCTC-A-T-G--A-G-TC-GGA-T-CG---TA--G-TA-AT-C-G-C----AGA-TC-A----A------AC--GCT--C-G-GT-G--AT-ACGT----CCCGGCCT-TGT-----CACACCG-CCC-GTC-----A---AG--TCA-TG-AA-A--G----TC-G-GA-AC-ACC--C-GA-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

input_seqs1_fasta = """>FAKE1 here is some desc.73602 tag1;tag2, tag3:tag4
AGGCGGCTACCTGGACCAACACTGACACTGAGGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGCGAACTGGATGTTGGGTGCAATTTGGCACGCAGTATCGAAGCTAACGCGTTAAGTTCGCCGCCTGGGGA
GTACGGTCGCAAGACTTAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGTATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAACTTTCCATAGATGGATTGGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCC
TTGTCCTTAGTTGCCAGCACGTAATGGTGGGAACTCTAAGGAGACCGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTAGGGGACCAGGGCTACACACGTACTACAATGGTAGGGACAGAGGGCTGCAAACCCGCGAGGGCAAGCCAATCCCAGAAACCCTATCTCAGTCCGGATTGGAGTTTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGCAGATCAGCATTGCTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTTGTTGCACCAGAA
GCAGGTAGCTTAACCTTCGGGAGGGCGCTCACGGTGTGGCCGATGACTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTTTGAGCATGACGTCATCGTCCTGTCGGGCGTCCTCACAAATTACCTGCATTCAGAGATGCGTATCGGCACAGGCCGGTATGCGAAAGTCCCATCATGGGGCCTTAGCTCAGCTGGGAGAGCACCTGCTTTGCAAGCAGGGGGTCGTCGGTTCGATCCCGACAGGCTCCACCATTTGAGTGAAACGACTTTGGGTCTGTAGCTCAGGTGGTTAGAGCGCACCCCTGATAAGGGTGAGGTCGGTGGTTCGAGTCCTCCCAGACCCACCACTCTGAATGTAGTGCACACTTAAGAATTTATATGGCTCAGCGTTGAGGCTGAGACATGTTCTTTTATAACTTGTGACGTAGCGAGCGTTTGAGATATCTATCTAAACGTGTCGTTGAGGCTAAGGCGGGGACTTCGAGTCCCTAAATAATTGAGTCGTATGTTCGCGTTGGTGGCTTTGTACCCCACACAACACGGCGTATGGCCCCGAGGCAACTTGGGGT
TATATGGTCAAGCGAATAAGCGCACACGGTGGATGCCTAGGCGGTCAGAGGCGATGAAGGACGTGGTAGCCTGCGAAAAGTGTCGGGGAGCTGGCAACAAGCTTTGATCCGGCAATATCCGAATGGGGAAACCCGG
>AKIW1129_fasta.screen.Contig1 description field
GAGTTTGATCATGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAATGACAGAGGAGCTTGCTCCTCTCGATTTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTTATAGCTTGGGATAACTCCGGGAAACCGGGGCTAATACCGAATAATACTTTTGGACACATGTTCGAAAGTTGAAAGATGGTTCTGCTATCACTATAAGATGGACCCGCGCTGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCCACGATGCATAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGAAAGTCTGATGGAGCAACGCCGCGTGAGTGAAGAAGGATTTCGGTTCGTAAAACTCTGTTGTAAGGGAAGAACAAGTACAGTAGTAACTGGCTGTACCTTGACGGTACCTTATTAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGTGGTCCTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGGGACTTGAGTGCAGAAGAGGATAGTGGAATTCCAAGTGTAGCGGTGAAATGCGTAGAGATTTGGAGGAACACCAGTGGCGAAGGCGACTGTCTGGTCTGTAACTGACACTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTGGGGGGTTTCCGCCCCTCAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCATTGACCACTGTAGAGATACAGTTTTCCCTTCGGGGACAACGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTTCTTAGTTGCCATCATTTAGTTGGGCACTCTAAGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGACGGTACAAACGGTTGCCAACCCGCGAGGGGGAGCTAATCCGATAAAACCGTTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTCGGTGAGGTAACCTTTTGGAGCCAGCCGCCGAAGGTGGGATAGATGATTGGGGTGAAGTCGTAACAAGGT
>AKIW521_fasta.screen.Contig1
GAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGCAGCGCGGGGCAACCTGGCGGCGAGTGGCGAACGGGTGAGTAATACATCGGAACGTACCCAGAAGTGGGGGATAACGTAGCGAAAGTTACGCTAATACCGCATACGTTCTACGGAAGAAAGTGGGGGATCTTCGGACCTCATGCTTTTGGAGCGGCCGATGTCTGATTAGCTAGTTGGTGAGGTAAAGGCTCACCAAGGCGACGATCAGTAGCTGGTCTGAGAGGACGACCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATTTTGGACAATGGGCGCAAGCCTGCTCCAGCAATGCCGCGTGAGTGAAGAAGGCCTTCGGGTTGTAAAGCTCTTTTGTCAGGGAAGAAACGGCTGAGGTTAATACCTTCGGCTAATGACGGTACCTGAAGAATAAGCGCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTTTGTAAGTCTGACGTGAAATCCCCGGGCTCAACCTGGGAATTGCGTTGGAGACTGCAAGGCTAGAGTCTGGCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAACACCGATGGGCGAAGGCAGCCCCCTGGGTCAAGACTGACGCTCATGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGTCTACTAGTTGTCGGGTCTTAATTGACTTGGTAACGCAGCTAACGCGTGAAGTAGACCGCCTGGGGAGTACGGTCACAAGATTAAAACTCAAAGGAATTGACGGGGACCCGCACAAGCGGTGGATGATGTGGATTAATTCGATGCAACGCGAAAAACCTTACCTACCCTTGACATGTCAGGAATCCTCGAGAGATTGAGGAGTGCCCGAAAGGGAACCTGAACACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCATTAGTTGCTACGAAAGGGCACTCTAATGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTATGGGTAGGGCTTCACACGTCATACAATGGTACATACAGAGGGCCGCCAACCCGCGAGGGGGAGCTAATCCCAGAAAGTGTATCGTAGTCCGGATCGCAGTCTGCAACTCGACTGCGTGAAGTTGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGTCTTGTACACACCGCCCGTCACACCATGGGAGCGGGTTTTACCAGAAGTAGGTAGCTTAACCGCAAGGGGGGCGCTTACCACGGTAGGATTCGTGACTGGGGTGAAGTCGTAACAAGGTAA
>modified_AKIW1129_both_ends_extended
CCGGAATTCCTTTTAAGAGTTTGATCATGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAATGACAGAGGAGCTTGCTCCTCTCGATTTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTTATAGCTTGGGATAACTCCGGGAAACCGGGGCTAATACCGAATAATACTTTTGGACACATGTTCGAAAGTTGAAAGATGGTTCTGCTATCACTATAAGATGGACCCGCGCTGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCCACGATGCATAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGAAAGTCTGATGGAGCAACGCCGCGTGAGTGAAGAAGGATTTCGGTTCGTAAAACTCTGTTGTAAGGGAAGAACAAGTACAGTAGTAACTGGCTGTACCTTGACGGTACCTTATTAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGTGGTCCTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGGGACTTGAGTGCAGAAGAGGATAGTGGAATTCCAAGTGTAGCGGTGAAATGCGTAGAGATTTGGAGGAACACCAGTGGCGAAGGCGACTGTCTGGTCTGTAACTGACACTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTGGGGGGTTTCCGCCCCTCAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCATTGACCACTGTAGAGATACAGTTTTCCCTTCGGGGACAACGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTTCTTAGTTGCCATCATTTAGTTGGGCACTCTAAGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGACGGTACAAACGGTTGCCAACCCGCGAGGGGGAGCTAATCCGATAAAACCGTTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTCGGTGAGGTAACCTTTTGGAGCCAGCCGCCGAAGGTGGGATAGATGATTGGGGTGAAGTCGTAACAAGGTCCGGAATTCCTTTTAA
>modified_AKIW1129_5_prime_end_extended
CCGGAATTCCTTTTAAGAGTTTGATCATGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAATGACAGAGGAGCTTGCTCCTCTCGATTTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTTATAGCTTGGGATAACTCCGGGAAACCGGGGCTAATACCGAATAATACTTTTGGACACATGTTCGAAAGTTGAAAGATGGTTCTGCTATCACTATAAGATGGACCCGCGCTGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCCACGATGCATAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGAAAGTCTGATGGAGCAACGCCGCGTGAGTGAAGAAGGATTTCGGTTCGTAAAACTCTGTTGTAAGGGAAGAACAAGTACAGTAGTAACTGGCTGTACCTTGACGGTACCTTATTAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGTGGTCCTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGGGACTTGAGTGCAGAAGAGGATAGTGGAATTCCAAGTGTAGCGGTGAAATGCGTAGAGATTTGGAGGAACACCAGTGGCGAAGGCGACTGTCTGGTCTGTAACTGACACTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTGGGGGGTTTCCGCCCCTCAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCATTGACCACTGTAGAGATACAGTTTTCCCTTCGGGGACAACGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTTCTTAGTTGCCATCATTTAGTTGGGCACTCTAAGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGACGGTACAAACGGTTGCCAACCCGCGAGGGGGAGCTAATCCGATAAAACCGTTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTCGGTGAGGTAACCTTTTGGAGCCAGCCGCCGAAGGTGGGATAGATGATTGGGGTGAAGTCGTAACAAGGT
>modified_AKIW1129_3_prime_end_extended
GAGTTTGATCATGGCTCAGGACGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAATGACAGAGGAGCTTGCTCCTCTCGATTTAGCGGCGGACGGGTGAGTAACACGTGGGTAACCTGCCTTATAGCTTGGGATAACTCCGGGAAACCGGGGCTAATACCGAATAATACTTTTGGACACATGTTCGAAAGTTGAAAGATGGTTCTGCTATCACTATAAGATGGACCCGCGCTGCATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCCACGATGCATAGCCGACCTGAGAGGGTGATCGGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTAGGGAATCTTCCACAATGGACGAAAGTCTGATGGAGCAACGCCGCGTGAGTGAAGAAGGATTTCGGTTCGTAAAACTCTGTTGTAAGGGAAGAACAAGTACAGTAGTAACTGGCTGTACCTTGACGGTACCTTATTAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGTGGTCCTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGGGGACTTGAGTGCAGAAGAGGATAGTGGAATTCCAAGTGTAGCGGTGAAATGCGTAGAGATTTGGAGGAACACCAGTGGCGAAGGCGACTGTCTGGTCTGTAACTGACACTGAGGCGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTGGGGGGTTTCCGCCCCTCAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCCATTGACCACTGTAGAGATACAGTTTTCCCTTCGGGGACAACGGTGACAGGTGGTGCATGGTTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTTCTTAGTTGCCATCATTTAGTTGGGCACTCTAAGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGACCTGGGCTACACACGTGCTACAATGGACGGTACAAACGGTTGCCAACCCGCGAGGGGGAGCTAATCCGATAAAACCGTTCTCAGTTCGGATTGTAGGCTGCAACTCGCCTACATGAAGCCGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAGAGTTTGTAACACCCGAAGTCGGTGAGGTAACCTTTTGGAGCCAGCCGCCGAAGGTGGGATAGATGATTGGGGTGAAGTCGTAACAAGGTGATTACACCGGAATTCCTTTTAA"""

input_seqs1_aligned_fasta = """>AKIW1129_fasta.screen.Contig1 description field 1..1507
-------------------------------------------------------------------------------------------------------------GAGTTT-GA--T-CA-T-G-GCTC-AG-GA-CGAA-C-GC--TGG-C--G-GC-G-TG--C----C-T--AATACA-T-GC-A-AGT-CGA-G-CGA---------A-T---GA-C---------------------------AGAGG---AG----------------------------------------------------CTT-G----------------------------------------------------------------------------------CTCCTCT-------------------CG--A--TT--T--AG-C-GGCG-G--A--C-------------GGG-TGAGT-A--AC-AC-G-T-G-GG---TAA--C-CTGC--C-T--TA--TA-G------------------------------------------------------------------C-TT----GGG-AT-AA-CTC-------------------------C-G-G-----------------------GAA-A---CCG-GGG-CTAATAC---CG-A----AT-A---------------------------------A-TA-C-T--T--T----------------TGG---AC-------------------------------------------------------------------------------------------------------------------------A-CA-T--------------------------------------------------------------------------------------------------------------------------------------G---T--TCG---------------A--A--AG-T-T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA--G-A-T-GG----------------------------------------------------------------------------------------------------------------------------------TTCT--------------------------------------------------------------------------------------------------------------------------G--C--TA--T---C-A--------------C----T-A---T-AA-G---AT---G-G-----A-CCC-GCG--C-TGC--A------TT--A--G-CT-A----G---TTGG-T-G-AG-G-T----AAC-GG-C-T-C-ACCA--A-GG-C-C--A-CG-A------------TGC-A-T------AG-CC-G-A-CCT-G-AG----A--GG-GT--G-AT-C-GG-CCAC-A-CTGGG--A-C-TG-A-GA-C-AC-G-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TA---GG-G-A-ATC-TTCCA-C-AA-T-GG--AC-GA-A----A-G-TC-T-GA-TG-GA-GCAA-CGCC-G-CG-T---G-A-G--T--GA-A-G--A--A-G-G-AT-----TT-CG---------G-T-T-C-G-T--A---AA-A-CTC--------TG-TT-G-T--A-AGG----GA-A--G---AACAAGT---ACAG-TA----G--T--AA-C---T----G-----G--C-TGT-ACC-TT-GA-CG-GT-A-C-CT-T-AT-T---------AG-----------AAAGC-CAC-GG-C-TAA---C--T-ACGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GT-AG-GTG-GCA-A-G-CG-TTGT-C-CGG-AA-TT-A--T-T--GGGC-GTA----AA-GCGC-GC--G-CA-G-G-T-G------------G--T-CC-T-T-T-AA----G-T-C-T---G-ATG-TG-A-AA-GC--CC-ACG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CC-G-T-GG-AG------G-GTC-A-T-T--------G--GA-A-A-C-T-G-GGG--G-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-T-G-----C-AG--AA-G-A------------G-GA-T-AG-T----GG--AATT-CCA-A-GT--GT-A-GCG-GTGAAA-TG-CGT-AGAG-A-TT-T-GGA--GG-A-AC-A-CC-AG--T--G--GC-GAA-G--G-C---G----A--C-T-GTCTG------G-TC-TG--------------------------------------------------------------TA-A-C-T--GA--CA-----CT-GA-GG--C-G-CGA--AA-G-C--------------G-TGGG-GAG-C-A-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCG-T-AAA--C-GATG-AG--TG-CT---------A-AG--T--G-T-TG-G-GG-G--G--T------------------------------------------------------------------------------------TT-CC----------------------------------------------------------------------------------------------------------------------------------------------G---C-C-C-CT--C-A-G-T-GC-T------GC--A----GC-TAA--CG-C-A-T--T--AA-GC--A----C-TCC-GCC-T-G-GG-GAG-TA---CGG-----T-C--G-C-A-A-GAC-T--GAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-AAG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-AGGTC-TT-G-AC-A-TCC--------------CAT-T-G-------------A-CC-A-C-T--GT--A-GA-G-A-T--A-C-A--G-T-T-T--T-C-----CC-------------------------------------T--TC-G------------------------------------------GG----G----A--CAA-CGG---T--GA---------------------------------------------------C-A-G-G-T-GGTG-CA-TGG-TT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TG--TT--C-TTAG--T-T-G-C-C---AT-C-A--T----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TTAG----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T----T-G------------G----G---C-A--CT---------------C-T-A-A-G-GA-G--AC-T-G-CCG--G-T------------------------------------G-A---CAA----------------------------------A-C-C-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAAT-C---ATC-A-T-G-C-C-C-CTT----AT-G--AC-C-T-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-ACGG-T-A--C-AAA-CG-GT--------------------------------------------------------------------------------------------------T-G-C-C-A--A-CCCG-C--G---------------------------------------A-GG-G-G-----------G--A-G-CT---A----------A--TCC-G------A-T-AAAAC-CG-T-T-C-T-CAG-TTC--------GGA-T-TGTAG-GC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-CCTAC-A-T-G-AA-G-CC-GGAAT-CG-C-TA--G-TA-AT-C-G-C----GGA-TC-A-G-C-------AT--GCC-GC-G-GT-G-AAT-ACGT-T-CCCGGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-CG-AG-A--G---TTT-G-TA-AC-ACC--C-GAA------G--T-CGG-TG-A-G-G-T-AA-C-C-T-----------------------------------------------------------T-TT--------------------------------------------------------------------------------------------------------GG-A-G-C-C--A---GC-CGC--CG--AAG-G----T-GGG-AT-AGA------------------------TG--ATT-GGGG-TG-AAG-TCGTAACAA-GGT---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
>AKIW521_fasta.screen.Contig1 1..1488
----------------------------------------------------------------------------------------------------------------GAGTTTGAT-CA-T-G-GCTC-AG-AT-TGAA-C-GC--TGG-C--G-GC-A-TG--C----C-T--TACACA-T-GC-A-AGT-CGA-A-CG----------G-CAG-C---------------------------------GC-G-GG----------------------------------------------------GCA-A----------------------------------------------------------------------------------CC---T------------------G-GC--G--GC--G--AG-T-GG-C-GA-A--C-------------GGG-TGAGT-A--AT-AC-A-T-C-GG---A-A--C-GT-A--C-C-CAG--AA-G------------------------------------------------------------------T-GG----GGG-AT-AA-CGT-------------------------A-G-C-----------------------GAA-A---GTT-ACG-CTAA-TA---CC-G--C-AT-A----------C--------------------G-------------------------------------TT-C-----------------------------------------------------------------------------------------------------------------------T-AC-G--------------------------------------------------------------------------------------------------------------------------------------G-A-A---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA--G-T-G-GG-----G--GA-T--C-------------------------------------------------------------------------------------------------------------------TTCG-G----------------------------------------------------------------------------------------------------------------------A----CC-TC--A---T-G--------------C----T-T---T-TG-G---AG---C-G-----G-CCG-ATG--T-CTG--A------TT--A--G-CT-A----G---TTGG-T-G-AG-G-T----AAA-GG-C-T-C-ACCA--A-GG-C-G--A-CG-A------------TCA-G-T------AG-CT-G-G-TCT-G-AG----A--GG-AC--G-AC-C-AG-CCAC-A-CTGGG--A-C-TG-A-GA-C-AC-G-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TG---GG-G-A-ATT-TTGGA-C-AA-T-GG--GC-GC-A----A-G-CC-T-GC-TC-CA-GCAA-TGCC-G-CG-T---G-A-G--T--GA-A-G--A--A-G-G-CC-----TT-CG---------G-G-T-T-G-T--A---AA-G-CTC--------TT-TT-G-T--C-AGG----GA-A--G---AA-ACGG---CTGA-GG----T--T--AA-T---A----C-----CT-T-CGGCTAA--T-GA-CG-GT-A-C-CT-G-AA-G---------AA-----------TAAGC-GCC-GG-C-TAA---C--T-ACGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GT-AG-GGT-GCA-A-G-CG-TTAA-T-CGG-AA-TT-A--C-T--GGGC-GTA----AA-GCGT-GC--G-CA-G-G-C-G------------G--T-TT-T-G-T-AA----G-T-C-T---G-ACG-TG-A-AA-TC--CC-CGG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CC-T-G-GG-AA-T----T-G-C-G-T-T--------G--GA-G-A-C-T-G-CAA--G-G-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-A-G-A-G-T-C-----T-GG--CA-G-A------------G-GG-G-GG-T----AG--AATT-CCA-C-GT--GT-A-GCA-GTGAAA-TG-CGT-AGAG-A-TG-T-GGA--GG-A-AC-A-CC-GA--T--G-GGC-GAA-G---GC---A----G--C-C-CCCTG------G-GT-CA--------------------------------------------------------------AG-A-C-T--GA--CG-----CT-CA-TG--C-A-CGA--AA-G-C--------------G-TGGG-GAG-C-A-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCC-T-AAA--C-GATG-TC--TA-CT---------A-GT--T--G-T-CG-G-GT-C--T---------------------------------------------------------------------------------------TA-AT--------------------------------------------------------------------------------------------------------------------------------------------------T-G-A-CT--T-G-G-T-AA-C------GC--A----GC-TAA--CG-C-G-T--G--AA-GT--A----G-ACC-GCC-T-G-GG-GAG-TA---CGG-----T-C--A-C-A-A-GAT-T--AAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-A-CCCG----C-A--C-A-A-GCG-GT-G--G--AT-GA-T--GT-GGA-TT-AATT-C-G-ATG-CAAC-G-CG-A-AA-A-A-CC-TT-A-CC-TACC-CTT-G-AC-A-T-G--------------TCA-G-G-------------A-AT-C-C-T--CG--A-GA-G-A-T--T-G-A--G-G-A-G--T-GC----CC-------------------------------------G--AA-A------------------------------------------GG---GA----A---CC-TGA---A--CA---------------------------------------------------C-A-G-G-T-GCTG-CA-TGG-CT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TG--TC--A-TTAG--T-T-G-C-T---A--C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------G------------G----G---C-A--CT---------------C-T-A-A-T-GA-G--AC-T-G-CCG--G-T------------------------------------G-A---CAA----------------------------------A-C-C-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAGT-C---CTC-A-T-G-G-C-C-CTT----AT-G--GG-T-A-GG-GC-TT-CAC-ACGTC-A--TA--CAATG---G-TACA-T-A--C-AGA--GGGC--------------------------------------------------------------------------------------------------C-G-C-C-A--A-CCCG-C--G---------------------------------------A-GG-G-G-----------G--A-G-CT---A----------A--TCC-C------A-G-AAAGT-GT-A-T-C-G-TAG-TCC--------GGA-T-CGCAG-TC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-ACTGC-G-T-G-AA-G-TT-GGAAT-CG-C-TA--G-TA-AT-C-G-C----GGA-TC-A-G-C-------AT--GCC-GC-G-GT-G-AAT-ACGT-T-CCCGGGTCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-TG-GG-A--G---CGG-G-TT-TT-ACC--A-GAA------G--T-AGG-TA-G-C-T-T-AA-C-C-------------------------------------------------------------G-CA-A------------------------------------------------------------------------------------------------------GG-G--GG-G--C---GC-TTA--CC--ACG-G----T-AGG-AT-TCG------------------------TG--ACT-GGGGTGAAGTCGTAACAAGGTAA-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
>modified_AKIW1129_both_ends_extended 16..1523
------------------------------------------------------------------------------------------------------------AGAGTTT-GA--T-CA-T-G-GCTC-AG-GA-CGAA-C-GC--TGG-C--G-GC-G-TG--C----C-T--AATACA-T-GC-A-AGT-CGA-G-CGA---------A-T---GA-C---------------------------AGAGG---AG----------------------------------------------------CTT-G----------------------------------------------------------------------------------CTCCTCT-------------------CG--A--TT--T--AG-C-GGCG-G--A--C-------------GGG-TGAGT-A--AC-AC-G-T-G-GG---TAA--C-CTGC--C-T--TA--TA-G------------------------------------------------------------------C-TT----GGG-AT-AA-CTC-------------------------C-G-G-----------------------GAA-A---CCG-GGG-CTAATAC---CG-A----AT-A---------------------------------A-TA-C-T--T--T----------------TGG---AC-------------------------------------------------------------------------------------------------------------------------A-CA-T--------------------------------------------------------------------------------------------------------------------------------------G---T--TCG---------------A--A--AG-T-T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA--G-A-T-GG----------------------------------------------------------------------------------------------------------------------------------TTCT--------------------------------------------------------------------------------------------------------------------------G--C--TA--T---C-A--------------C----T-A---T-AA-G---AT---G-G-----A-CCC-GCG--C-TGC--A------TT--A--G-CT-A----G---TTGG-T-G-AG-G-T----AAC-GG-C-T-C-ACCA--A-GG-C-C--A-CG-A------------TGC-A-T------AG-CC-G-A-CCT-G-AG----A--GG-GT--G-AT-C-GG-CCAC-A-CTGGG--A-C-TG-A-GA-C-AC-G-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TA---GG-G-A-ATC-TTCCA-C-AA-T-GG--AC-GA-A----A-G-TC-T-GA-TG-GA-GCAA-CGCC-G-CG-T---G-A-G--T--GA-A-G--A--A-G-G-AT-----TT-CG---------G-T-T-C-G-T--A---AA-A-CTC--------TG-TT-G-T--A-AGG----GA-A--G---AACAAGT---ACAG-TA----G--T--AA-C---T----G-----G--C-TGT-ACC-TT-GA-CG-GT-A-C-CT-T-AT-T---------AG-----------AAAGC-CAC-GG-C-TAA---C--T-ACGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GT-AG-GTG-GCA-A-G-CG-TTGT-C-CGG-AA-TT-A--T-T--GGGC-GTA----AA-GCGC-GC--G-CA-G-G-T-G------------G--T-CC-T-T-T-AA----G-T-C-T---G-ATG-TG-A-AA-GC--CC-ACG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CC-G-T-GG-AG------G-GTC-A-T-T--------G--GA-A-A-C-T-G-GGG--G-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-T-G-----C-AG--AA-G-A------------G-GA-T-AG-T----GG--AATT-CCA-A-GT--GT-A-GCG-GTGAAA-TG-CGT-AGAG-A-TT-T-GGA--GG-A-AC-A-CC-AG--T--G--GC-GAA-G--G-C---G----A--C-T-GTCTG------G-TC-TG--------------------------------------------------------------TA-A-C-T--GA--CA-----CT-GA-GG--C-G-CGA--AA-G-C--------------G-TGGG-GAG-C-A-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCG-T-AAA--C-GATG-AG--TG-CT---------A-AG--T--G-T-TG-G-GG-G--G--T------------------------------------------------------------------------------------TT-CC----------------------------------------------------------------------------------------------------------------------------------------------G---C-C-C-CT--C-A-G-T-GC-T------GC--A----GC-TAA--CG-C-A-T--T--AA-GC--A----C-TCC-GCC-T-G-GG-GAG-TA---CGG-----T-C--G-C-A-A-GAC-T--GAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-AAG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-AGGTC-TT-G-AC-A-TCC--------------CAT-T-G-------------A-CC-A-C-T--GT--A-GA-G-A-T--A-C-A--G-T-T-T--T-C-----CC-------------------------------------T--TC-G------------------------------------------GG----G----A--CAA-CGG---T--GA---------------------------------------------------C-A-G-G-T-GGTG-CA-TGG-TT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TG--TT--C-TTAG--T-T-G-C-C---AT-C-A--T----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TTAG----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T----T-G------------G----G---C-A--CT---------------C-T-A-A-G-GA-G--AC-T-G-CCG--G-T------------------------------------G-A---CAA----------------------------------A-C-C-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAAT-C---ATC-A-T-G-C-C-C-CTT----AT-G--AC-C-T-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-ACGG-T-A--C-AAA-CG-GT--------------------------------------------------------------------------------------------------T-G-C-C-A--A-CCCG-C--G---------------------------------------A-GG-G-G-----------G--A-G-CT---A----------A--TCC-G------A-T-AAAAC-CG-T-T-C-T-CAG-TTC--------GGA-T-TGTAG-GC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-CCTAC-A-T-G-AA-G-CC-GGAAT-CG-C-TA--G-TA-AT-C-G-C----GGA-TC-A-G-C-------AT--GCC-GC-G-GT-G-AAT-ACGT-T-CCCGGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-CG-AG-A--G---TTT-G-TA-AC-ACC--C-GAA------G--T-CGG-TG-A-G-G-T-AA-C-C-T-----------------------------------------------------------T-TT--------------------------------------------------------------------------------------------------------GG-A-G-C-C--A---GC-CGC--CG--AAG-G----T-GGG-AT-AGA------------------------TG--ATT-GGGG-TG-AAG-TCGTAACAA-GGT---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
>modified_AKIW1129_5_prime_end_extended 16..1523
------------------------------------------------------------------------------------------------------------AGAGTTT-GA--T-CA-T-G-GCTC-AG-GA-CGAA-C-GC--TGG-C--G-GC-G-TG--C----C-T--AATACA-T-GC-A-AGT-CGA-G-CGA---------A-T---GA-C---------------------------AGAGG---AG----------------------------------------------------CTT-G----------------------------------------------------------------------------------CTCCTCT-------------------CG--A--TT--T--AG-C-GGCG-G--A--C-------------GGG-TGAGT-A--AC-AC-G-T-G-GG---TAA--C-CTGC--C-T--TA--TA-G------------------------------------------------------------------C-TT----GGG-AT-AA-CTC-------------------------C-G-G-----------------------GAA-A---CCG-GGG-CTAATAC---CG-A----AT-A---------------------------------A-TA-C-T--T--T----------------TGG---AC-------------------------------------------------------------------------------------------------------------------------A-CA-T--------------------------------------------------------------------------------------------------------------------------------------G---T--TCG---------------A--A--AG-T-T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA--G-A-T-GG----------------------------------------------------------------------------------------------------------------------------------TTCT--------------------------------------------------------------------------------------------------------------------------G--C--TA--T---C-A--------------C----T-A---T-AA-G---AT---G-G-----A-CCC-GCG--C-TGC--A------TT--A--G-CT-A----G---TTGG-T-G-AG-G-T----AAC-GG-C-T-C-ACCA--A-GG-C-C--A-CG-A------------TGC-A-T------AG-CC-G-A-CCT-G-AG----A--GG-GT--G-AT-C-GG-CCAC-A-CTGGG--A-C-TG-A-GA-C-AC-G-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TA---GG-G-A-ATC-TTCCA-C-AA-T-GG--AC-GA-A----A-G-TC-T-GA-TG-GA-GCAA-CGCC-G-CG-T---G-A-G--T--GA-A-G--A--A-G-G-AT-----TT-CG---------G-T-T-C-G-T--A---AA-A-CTC--------TG-TT-G-T--A-AGG----GA-A--G---AACAAGT---ACAG-TA----G--T--AA-C---T----G-----G--C-TGT-ACC-TT-GA-CG-GT-A-C-CT-T-AT-T---------AG-----------AAAGC-CAC-GG-C-TAA---C--T-ACGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GT-AG-GTG-GCA-A-G-CG-TTGT-C-CGG-AA-TT-A--T-T--GGGC-GTA----AA-GCGC-GC--G-CA-G-G-T-G------------G--T-CC-T-T-T-AA----G-T-C-T---G-ATG-TG-A-AA-GC--CC-ACG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CC-G-T-GG-AG------G-GTC-A-T-T--------G--GA-A-A-C-T-G-GGG--G-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-T-G-----C-AG--AA-G-A------------G-GA-T-AG-T----GG--AATT-CCA-A-GT--GT-A-GCG-GTGAAA-TG-CGT-AGAG-A-TT-T-GGA--GG-A-AC-A-CC-AG--T--G--GC-GAA-G--G-C---G----A--C-T-GTCTG------G-TC-TG--------------------------------------------------------------TA-A-C-T--GA--CA-----CT-GA-GG--C-G-CGA--AA-G-C--------------G-TGGG-GAG-C-A-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCG-T-AAA--C-GATG-AG--TG-CT---------A-AG--T--G-T-TG-G-GG-G--G--T------------------------------------------------------------------------------------TT-CC----------------------------------------------------------------------------------------------------------------------------------------------G---C-C-C-CT--C-A-G-T-GC-T------GC--A----GC-TAA--CG-C-A-T--T--AA-GC--A----C-TCC-GCC-T-G-GG-GAG-TA---CGG-----T-C--G-C-A-A-GAC-T--GAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-AAG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-AGGTC-TT-G-AC-A-TCC--------------CAT-T-G-------------A-CC-A-C-T--GT--A-GA-G-A-T--A-C-A--G-T-T-T--T-C-----CC-------------------------------------T--TC-G------------------------------------------GG----G----A--CAA-CGG---T--GA---------------------------------------------------C-A-G-G-T-GGTG-CA-TGG-TT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TG--TT--C-TTAG--T-T-G-C-C---AT-C-A--T----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TTAG----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T----T-G------------G----G---C-A--CT---------------C-T-A-A-G-GA-G--AC-T-G-CCG--G-T------------------------------------G-A---CAA----------------------------------A-C-C-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAAT-C---ATC-A-T-G-C-C-C-CTT----AT-G--AC-C-T-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-ACGG-T-A--C-AAA-CG-GT--------------------------------------------------------------------------------------------------T-G-C-C-A--A-CCCG-C--G---------------------------------------A-GG-G-G-----------G--A-G-CT---A----------A--TCC-G------A-T-AAAAC-CG-T-T-C-T-CAG-TTC--------GGA-T-TGTAG-GC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-CCTAC-A-T-G-AA-G-CC-GGAAT-CG-C-TA--G-TA-AT-C-G-C----GGA-TC-A-G-C-------AT--GCC-GC-G-GT-G-AAT-ACGT-T-CCCGGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-CG-AG-A--G---TTT-G-TA-AC-ACC--C-GAA------G--T-CGG-TG-A-G-G-T-AA-C-C-T-----------------------------------------------------------T-TT--------------------------------------------------------------------------------------------------------GG-A-G-C-C--A---GC-CGC--CG--AAG-G----T-GGG-AT-AGA------------------------TG--ATT-GGGG-TG-AAG-TCGTAACAA-GGT---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
>modified_AKIW1129_3_prime_end_extended 12..1520
-----------------------------------------------------------------------------------------------------------------------------T-G-GCTC-AG-GA-CGAA-C-GC--TGG-C--G-GC-G-TG--C----C-T--AATACA-T-GC-A-AGT-CGA-G-CGA---------A-T---GA-C-------------------------------A-G-AG------------------------------------------------GAGCTTGCT----------------------------------------------------------------------------------CCTC--------------------T-CG--A--TT--T--AG-C-GG-C-GG-A--C-------------GGG-TGAGT-A--AC-AC-G-T-G-GG---TAA--C-CT-G--C-C-TTA--TA-G------------------------------------------------------------------C-TT----GGG-AT-AA-CTC-------------------------C-G-G-----------------------GAA-A---CCG-GGG-CTAA-TA---CC-G--A-AT-A----------A--------------------T-A--C-T-T--T--T-----------------GG---AC-A-----------------------------------------------------------------------------------------------------------------------C-AT-G--------------------------------------------------------------------------------------------------------------------------------------T---T--C-G---------------A--A-A-G-T-T-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GAAA--G-A-T-GG--------TTCT-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------GC-TA--T---C-A--------------C----T-A---T-AA-G---AT---G-G-----A-CCC-GCG--C-TGC--A------TT--A--G-CT-A----G---TTGG-T-G-AG-G-T----AAC-GG-C-T-C-ACCA--A-GG-C-C--A-CG-A------------TGC-A-T------AG-CC-G-A-CCT-G-AG----A--GG-GT--G-AT-C-GG-CCAC-A-CTGGG--A-C-TG-A-GA-C-AC-G-G-CCCAGA-CTCC-TAC-G--G-G-A-G-GC-A-GC-A-G-TA---GG-G-A-ATC-TTCCA-C-AA-T-GG--AC-GA-A----A-G-TC-T-GA-TG-GA-GCAA-CGCC-G-CG-T---G-A-G--T--GA-A-G--A--A-G-G-AT-----TT-CG---------G-T-T-C-G-T--A---AA-A-CTC--------TG-TT-G-T--A-AGG----GA-A--G---AACAAGT---ACAG-TA----G--T--AA-C---T----G-----G--C-TGT-ACC-TT-GA-CG-GT-A-C-CT-T-AT-T---------AG-----------AAAGC-CAC-GG-C-TAA---C--T-ACGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GT-AG-GTG-GCA-A-G-CG-TTGT-C-CGG-AA-TT-A--T-T--GGGC-GTA----AA-GCGC-GC--G-CA-G-G-T-G------------G--T-CC-T-T-T-AA----G-T-C-T---G-ATG-TG-A-AA-GC--CC-ACG-G--------------------------------------------------------------------CT-C-AA-------------------------------------------------------------------------CC-G-T-GG-AG-G----G-T-C-A-T-T--------G--GA-A-A-C-T-G-GGG--G-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-T-G-----C-AG--AA-G-A------------G-GA-T-AG-T----GG--AATT-CCA-A-GT--GT-A-GCG-GTGAAA-TG-CGT-AGAG-A-TT-T-GGA--GG-A-AC-A-CC-AG--T--G--GC-GAA-G--G-C---G----A--C-T-GTCTG------G-TC-TG--------------------------------------------------------------TA-A-C-T--GA--CA-----CT-GA-GG--C-G-CGA--AA-G-C--------------G-TGGG-GAG-C-A-AACA--GG-ATTA-G-ATA-C-----CC-T-G-GTA-G-T----C-CA--C-G-CCG-T-AAA--C-GATG-AG--TG-CT---------A-AG--T--G-T-TG-G-GG-G----------------------------------------------------------------------------------------GTTT-CCGC--------------------------------------------------------------------------------------------------------------------------------------------------C-C-CT--C-A-G-T-GC-T------GC--A----GC-TAA--CG-C-A-T--T--AA-GC--A----C-TCC-GCC-T-G-GG-GAG-TA---CGG-----T-C--G-C-A-A-GAC-T--GAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-AAG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-AGGTC-TT-G-AC-A-T-C--------------CCATT-G-------------A----C-C-A--CT--GTAG-A-G-A--T------A-C-A-G--T-TT----TC-----------------------------------CCTT-CG-G------------------------------------------GG---AC----A----A-CGG---T--GA---------------------------------------------------C-A-G-G-T-GGTG-CA-TGG-TT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TG--TT--C-TTAG--T-T-G-C-C---AT-C-A--T----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TTAG----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T----T-G------------G----G---C-A--CT---------------C-T-A-A-G-GA-G--AC-T-G-CCG--G-T------------------------------------G-A---CAA----------------------------------A-C-C-G--G-A-GG-A--AGG-T--GGGG-A-TGAC-GTC--AAAT-C---ATC-A-T-G-C-C-C-CTT----AT-G--AC-C-T-GG-GC-TA-CAC-ACGTG-C--TA--CAATG---G-ACGG-T-A--C-AAA-CG-GT--------------------------------------------------------------------------------------------------T-G-C-C-A--A-CCCG-C--G---------------------------------------A-GG-G-G-----------G--A-G-CT---A----------A--TCC-G------A-T-AAAAC-CG-T-T-C-T-CAG-TTC--------GGA-T-TGTAG-GC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-CCTAC-A-T-G-AA-G-CC-GGAAT-CG-C-TA--G-TA-AT-C-G-C----GGA-TC-A-G-C-------AT--GCC-GC-G-GT-G-AA-TACGT-T-CCCGGGCCT-TGTA----CACACCG-CCC-GTC-----A---CA--CCA-CG-AG-A--G---TTT-G-TA-AC-ACC--C-GAA------G--T-CGG-TG-A-G-G-T-AA-C-C---------------------------------------------------------------TT-T---------------------------------------------------------------------------------------------------T--GG-A--GC-C--A---GC-CGC--CG--AAG-G----T-GGG-AT-AGA------------------------TG--ATT-GGGG-TG-AAG-TCGTAACAA-GGTGA-TTAC-ACCGGAA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

input_seqs1_fail_fasta = """>FAKE1 here is some desc.73602 tag1;tag2, tag3:tag4
AGGCGGCTACCTGGACCAACACTGACACTGAGGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGCGAACTGGATGTTGGGTGCAATTTGGCACGCAGTATCGAAGCTAACGCGTTAAGTTCGCCGCCTGGGGAGTACGGTCGCAAGACTTAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCGGTGGAGTATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAACTTTCCATAGATGGATTGGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCCTTAGTTGCCAGCACGTAATGGTGGGAACTCTAAGGAGACCGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTAGGGGACCAGGGCTACACACGTACTACAATGGTAGGGACAGAGGGCTGCAAACCCGCGAGGGCAAGCCAATCCCAGAAACCCTATCTCAGTCCGGATTGGAGTTTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGCAGATCAGCATTGCTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTTTGTTGCACCAGAAGCAGGTAGCTTAACCTTCGGGAGGGCGCTCACGGTGTGGCCGATGACTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTTTGAGCATGACGTCATCGTCCTGTCGGGCGTCCTCACAAATTACCTGCATTCAGAGATGCGTATCGGCACAGGCCGGTATGCGAAAGTCCCATCATGGGGCCTTAGCTCAGCTGGGAGAGCACCTGCTTTGCAAGCAGGGGGTCGTCGGTTCGATCCCGACAGGCTCCACCATTTGAGTGAAACGACTTTGGGTCTGTAGCTCAGGTGGTTAGAGCGCACCCCTGATAAGGGTGAGGTCGGTGGTTCGAGTCCTCCCAGACCCACCACTCTGAATGTAGTGCACACTTAAGAATTTATATGGCTCAGCGTTGAGGCTGAGACATGTTCTTTTATAACTTGTGACGTAGCGAGCGTTTGAGATATCTATCTAAACGTGTCGTTGAGGCTAAGGCGGGGACTTCGAGTCCCTAAATAATTGAGTCGTATGTTCGCGTTGGTGGCTTTGTACCCCACACAACACGGCGTATGGCCCCGAGGCAACTTGGGGTTATATGGTCAAGCGAATAAGCGCACACGGTGGATGCCTAGGCGGTCAGAGGCGATGAAGGACGTGGTAGCCTGCGAAAAGTGTCGGGGAGCTGGCAACAAGCTTTGATCCGGCAATATCCGAATGGGGAAACCCGG
"""


input_seqs2_fasta = """>2855189 SLEpi20M_15561395
TACGAAAGATCCAAGCGTTATTCGAAATGATTGGGCNTAAANAGTTTGTAGGCGGTATTTGTACTCACTTCTAAAAAACTAAGATTATCTCTTAGTATGG
"""

pynast_test_template_fasta2 = """>26799
-----------------------------------------------------------------------------------------------------AAATGGAGAGGTTT-GA--T-CC-T-G-GCTC-AG-GA-TGAA-C-GC--TGG-C--G-AT-A-TG--C----T-T--AACACA-T-GC-A-AGT-CGA-A-CGA---------A-T---AT-T--------------------------AAGTTTTCTTAAA--------------------------------------------------TTT-G----------------------------------------------------------------------------------TAG-AAA-------------------TT--TA-AT--ATTAG-T-GG-C-GA-A--C-------------GGG-TGAGT-A--AC-GC-G-T-A-AG---A-A--T-CT-G--C-T-TTT--GG-G------------------------------------------------------------------T-AA----AGA-AT-AA-CAA-------------------------T-T-G-----------------------GAA-A---CGA-TTG-CTAA-TA---CT-T--T-AT-A----------G----------------------------------------------------------GC-T-----------------------------------------------------------------------------------------------------------------------G-AG-G--------------------------------------------------------------------------------------------------------------------------------------A-G-T---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TAAA--G-G-T--------------------------------------------------------------------------------------------------------------------------------------TTT-A-------------------------------------------------------------------------------------------------------------------------------T--T-TCC-G--------------C----C-C---A-GA-A---AT---G-A-----G-CTT-GCG--T-CTG--A------TT--A--G-CT-A----G---TTGG-T-A-AG-A-T----AAA-AG-C-T-T-ACCA--A-GG-C-A--A-TG-A------------TCA-G-T------AG-TT-G-G-TCT-G-AG----A--GG-AT--G-AT-C-AA-CCAC-A-CTGGG--A-C-TG-A-GA-T-AC-G-G-CCCAGA-CCTT-TAC-G--G-A-G-G-GC-A-GC-A-G-TG---AG-G-A-ATT-TTCCG-C-AA-T-GG--GC-GA-A----A-G-CC-T-GA-CG-GA-GCAA-TATC-G-CG-T---G-A-A--G--GA-T-G--A--C-G-G-CC-----TG-TG---------G-G-T-T-G-T--A---AA-C-TTC--------TT-TT-C-T--T-AAG----AA-A--G---A--------------------A--T--TC------------------------------T-GA-CG-GT-A-C-TT-A-AG-G---------AA-----------TAAGC-ATC-GG-C-TAA---C--T-CCGT--GCCA--G-C---A--GCCG---C-GG--TA-AT--AC---GG-AG-GAT-GCA-A-G-CG-TTAT-C-CGA-AA-TT-A--T-T--GGGC-GTA----AA-GAGT-TT--G-TA-G-G-T-G------------G--T-TT-T-T-T-AA----G-T-C-T---A-CTG-TT-A-AA-TA--TC-AGA-G--------------------------------------------------------------------CT-T-AA-------------------------------------------------------------------------CT-T-T-GA-AC-A----A-G-C-A-G-T--------A-TGA-A-A-C-T-A-ATT--A-A-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------T-T-G-A-G-T-T-----T-GG--TA-G-A------------G-GC-A-GA-G----GG--AACT-CTC-G-AT--GT-A-GTG-GTGAAA-TA-CGT-AGAT-A-TC-G-GGG--GG-A-AC-A-CC-AG--T--A--GC-GAA-A--G-C---G----C--T-C-TGCTG------G-GC-CA--------------------------------------------------------------TA-A-C-T--GA--CA-----CT-GA-GA--A-A-CGA--AA-G-C--------------T-AGGG-GAG-C-A-AATA--GG-ATTA-G-ATA-C-----CC-T-A-GTA-G-T----C-CT--A-G-CTG-T-AAA--C-GATG-GA--TA-CT---------A-AG--T--A-T-TG-G-GC------------------------------------------------------------------------------------------TTTTTGAAG------------------------------------------------------------------------------------------------------------------------------------------------------TT--C-A-G-T-GT-T------GA--A----GC-TAA--CG-C-G-T--T--AA-GT--A----T-CCC-GCC-T-G-GG-GAG-TA---CGT-----T-C--G-C-A-A-GAA-T--GAA-ACTC-AAA---------GGAA-TTG-ACGGG-G-G-CCCG----C-A--C-A-A-GCG-GT-G--G--AG-CA-T--GT-GGT-TT-AATT-C-G-ATG-CAAC-G-CG-A-AG-A-A-CC-TT-A-CC-AGGAA-TT-G-AC-A-T-A--------------CTC-G-T--------------TGGTT-T-T--TT--A-GA-A-A-T--A-A-A--A-A-A-------------C-------------------------------------T--GT-T------------------------------------------A--------------AA-GAG---A--TA---------------------------------------------------C-A-G-G-T-GGTG-CA-TGG-CT--GTC-GTC-A-GC-TC---G-TG-TC-G--TGA-GA-TGT-T-GG-G-TT-AA-GT-CCCGC-AA--------C-GAG-CGC-A-ACC-C-T-TG--TC--T-TTAG--T-T-G-T-T---AT-C---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TA---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------G-A-GA-G--AC-T-G-CCG--G-T------------------------------------G-A---TAA----------------------------------A-C-C-G--G-A-GG-A--AGG-T--GAGG-A-TGAC-GTC--AAGT-C---AGC-A-T-G-C-C-C-CTT----AA-G--TC-C-T-GG-GC-GA-CAC-ACGTG-C--TA--CAATG---G-TATA-G-A--C-AAA-GG-GA--------------------------------------------------------------------------------------------------A-G-C-A-A--A-TCTG-C--G---------------------------------------A-AG-A-G-----------T--A-G-CA---A----------A--TCT-C------A---AAAAC-TATA-T-C-T-CAG-TTC--------GGA-T-TGCAG-GC--T-GCAA-CT-C-------------------------------------------------------------------------------------------------G-CCTGC-A-T-G-AA-G-TC-GGAAT-CG-C-TA--G-TA-AT-C-G-C----TGG-TC-A-G-CC------AT--ACA-GC-G-GT-G-AAT-ATGT-T-CTCGGGCCT-TGTA----CACACCG-CCC-ATC-----A---CG--CTC-GA-GA-A--A---TTG-G-AA-AT-ACC--C-AAA------G--T-CAT-CA-T-T-C-T-AA-CCATATT---------------------------------------------------------T-TT-T---------------------------------------------------------------------------------------------------G---G-A--AG-A--T---AA-TGC--CA--AAG-G----T-AGA-GC-TAG------------------------TG--ACT-CAAG-CG-AAG-TTGTAACAA-GGTAA-CCGT-ACTGGAA-GGTG-CGGT-TGGATCACCTCCTTA----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

if __name__ == "__main__":
    main()
