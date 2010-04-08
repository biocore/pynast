#!/usr/bin/env python

from __future__ import division
from os import system, remove, popen
from os.path import exists
from shutil import copy as copy_file
from glob import glob
from cogent import DNA, LoadSeqs
from cogent.util.misc import remove_files
from cogent.core.alignment import SequenceCollection, DenseAlignment
from cogent.align.align import make_dna_scoring_dict, global_pairwise
from cogent.app.blast import blastn
from cogent.app.formatdb import build_blast_db_from_seqs, \
 build_blast_db_from_fasta_path
from cogent.app.muscle import align_unaligned_seqs as muscle_align_unaligned_seqs
from cogent.app.mafft import align_unaligned_seqs as mafft_align_unaligned_seqs
from cogent.app.clustalw import align_unaligned_seqs as clustal_align_unaligned_seqs
from cogent.app.util import get_tmp_filename
from cogent.app.uclust import uclust_search_and_align_from_fasta_filepath
from cogent.parse.blast import BlastResult
from cogent.parse.fasta import MinimalFastaParser
from pynast.logger import NastLogger

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The PyNAST Project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

""" PyNAST is a complete rewrite of the NAST algorithm written in python.
 While PyNAST 1.0 strived to exactly match the results of the original 
 NAST algorithm, the later version (beginning with the post-1.0 development
 code) no longer exactly matches the the original NAST algorithm, hopefully
 in favor of better results.
 
 PyNAST depends on PyCogent, NumPy, Python, and uclust. The versions 
 used for development are:
 
 PyCogent 1.4.1
 NumPy 1.3.0
 Python 2.5.1
 uclust 1.1.579
 
The PyNAST algorithm works as follows:

    (1) Using uclust, identify the closest match to a sequence in a template
     alignment. 
    (2) Pairwise align the candidate sequence and template match identified
     in step 1 (default uses the uclust result, but users can specify an
     alternative pairwise aligner).
    (3) Reintroduce gap pattern from the template sequence.
    (4) Identify insertions which expand the template length. For each 
     'template-expanding' insertion, find the nearest gap character in the
     candidate sequence and remove it. 
    (5) Return the aligned candidate sequence.

"""

class UnalignableSequenceError(Exception):
    pass

def pair_hmm_align_unaligned_seqs(seqs,moltype,params={}):
    """
        This needs to be moved to cogent.align.align
    """
    
    seqs = LoadSeqs(data=seqs,moltype=moltype,aligned=False)
    try:
        s1, s2 = seqs.values()
    except ValueError:
        raise ValueError,\
         "Pairwise aligning of seqs requires exactly two seqs."
    
    try:
        gap_open = params['gap_open']
    except KeyError:
        gap_open = 5
    try:
        gap_extend = params['gap_extend']
    except KeyError:
        gap_extend = 2
    try:
        score_matrix = params['score_matrix']
    except KeyError:
        score_matrix = make_dna_scoring_dict(\
         match=1,transition=-1,transversion=-1)
    
    return global_pairwise(s1,s2,score_matrix,gap_open,gap_extend)

def blast_align_unaligned_seqs(seqs,moltype,params={}):
    """ Pairwise align two seqs using bl2seq
    
        This needs to be moved to the blast application controller.
    
    """
    seqs = dict(LoadSeqs(data=seqs,moltype=moltype,aligned=False).items())
    seq_ids = seqs.keys()
    query_id = seq_ids[0]
    subject_id = seq_ids[1]
    if len(seq_ids) != 2:
        raise ValueError,\
         "Pairwise aligning of seqs with blast requires exactly two seqs."
    
    in_filepath1 = get_tmp_filename(tmp_dir='/tmp/',\
        prefix='bl2seq_input1_',suffix='.fasta')
    in_filepath2 = get_tmp_filename(tmp_dir='/tmp/',\
        prefix='bl2seq_input2_',suffix='.fasta')
    in_filepaths = [in_filepath1,in_filepath2]
    out_filepath = get_tmp_filename(tmp_dir='/tmp/',\
        prefix='bl2seq_output_',suffix='.fasta')
     
    for n,in_filepath in zip(seq_ids,in_filepaths):
        f = open(in_filepath,'w')
        f.write('>%s\n' % n)
        f.write(str(seqs[n]))
        f.write('\n')
        f.close()
    
    # Note: -S 1 indicated that we don't want to blast both orientations -- at
    # this would be different behavior than other pairwise aligners.
    bl2seq_res = system('bl2seq -i %s -j %s -o %s -F F -S 1 -q -1 -p blastn -VT' %\
     (in_filepath1,in_filepath2,out_filepath))
    if bl2seq_res != 0:
        raise RuntimeError, "bl2seq failed:\n %s" % bl2seq_res 
    
    query_seq = []
    subject_seq = []
    blast_res = open(out_filepath)
    in_result = False
    for line in blast_res:
        if line.strip().startswith('Score'):
            if in_result:
                break
            else:
                in_result = True
            
        if line.startswith('Query: '):
            fields = line.split()
            query_seq.append(fields[2].upper())
        elif line.startswith('Sbjct: '):
            fields = line.split()
            subject_seq.append(fields[2].upper())
        else:
            continue
     
    remove(in_filepath1)
    remove(in_filepath2)
    remove(out_filepath)
    
    # reintroduce terminal characters which were not aligned -- this
    # needs to be split out to another function to facilitate easier testing       
    q = ''.join(query_seq)
    q = q.replace('-','')
    s = ''.join(subject_seq)
    s = s.replace('-','')
    query_in = str(seqs[query_id])
    subject_in = str(seqs[subject_id])
    q_start = query_in.index(q[:100])
    q_end = q_start + len(q)
    s_start = subject_in.index(s[:100])
    s_end = s_start + len(s)
    
    five_prime_bases_to_add = max(q_start,s_start)
    three_prime_bases_to_add = max(len(query_in)-q_end, len(subject_in)-s_end)
    
    if five_prime_bases_to_add:
        leading_bases = query_in[:q_start]
        query_seq = '%s%s%s' % \
         ('-'*(five_prime_bases_to_add-len(leading_bases)),\
          leading_bases, 
          ''.join(query_seq))
         
        leading_bases = subject_in[:s_start]
        subject_seq = '%s%s%s' % \
         ('-'*(five_prime_bases_to_add-len(leading_bases)),\
          leading_bases,\
          ''.join(subject_seq))
         
    if three_prime_bases_to_add:
        trailing_bases = query_in[q_end:]
        query_seq = '%s%s%s' %\
         (''.join(query_seq),\
          trailing_bases,\
          '-'*(three_prime_bases_to_add-len(trailing_bases)))
          
        trailing_bases = subject_in[s_end:]
        subject_seq = '%s%s%s' %\
         (''.join(subject_seq),\
          trailing_bases,\
          '-'*(three_prime_bases_to_add-len(trailing_bases)))

    result = [(query_id,query_seq),\
              (subject_id,subject_seq)]
    
    return LoadSeqs(data=result,moltype=moltype)


def align_two_seqs(template, candidate,
    align_unaligned_seqs_f=muscle_align_unaligned_seqs,
    params={},moltype=DNA):
    """ Align the two sequences with an arbitrary aligner function
    
        template: the template sequence to align (string)
        candidate: the candidate sequence to align (string)
        align_unaligned_seqs_f: function to be applied to aligned the
         candidate and template sequences -- function must be of the form
         align_unaligned_seqs_f(seqs,moltype,params=params)
        params: params to be passed to align_unaligned_seqs
        moltype: moltype to be passed to align_unaligned_seqs
    """
    # Load the sequences into a form useful to align_unaligned_seq_f
    seqs = [('template',str(template)), ('candidate',str(candidate))]
    # Align the sequences
    aln = align_unaligned_seqs_f(seqs,moltype,params=params)
    # Extract the sequences from the alignment object and return them
    return aln.getGappedSeq('template'), aln.getGappedSeq('candidate')
 
def reintroduce_template_spacing(template,
    pw_aligned_template,pw_aligned_candidate):
    """ reintroduce template gap spacing into pairwise aligned sequences
    """
    # Check for the simple case where the alignment reproduced the
    # template spacing
    if template == pw_aligned_template:
        return (pw_aligned_template, pw_aligned_candidate,[])   
    
    # get gap maps to help with relating the aligned template sequence
    # to the pairwise aligned template and candidate sequences
    template_seq_to_aln = template.gapMaps()[0]
    pw_template_seq_to_aln, pw_template_aln_to_seq = \
     pw_aligned_template.gapMaps()
     
    # build a list to keep track of gaps that were introduced in
    # the pairwise alignment but which were not present in the template
    # alignment 
    new_gaps_in_pw_alignment = [] 
    # create variable to keep track of how many gaps have been 
    # reintroduced so far from the template to the pw_aligned_template -
    # this is necessary to efficently compute new_gaps_in_pw_alignment
    total_reintroduced_gaps = 0

    template_result = list(pw_aligned_template)
    candidate_result = list(pw_aligned_candidate) 
    # begin iteration over the alignment positions
    for aln_curr_pos in range(len(pw_aligned_template)):
        try:
            # map the current alignment position to the 
            # corresponding sequence (ie. ungapped) position
            seq_curr_pos = \
             pw_template_aln_to_seq[aln_curr_pos]
        except KeyError:
            # if the current alignment position is a gap, move
            # on to the next alignment position
            continue
        # store the next sequence position as it is used in several places
        seq_next_pos = seq_curr_pos + 1
        try:
            # Get the number of gaps between the next and current 
            # alignment positions in the template alignment
            template_post_char_gaps = \
             template_seq_to_aln[seq_next_pos] - \
             template_seq_to_aln[seq_curr_pos] - 1
        except KeyError:
            # at the end of the sequence
            break
        # Get the number of gaps between the next and current 
        # alignment positions in the template sequence in the
        # pairwise alignment
        pw_template_post_char_gaps = \
         pw_template_seq_to_aln[seq_next_pos] -\
         aln_curr_pos - 1

        # compute the difference in the number of gaps following the
        # current position in the two alignments
        addl_gaps = template_post_char_gaps - pw_template_post_char_gaps
        if addl_gaps > 0:
            # if the additional gaps is greater than zero, additional
            # gap characters need to be added to the pairwise alignment
            insertion_point = aln_curr_pos + 1 + total_reintroduced_gaps 
            template_result[insertion_point:insertion_point] = ['-'] * addl_gaps
            candidate_result[insertion_point:insertion_point] = ['-'] * addl_gaps
            # update the tally of reintroduced gaps
            total_reintroduced_gaps += addl_gaps
        elif addl_gaps < 0:
            # if the additional gaps is less than zero, the pairwise 
            # alignment introduced new gaps -- store these positions to be 
            # dealt with later. Note that first_new_gap_pos is
            # adjusted by adding the number of the gap characters
            # reintroduced to the current point. Positions 
            # in new_gaps_in_pw_alignment therefore refer to positions in 
            # the alignments being returned from this function
            first_new_gap_pos = aln_curr_pos + total_reintroduced_gaps + 1
            # add the positions of the new gaps chars to the list
            # of new gaps
            new_gaps_in_pw_alignment += \
             range(first_new_gap_pos,first_new_gap_pos + (-1*addl_gaps))
        else:
            # gap pattern is the same following the current sequence 
            # position
            pass

    return (DNA.makeSequence(''.join(template_result)), \
            DNA.makeSequence(''.join(candidate_result)),\
            new_gaps_in_pw_alignment)           
            
def nearest_gap(seq,pos):
    """ Returns the position of the nearest gap to pos in seq
    """
    # Catch negative sequence positions
    if pos < 0:
        raise IndexError, "Sequence positions cannot be negative: %d" % pos
    
    # If pos contains a gap, that's the closest gap
    if seq[pos] == '-':
        return pos
        
    # create a list to store the nearest gap character in the 5' and
    # 3' directions
    choices = []
    # find the nearest gap 5' of pos
    try:
        gap_index = ''.join(seq[:pos]).rindex('-')
        distance = pos - gap_index
        choices.append((distance,gap_index))
    except ValueError:
        pass
        
    # find the nearest gap 3' of pos
    try:
        gap_index = pos + ''.join(seq[pos:]).index('-')
        distance = gap_index - pos
        choices.append((distance,gap_index))
    except ValueError:
        pass
    
    # error if there are no gaps in the sequence
    if not choices:
        raise UnalignableSequenceError,\
         "Can't adjust alignment because there are too few gaps to "+\
         "remove in the aligned candidate to reduce to the length of "+\
         "the template alignment (i.e., candidate adds too many insertions "+\
         "during pairwise alignment)."
        
    # return the gap_index of the choice with the smaller distance -- if there
    # is a tie, will delete the 5' gap (which is what original NAST does)
    return min(choices)[1]
            
def adjust_alignment(template,candidate,new_gaps):
    """adjust template/candidate aln to remove gaps added by pairwise alignment
    
        This step adjusts the alignment to reduce the length back to the 
         template alignment length by introducing local misalignments to
         remove gap characters that are present in the pairwise alignment
         but not in the template alignment.
    
    """
    template_l = list(template)
    candidate_l = list(candidate)
    new_gaps.reverse()
    for pos in new_gaps:
        del template_l[pos]
        del candidate_l[nearest_gap(candidate_l,pos)]
        
    return (DNA.makeSequence(''.join(template_l)), \
            DNA.makeSequence(''.join(candidate_l)))
        
def introduce_terminal_gaps(template,aligned_template,aligned_candidate):
    """ introduce terminal gaps from template into the aligned candidate seq
    """
    
    # count the 5' gaps in the original aligned template
    original_five_prime_gaps = 0
    for c in template:
        if c == '-':
            original_five_prime_gaps +=1
        else:
            break
            
    # count the 5' gaps already existing in the pairwise aligned template
    # (because we don't need to add these)
    aligned_template_five_prime_gaps = 0
    for c in aligned_template:
        if c == '-':
            aligned_template_five_prime_gaps += 1
        else:
            break
            
    # compute the number of 5' gaps that need to be added to get to the
    # original alignment length
    five_prime_gaps_to_add = \
     original_five_prime_gaps - aligned_template_five_prime_gaps
            
    # count the 3' gaps in the original aligned template
    original_three_prime_gaps = 0
    for c in reversed(template):
        if c == '-':
            original_three_prime_gaps +=1
        else:
            break
            
    # count the 3' gaps already existing in the pairwise aligned template
    # (because we don't need to add these)
    aligned_template_three_prime_gaps = 0
    for c in reversed(aligned_template):
        if c == '-':
            aligned_template_three_prime_gaps += 1
        else:
            break
            
    # compute the number of 3' gaps that need to be added to get to the
    # original alignment length
    three_prime_gaps_to_add = \
     original_three_prime_gaps - aligned_template_three_prime_gaps

    # return the sequence with the 5' and 3' gaps added
    return DNA.makeSequence(''.join([\
     '-'*five_prime_gaps_to_add,\
     str(aligned_candidate),\
     '-'*three_prime_gaps_to_add]),\
     Name=aligned_candidate.Name)
    
def remove_template_terminal_gaps(candidate,template):
    """Remove template terminal gaps and corresponding bases in candidate 
    """
    if len(template) != len(candidate):
        raise ValueError, \
         "Sequences must be aligned, but their "+\
         "lengths aren't equal. %d != %d" % (len(candidate),len(template))
         
    if len(template) == 0:
        return candidate, template
    
    degapped_candidate_len = len(candidate.degap())
    
    candidate = DNA.makeSequence(candidate)
    template = DNA.makeSequence(template)
    
    template_gap_vector = template.gapVector()
    first_non_gap = template_gap_vector.index(False)
    num_three_prime_gaps = template_gap_vector[::-1].index(False)
    last_non_gap = len(template_gap_vector) - num_three_prime_gaps
    
    # Construct the candidate name, which will include the range of bases
    # from the original sequence
    candidate = candidate[first_non_gap:last_non_gap]
    template = template[first_non_gap:last_non_gap]
    candidate_start_pos = first_non_gap + 1
    candidate_end_pos = degapped_candidate_len - num_three_prime_gaps
    candidate_name = candidate.Name
    if candidate_name.endswith('RC'):
        name_delimiter = ':'
    else:
        name_delimiter = ' '
    candidate_name = '%s%s%d..%d' %\
     (candidate_name,name_delimiter,candidate_start_pos,candidate_end_pos)
    
    return DNA.makeSequence(candidate,Name=candidate_name), template

def depreciation_warning(d):
    if d:
        print "Unsupported or depreciated options "+\
         "passed to pynast: %s\n" % ' '.join(d.keys()) +\
         "  blast_db, max_e_value, and addl_blast_params are depreciated " +\
         "and will be removed in PyNAST 1.2."

def pynast_seq(candidate_sequence, template_alignment,
    max_hits=30, min_pct=75.0, min_len=1000, align_unaligned_seqs_f=None,
    **kwargs):
    """ Apply PyNAST to a single sequence 
    
    candidate_sequence
        a single DNA sequence object
    template_alignment
        a PyCogent alignment object containing the template alignment
        or a fasta filepath
    max_hits
      Maximum number of uclust hits to return
    min_pct
      minimum % identity for best database match
    min_len
      minimum length of match for alignment     
    align_unaligned_seqs_f
      Function to align sequences. Must be of the form:
       align_unaligned_seqs(seqs, moltype, params=None)
       see cogent.app.muscle.align_unaligned_seqs
    """
    depreciation_warning(kwargs)
    class SingleSeqLogger(object):
        """ A simple object to store results of a single pynast run """
        
        def setUp(self):
            self.Data = None
        
        def record(self,*args):
            self.Data = tuple(args)
    
    l = SingleSeqLogger()
    candidate_sequences = [(candidate_sequence.Name,str(candidate_sequence))]
    
    aligned_seq, exit_status = list(ipynast_seqs(candidate_sequences,
     template_alignment, max_hits=max_hits, min_pct=min_pct, min_len=min_len,
     align_unaligned_seqs_f=align_unaligned_seqs_f,
     log_fp=None, logger=l))[0]
    
    if exit_status == 0:
        return l.Data[3], aligned_seq
    else:
        raise UnalignableSequenceError, l.Data[2]

def ipynast_seqs(candidate_sequences, template_alignment,
    max_hits=30, min_pct=75.0, min_len=1000, align_unaligned_seqs_f=None,
    log_fp=None, logger=None,**kwargs):
    """Iterator that yields results of pynast on candidate_sequences
    
    This function yields the sequence and exit status of the alignment step,
     as (sequence, exit status) tuples.
     Status values can be:
       0 : indicates a sucessful alignment, in which case the sequence will be
            aligned
       1 : indicates unsucessful sequence search, in which case the sequence 
            will be unaligned
       2 : indicates alignment did not meet minimum requirements, in which case 
            the sequence will be unaligned
            
     All sequences are returned as DNA sequence objects.
    
    candidate_sequences
        an iterable object (e.g., a list) containing tuples of
        (seq_id, sequence) pairs (e.g., as returned by MinimalFastaParser)
        or a fasta filepath
    template_alignment
        a PyCogent alignment object containing the template alignment
        or a fasta filepath
    max_hits
      Maximum number of uclust hits to return
    min_pct
      minimum % identity for best database match
    min_len
      minimum length of match for alignment     
    align_unaligned_seqs_f
      Function to align sequences. Must be of the form:
       align_unaligned_seqs(seqs, moltype, params=None)
       see cogent.app.muscle.align_unaligned_seqs
    log_fp
      Optional path to log file
    logger
      Optional NastLogger object, takes precedence over log_fp
      
    """
    depreciation_warning(kwargs)
    
    files_to_remove = []
    if type(candidate_sequences) == str:
        # filepath provided for candidate sequences
        candidate_fasta_filepath = candidate_sequences
    else:
        # sequence list provided for candidate sequence -- write 
        # the seqs to a temp file to pass to uclust
        candidate_fasta_filepath = \
         get_tmp_filename(prefix='pynast_candidate',suffix='.fasta')
        candidate_fasta_f = open(candidate_fasta_filepath,'w')
        for s in candidate_sequences:
            candidate_fasta_f.write('>%s\n%s\n' % (s))
        candidate_fasta_f.close()
        files_to_remove.append(candidate_fasta_filepath)

    # degap the template alignment for the sequence searching step
    if type(template_alignment) == str:
        # template alignment provided as filepath -- load it into
        # an Alignment object
        template_alignment = LoadSeqs(template_alignment,moltype=DNA)
    # degap the alignment, and write it to a temp file
    template_fasta_filepath = \
     get_tmp_filename(prefix='pynast_template',suffix='.fasta')
    template_fasta_f = open(template_fasta_filepath,'w')
    template_fasta_f.write(template_alignment.degap().toFasta())
    template_fasta_f.close()
    files_to_remove.append(template_fasta_filepath)
         
    # Set up logging.  NastLogger object takes precedence over log
    # file path, if both are provided.
    if logger is not None:
        logger = logger
    elif log_fp is not None:
        logger = NastLogger(log_fp)
    else:
        logger = NastLogger()
    
    min_pct /= 100.
    # get the alignment iterator
    pw_alignment_iterator = uclust_search_and_align_from_fasta_filepath(
            candidate_fasta_filepath,
            template_fasta_filepath,
            percent_ID=min_pct,
            enable_rev_strand_matching=True)

    try:
        current_result = pw_alignment_iterator.next()
    except StopIteration:
        current_result = None
        
    for seq_id, seq in MinimalFastaParser(open(candidate_fasta_filepath)):
        seq_len = len(seq)
        if '-' in seq:
            # clean-up temporary blast database files if any were created
            pw_alignment_iterator.close()
            remove_files(files_to_remove,error_on_missing=False)
            raise ValueError, "Candidate sequence contains gaps. This is not supported."
        
        try:
            candidate_seq_id, template_seq_id, pw_aligned_candidate,\
             pw_aligned_template, pct_identity = current_result
        except TypeError:
            pass
        
        if not current_result or seq_id.split()[0] != candidate_seq_id.split()[0]:
            # a suitable match was not found - don't align the sequence
            # log the failure
            logger.record(
                seq_id, # input sequence identifier
                len(seq), # input sequence length
                "No search results.")
            # yield the unaligned sequence and failure code
            yield DNA.makeSequence(seq,Name=seq_id), 1
        else:
            # this sequence was aligned
            if align_unaligned_seqs_f:
                # if an alternate pairwise aligner was specified, unalign
                # and re-align the sequences.
                pw_aligned_template, pw_aligned_candidate =\
                 align_two_seqs(pw_aligned_template.replace('-',''),
                                pw_aligned_candidate.replace('-',''),
                                align_unaligned_seqs_f)
                                    
            # Cast the pairwise alignments to DNA sequence objects
            pw_aligned_candidate = \
             DNA.makeSequence(pw_aligned_candidate,Name=candidate_seq_id)
            pw_aligned_template = \
             DNA.makeSequence(pw_aligned_template,Name=template_seq_id)
    
            # Remove any terminal gaps that were introduced into the template
            # sequence
            pw_aligned_candidate, pw_aligned_template = \
                remove_template_terminal_gaps(
                pw_aligned_candidate, pw_aligned_template)
            candidate_seq_id = pw_aligned_candidate.Name
    
            # get the aligned template sequence from the template alignment
            template_aligned_seq = \
             template_alignment.getGappedSeq(template_seq_id)
    
            # reintroduce the gap spacing from the template alignment
            pw_aligned_template, pw_aligned_candidate, new_gaps =\
              reintroduce_template_spacing(template_aligned_seq,\
              pw_aligned_template,pw_aligned_candidate)
    
            # delete any new gaps that were introduced during the 
            # pairwise alignment step
            pw_aligned_template, pw_aligned_candidate = adjust_alignment(\
             pw_aligned_template,pw_aligned_candidate,new_gaps)
     
            # reintroduce any terminal gaps that were present in the template
            result = introduce_terminal_gaps(\
                template_aligned_seq,pw_aligned_template,pw_aligned_candidate)
        
            unaligned_length = len(result.degap())
            if unaligned_length < min_len:
                # alignment is too short - log this as a failure
                error = "Alignment does not meet minimum length "+\
                            "requirement for alignment (%d < %d)"\
                             % (seq_len,min_len)
                logger.record(
                    seq_id, # input sequence identifier
                    len(seq), # input sequence length
                    "No search results.")
                # yield the unaligned sequence and failure code
                yield DNA.makeSequence(seq,Name=seq_id), 2
            else:        
                # log the alignment
                logger.record(
                    seq_id, # input sequence identifier
                    len(seq), # input sequence length
                    '',                  # Errors
                    template_seq_id, # best template match id
                    '%3.2f' % pct_identity, # pct id to template
                    unaligned_length, # post alignment sequence length
                    )

                # yield the aligned sequence and sucess code
                yield DNA.makeSequence(result,Name=candidate_seq_id), 0
                
            # get the next alignment
            try:
                current_result = pw_alignment_iterator.next()
            except StopIteration:
                # end of the input fasta file indicates completion,
                # not end of the aligned sequences
                continue

    # clean-up temporary blast database files if any were created
    remove_files(files_to_remove,error_on_missing=False)


def null_status_callback_f(x):
    """Dummy function to pass as default status_callback_f"""
    pass

def pynast_seqs(candidate_sequences, template_alignment, max_hits=30,
    min_pct=75.0, min_len=1000, align_unaligned_seqs_f=None, log_fp=None,
    logger=None, status_callback_f=null_status_callback_f,**kwargs):
    """Function which runs pynast_seq on candidate_sequences.
    
    Results are returned as a tuple of lists:
     (aligned_sequences, failed_to_align_sequences)
     where all sequences are DNA sequence objects.
   
    candidate_sequences
        an iterable object (e.g., a list) containing tuples of
        (seq_id, sequence) pairs (e.g., as returned by MinimalFastaParser)
        or a fasta filepath
    template_alignment
        a PyCogent alignment object containing the template alignment
        or a fasta filepath
    max_hits
      Maximum number of uclust hits to return
    min_pct
      minimum % identity for best database match
    min_len
      minimum length of match for alignment     
    align_unaligned_seqs_f
      Function to align sequences. Must be of the form:
       align_unaligned_seqs(seqs, moltype, params=None)
       see cogent.app.muscle.align_unaligned_seqs
    log_fp
      Optional path to log file
    logger
      Optional NastLogger object, takes precedence over log_fp
    status_callback_f:
      Callback function to provide status updates to callers of pynast_seqs.
      This function must take a single parameter.
    """
    depreciation_warning(kwargs)
    # create lists to keep track of the aligned candidate sequences 
    # and the sequences which fail to align
    aligned = []
    failed_to_align = []
    
    pynast_iterator = ipynast_seqs(
     candidate_sequences, template_alignment,
     max_hits=max_hits, min_pct=min_pct, min_len=min_len,
     align_unaligned_seqs_f=align_unaligned_seqs_f, log_fp=log_fp,
     logger=logger)
    
    for seq, status in pynast_iterator:
        if status == 0:
            aligned.append(seq)
            status_callback_f(seq)
        else:
            failed_to_align.append(seq)
            status_callback_f(seq)
    
    return aligned, failed_to_align

pairwise_alignment_methods = {\
     'muscle':muscle_align_unaligned_seqs,\
     'mafft':mafft_align_unaligned_seqs,\
     'clustal':clustal_align_unaligned_seqs,\
     'blast':blast_align_unaligned_seqs,\
     'pair_hmm':pair_hmm_align_unaligned_seqs,\
     'uclust':None}
