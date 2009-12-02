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
from cogent.parse.blast import BlastResult
from cogent.parse.fasta import MinimalFastaParser
from pynast.logger import NastLogger

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, the PyNAST project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Beta"

""" PyNAST is a complete rewrite of the NAST algorithm written in python. 
 The dependencies are PyCogent, NumPy, Python, BLAST, and muscle. The versions 
 used for development are:
 
 PyCogent 1.4.0.dev
 NumPy 1.3.0
 Python 2.5.1
 blastall v2.2.20
 muscle v3.6
 
The NAST algorithm, reimplemented here, works as follows:

    (1) Using BLAST, identify the closest match to a sequence in a template
     alignment. 
    (2) Trim the candidate sequence to the beginning and end points of the 
     best BLAST hit. 
    (3) Perform a pairwise alignment between the candidate sequence and
     template alignment.
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
     

def blast_sequence(candidate_sequence,blast_db,max_hits,\
    max_e_value=1e-1,addl_blast_params={}):
 
    # set default blast parameters
    params = {
        # max procs
        "-a":"1",
        # do not filter query sequences
        "-F":"F",
        # mismatch penalty
        "-q":"-1"}
    # override blast params with any additional passed params
    params.update(addl_blast_params)
    if not(exists(blast_db + ".nsq")):
        raise IOError, \
         "BLAST database: '%s' probably does not exist." % blast_db
    # blast sequence
    
    result = blastn(\
        ['>query_%s' % candidate_sequence.Name,str(candidate_sequence)],\
        blast_db=blast_db,
        e_value=str(max_e_value), 
        max_hits=max_hits, working_dir="/tmp/", blast_mat_root=None, 
        extra_params=params)
    return result
    
def process_blast_result(blast_res, templ_align, templ_align_seq, cur_seq, 
    min_len, min_pct, num_hits=30):
    """ Post process blast results 
     Originally derived from nast_farm.py
    """
    # if no hits, bail out
    if not (blast_res and blast_res.keys() != ['']):
        raise UnalignableSequenceError, "No blast results."
    # if hits for multiple queries, bail out
    if len(blast_res) > 1: 
        raise ValueError,\
         "Blast results must processed one query at a time. "+\
         "Results provided for multiple queries: %s " \
         % ' '.join(blast_res.keys())
        
    ## variables to track best hit
    # will store the best blast hit list
    best_hit = None
    # will store the best % idnetity score
    best_pct_identity = None
    # will store the id of the best blast hit sequence
    best_templ_id = None
    # will store the number of aligned query positions for the
    # best blast hit
    best_aligned_positions = None
    # will store the alignment length of the best hit
    best_aln_len = None
    
    # find the best hit that meets threshold values 
    for query_id, hits in blast_res.bestHitsByQuery(n=num_hits):
        
        for hit in hits:
            # check current hit
            cur_pct_identity = float(hit[BlastResult.PERCENT_IDENTITY])
            cur_aln_len = int(hit[BlastResult.ALIGNMENT_LENGTH])
            # count the number of positions which are aligned --
            # this is used to identify the best blast hit in original
            # NAST, may want to parameterize how that is determined
            cur_aligned_positions = \
             int(hit[BlastResult.QUERY_END]) - \
             int(hit[BlastResult.QUERY_START]) + 1
            # check if current hit is the best blast hit
            if (best_aligned_positions is None) or \
                    (cur_aligned_positions > best_aligned_positions) or\
                    (cur_aligned_positions == best_aligned_positions \
                         and cur_pct_identity > best_pct_identity):
                best_pct_identity = cur_pct_identity 
                best_aln_len = cur_aln_len
                best_templ_id = hit[BlastResult.SUBJECT_ID] 
                best_hit = hit
                best_aligned_positions = cur_aligned_positions

    # raise error if best hit does not meet user specified criteria 
    if best_pct_identity < min_pct or best_aln_len < min_len:
        raise UnalignableSequenceError,\
         ''.join(["Best alignment did not meet user requirements. ",
            "Length: %d " % best_aln_len,
            "Percent Indentity: %.2f " % best_pct_identity,
            "Template: %s" % best_templ_id])
  
    # get ids
    query_id = best_hit[BlastResult.QUERY_ID]
    subj_id = best_hit[BlastResult.SUBJECT_ID]
    q_start = int(best_hit[BlastResult.QUERY_START])
    q_end = int(best_hit[BlastResult.QUERY_END])
    s_start = int(best_hit[BlastResult.SUBJECT_START])
    s_end = int(best_hit[BlastResult.SUBJECT_END])
    
    # check orientation of blast results
    query_rev = False
    if  s_start > s_end:
        templ_seq = templ_align_seq.getSeq(subj_id)
        len_cur_seq = len(cur_seq)
        cur_seq = DNA.makeSequence(cur_seq[q_start-1:q_end].rc(),\
         Name='%s RC:%s..%s' %\
          (cur_seq.Name,str(len_cur_seq-q_end+1),str(len_cur_seq-q_start+1)))
    elif q_start > q_end:
        # is it possible to get here???
        raise NotImplementedError,\
         "Seqs match is to a rc of database seq."
    else:
        # no orientation adjustments are necessary
        templ_seq = templ_align_seq.getSeq(subj_id)
        cur_seq = DNA.makeSequence(cur_seq[q_start-1:q_end],\
            Name='%s %s..%s' % (cur_seq.Name,str(q_start),str(q_end)))

    pct_identity = best_hit[BlastResult.PERCENT_IDENTITY]

    return templ_seq, cur_seq, pct_identity

def align_two_seqs(template, candidate,\
    align_unaligned_seqs_f=muscle_align_unaligned_seqs,\
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
 
def reintroduce_template_spacing(template,\
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
    
def pynast_seq(candidate_sequence,template_alignment,\
    degapped_template_alignment=None,blast_db=None,max_hits=30,\
    max_e_value=1e-1,addl_blast_params={},min_pct=75.0,min_len=1000,\
    align_unaligned_seqs_f=blast_align_unaligned_seqs, log_fp=None, logger=None):
    """ align candidate sequence to template aln, preserving the aln length
    """
    if candidate_sequence.isGapped():
        raise ValueError,\
         "Candidate sequence cannot contain gap characters: %s"\
         % candidate_sequence.Name
         
    # Set up logging.  NastLogger object takes precedence over log
    # file path, if both are provided.
    if logger is not None:
        logger = logger
    elif log_fp is not None:
        logger = NastLogger(log_fp)
    else:
        logger = NastLogger()
        
    degapped_template_alignment = degapped_template_alignment or\
     template_alignment.degap()

    # if a blast_db wasn't passed, create one from the (degapped) 
    # template alignment
    if not blast_db:
        blast_db, db_files_to_remove = \
         build_blast_db_from_seqs(template_alignment,output_dir='/tmp/')
    else:
        db_files_to_remove = []
        
    # blast the candidate sequence against blast_db
    blast_result = \
     blast_sequence(candidate_sequence,blast_db,max_hits,\
     max_e_value,addl_blast_params)
    
    # get the unaligned template sequence (of the best blast hit) and 
    # candidate sequence trimmed to only contain the bases which aligned to
    # the template sequence
    template_seq, candidate_seq, pct_identity = \
     process_blast_result(blast_result,template_alignment,\
      degapped_template_alignment, candidate_sequence,\
      min_len=min_len, min_pct=min_pct, num_hits=max_hits)
    candidate_seq_id = candidate_seq.Name
    
    # store the id of the best blast hit
    template_seq_id = template_seq.Name
    
    # get the aligned template sequence from the template alignment
    template_aligned_seq = template_alignment.getGappedSeq(template_seq_id)
        
    # aligned template_seq and candidate_seq
    pw_aligned_template, pw_aligned_candidate = \
     align_two_seqs(template_seq, candidate_seq, align_unaligned_seqs_f)
    
    # reintroduce the gap spacing from the template alignment
    pw_aligned_template, pw_aligned_candidate, new_gaps =\
      reintroduce_template_spacing(template_aligned_seq,\
      pw_aligned_template,pw_aligned_candidate)
    
    # delete any new gaps that were introduced during the pairwise alignment
    # step
    pw_aligned_template, pw_aligned_candidate = adjust_alignment(\
     pw_aligned_template,pw_aligned_candidate,new_gaps)
     
    # reintroduce any terminal gaps that were present in the template
    result = introduce_terminal_gaps(\
        template_aligned_seq,pw_aligned_template,pw_aligned_candidate)
     
    # clean-up temporary blast database files if any were created
    remove_files(db_files_to_remove,error_on_missing=False)

    # log the alignment
    # TODO: Fill in or remove missing fields
    logger.record(
        candidate_seq.Name,
        len(candidate_seq),
        '',                  # Errors
        template_seq_id,
        pct_identity,
        len(result.degap()), # Candidate Sequence length post-Nast
        )

    # return the id of the best blast hit and the aligned candidate sequence
    result = DNA.makeSequence(result,Name=candidate_seq_id)
    return template_seq_id, result

def ipynast_seqs(candidate_sequences,template_alignment,\
    blast_db=None,max_hits=30,\
    max_e_value=1e-1,addl_blast_params={},min_pct=75.0,min_len=1000,\
    align_unaligned_seqs_f=blast_align_unaligned_seqs, log_fp=None,\
    logger=None):
    """Iterator that yields results of pynast_seq on candidate_sequences
    
    This function yields the sequence and exit status of the alignment step,
     as (sequence, exit status) tuples.
     Status values can be:
       0 : indicates a sucessful alignment, in which case the sequence will be
            aligned
       1 : indicates unsucessful alignment (i.e., an UnalignableSequenceError
            was caught), in which case the sequence will be unaligned
            
     Statuses are returned as ints (rather than bools) to support new statuses
      in the future. 0 is used to indicate success because there is only one
      way the call can result in a sucessful alignment, but we may want to
      allow for alternate failure types in the future.
     All sequences are returned as DNA sequence objects.
    
    candidate_sequences
        an iterable object (e.g., a list) containing tuples of
        (seq_id, sequence) pairs (e.g., as returned by MinimalFastaParser)
    template_alignment
        a PyCogent alignment object containing the template alignment
    blast_db
      Database to BLAST against (default: derived from template)
    max_hits
      Maximum number of BLAST hits to return (passed to 
      cogent.app.blast.blastn())
    max_e_value
      expectation value passed to BLAST application controller
    addl_blast_params
      additional parameters to pass to the BLAST application controller 
      (see documentation for cogent.app.blast.blastn()) 
    min_pct
      minimum % identity for BLAST hit
    min_len
      minimum length of match for BLAST hit      
    align_unaligned_seqs_f
      Function to align sequences. Must be of the form:
       align_unaligned_seqs(seqs, moltype, params=None)
       see cogent.app.muscle.align_unaligned_seqs
    log_fp
      Optional path to log file
    logger
      Optional NastLogger object, takes precedence over log_fp
    """
    
    # Set up logging.  NastLogger object takes precedence over log
    # file path, if both are provided.
    if logger is not None:
        logger = logger
    elif log_fp is not None:
        logger = NastLogger(log_fp)
    else:
        logger = NastLogger()
    
    # if a blast database was not passed in, create one from the
    # (degapped) template alignment
    if not blast_db:
        blast_db, db_files_to_remove = \
         build_blast_db_from_seqs(template_alignment,output_dir='/tmp/')
    else:
        db_files_to_remove = []
        
    # cache a copy of the degapped template alignment
    degapped_template_alignment = template_alignment.degap()
    
    # iterate over the candidate sequences
    for seq_id, seq in candidate_sequences:
        candidate_sequence = DNA.makeSequence(seq,Name=seq_id)
        try:
            # align the candidate sequence to the template alignment
            # and store the best blast hit id and aligned sequence
            template_seq_id, aligned_seq = \
             pynast_seq(candidate_sequence,template_alignment,\
             degapped_template_alignment=degapped_template_alignment,\
             blast_db=blast_db,max_hits=max_hits,max_e_value=max_e_value,\
             addl_blast_params=addl_blast_params,min_pct=min_pct,\
             min_len=min_len,align_unaligned_seqs_f=align_unaligned_seqs_f,\
             logger=logger)
            yield aligned_seq, 0
            
        except UnalignableSequenceError,e:
            # if the sequence could not be aligned, store the sequence
            # in the failures list
            yield DNA.makeSequence(seq,Name=seq_id), 1
            logger.record(seq_id, len(seq), e)
            
        except ValueError,e:
            # clean-up temporary blast database files if any were created
            remove_files(db_files_to_remove,error_on_missing=False)
            # Re-raise the error
            e = str(e) + "\nseq_id:%s\nseq:%s\n" % (seq_id, seq)
            raise ValueError, e

    # clean-up temporary blast database files if any were created
    remove_files(db_files_to_remove,error_on_missing=False)

def null_status_callback_f(x):
    """Dummy function to pass as default status_callback_f"""
    pass

def pynast_seqs(candidate_sequences,template_alignment,blast_db=None,max_hits=30,\
    max_e_value=1e-1,addl_blast_params={},min_pct=75.0,min_len=1000,\
    align_unaligned_seqs_f=blast_align_unaligned_seqs, log_fp=None,\
    logger=None, status_callback_f=null_status_callback_f):
    """Function which runs pynast_seq on candidate_sequences.
    
    Results are returned as a tuple of lists:
     (aligned_sequences, failed_to_align_sequences)
     where all sequences are DNA sequence objects.
    
    candidate_sequences
        an iterable object (e.g., a list) containing tuples of
        (seq_id, sequence) pairs (e.g., as returned by MinimalFastaParser)
    template_alignment
        a PyCogent alignment object containing the template alignment
    blast_db
      Database to BLAST against (default: derived from template)
    max_hits
      Maximum number of BLAST hits to return (passed to 
      cogent.app.blast.blastn())
    max_e_value
      expectation value passed to BLAST application controller
    addl_blast_params
      additional parameters to pass to the BLAST application controller 
      (see documentation for cogent.app.blast.blastn()) 
    min_pct
      minimum % identity for BLAST hit
    min_len
      minimum length of match for BLAST hit      
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

    # create lists to keep track of the aligned candidate sequences 
    # and the sequences which fail to align
    aligned = []
    failed_to_align = []
    
    pynast_iterator = ipynast_seqs(\
     candidate_sequences,template_alignment,blast_db=blast_db,\
     max_hits=max_hits,max_e_value=max_e_value,\
     addl_blast_params=addl_blast_params,min_pct=min_pct,min_len=min_len,\
     align_unaligned_seqs_f=align_unaligned_seqs_f, log_fp=log_fp,\
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
     'pair_hmm':pair_hmm_align_unaligned_seqs}
