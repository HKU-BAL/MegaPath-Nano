import os
from pyssw import SSW as ssw

align_param = {
"match": 4,
"mismatch": 6,
"gap_open": 8,
"gap_extend": 2,
"k": 23,
"error_rate": 0.009999999776482582,
"read_size": 250,
"kmer_size": 32,
"max_num_of_mismatches": 2,
"realignment_similarity_threshold": 0.16934,
}
transform = ['?', 'X', "I", "D", 'Skip', 'S', 'SH']
lib_path = os.path.join(*(os.path.split(__file__)[:-1] + ('realign/libssw.so',)))

class FastPassAligner(object):
  def __init__(self, ctg_name=None, consensus=None, reference=None, reference_start=0, reference_end=0, reference_prefix=0, reference_suffix=0, read_name_list=None):
    self.reference = reference
    self.reference_start = reference_start
    self.reference_end = reference_end
    self.align_param = align_param
    self.ctg_name = ctg_name
    self.consensus = consensus
    self.read_size = 250
    self.kmer_size = align_param['kmer_size']
    self.kmer_index_ = {}
    self.read_name_list = read_name_list
    self.match = align_param['match']
    self.mismatch = align_param['mismatch']
    self.gap_open = align_param['gap_open']
    self.gap_extend = align_param['gap_extend']
    self.aligner = ssw(match=self.match, mismatch=self.mismatch, gap_open=self.gap_open, gap_extend=self.gap_extend, lib_path=lib_path)

  def AlignHaplotypesToReference(self):
    self.consensus_info = [[]] * len(self.consensus)
    for cons_id, (read_name, consensus) in enumerate(zip(self.read_name_list,self.consensus)):
      sw_score, cigar, ref_begin = self.aligner.align(consensus)
      if sw_score > 0:
        ref_pos = ref_begin + self.reference_start

      self.consensus_info[cons_id] = [consensus, cigar, ref_pos, read_name]

  def align_reads(self):
    self.aligner.set_reference_sequence(self.reference)
    self.AlignHaplotypesToReference()
    return self.consensus_info



