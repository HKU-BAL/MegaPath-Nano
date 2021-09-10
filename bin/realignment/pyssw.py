#!/usr/bin/env python
import ctypes as ct
import os
lib_path = os.path.join(*(os.path.split(__file__)[:-1] + ('realign/libssw.so',)))


class CAlignRes(ct.Structure):
    _fields_ = [('nScore', ct.c_uint16),
                ('nScore2', ct.c_uint16),
                ('nRefBeg', ct.c_int32),
                ('nRefEnd', ct.c_int32),
                ('nQryBeg', ct.c_int32),
                ('nQryEnd', ct.c_int32),
                ('nRefEnd2', ct.c_int32),
                ('sCigar', ct.POINTER(ct.c_uint32)),
                ('nCigarLen', ct.c_int32)]

class CProfile(ct.Structure):

    _fields_ = [('pByte', ct.POINTER(ct.c_int32)),
                ('pWord', ct.POINTER(ct.c_int32)),
                ('pRead', ct.POINTER(ct.c_int8)),
                ('pMat', ct.POINTER(ct.c_int8)),
                ('nReadLen', ct.c_int32),
                ('nN', ct.c_int32),
                ('nBias', ct.c_uint8)]



class CSsw(object):

    def __init__(self, sLibPath):
        self.ssw = ct.cdll.LoadLibrary(sLibPath)

        self.ssw_init = self.ssw.ssw_init
        self.ssw_init.argtypes = [ct.POINTER(ct.c_int8), ct.c_int32, ct.POINTER(ct.c_int8), ct.c_int32, ct.c_int8]
        self.ssw_init.restype = ct.POINTER(CProfile)

        self.init_destroy = self.ssw.init_destroy
        self.init_destroy.argtypes = [ct.POINTER(CProfile)]
        self.init_destroy.restype = None
        self.ssw_align = self.ssw.ssw_align
        self.ssw_align.argtypes = [ct.c_void_p, ct.POINTER(ct.c_int8), ct.c_int32, ct.c_uint8, ct.c_uint8, ct.c_uint8, ct.c_uint16, ct.c_int32, ct.c_int32]
        self.ssw_align.restype = ct.POINTER(CAlignRes)

        self.align_destroy = self.ssw.align_destroy
        self.align_destroy.argtypes = [ct.POINTER(CAlignRes)]
        self.align_destroy.restype = None


class SSW(object):
    def __init__(self, match=4, mismatch=6, gap_open=8, gap_extend=2, lib_path=lib_path):
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.lib_path = lib_path
        self.mat = self.build_matrix()
        self.ssw = CSsw(self.lib_path)

    def build_matrix(self):
        self.dEle2Int = {}
        dInt2Ele = {}
        self.lEle = ['A', 'C', 'G', 'T', 'N']
        for i, ele in enumerate(self.lEle):
            self.dEle2Int[ele] = i
            self.dEle2Int[ele.lower()] = i
            dInt2Ele[i] = ele
        nEleNum = len(self.lEle)
        lScore = [0 for _ in range(nEleNum ** 2)]
        for i in range(nEleNum - 1):
            for j in range(nEleNum - 1):
                if self.lEle[i] == self.lEle[j]:
                    lScore[i * nEleNum + j] = self.match
                else:
                    lScore[i * nEleNum + j] = -self.mismatch
        mat = (len(lScore) * ct.c_int8)()
        mat[:] = lScore
        return mat

    def set_reference_sequence(self, reference):
        self.reference = reference
        self.rNum = self.to_int(self.reference)
        self.reference_len = len(self.reference)

    def to_int(self, seq,):

        num_decl = len(seq) * ct.c_int8
        num = num_decl()
        for i, ele in enumerate(seq):
            try:
                n = self.dEle2Int[ele]
            except KeyError:
                n = self.dEle2Int[self.lEle[-1]]
            finally:
                num[i] = n

        return num

    def get_cigar(self, cigar, ref_position_start, ref_position_end, query_position_start, query_position_end, query):
        sCigarInfo = 'MIDNSHP=X'
        sCigar = ''
        if query_position_start > 0:
            sCigar += str(query_position_start) + 'S'
        for x in cigar:
            n = x >> 4
            m = x & 15
            if m > 8:
                c = 'M'
            else:
                c = sCigarInfo[m]
            if c == 'M': c = '='
            sCigar += str(n) + c
        cigar_l = query_position_end - query_position_start + 1
        if cigar_l < len(query):
            sCigar += str(len(query) - cigar_l) + 'S'

        return sCigar

    def align_one(self, qProfile, rNum, nRLen, nOpen, nExt, nFlag, nMaskLen):

        res = self.ssw.ssw_align(qProfile, rNum, ct.c_int32(nRLen), nOpen, nExt, nFlag, 0, 0, nMaskLen)

        nScore = res.contents.nScore
        nScore2 = res.contents.nScore2
        nRefBeg = res.contents.nRefBeg
        nRefEnd = res.contents.nRefEnd
        nQryBeg = res.contents.nQryBeg
        nQryEnd = res.contents.nQryEnd
        nRefEnd2 = res.contents.nRefEnd2
        lCigar = [res.contents.sCigar[idx] for idx in range(res.contents.nCigarLen)]
        nCigarLen = res.contents.nCigarLen
        self.ssw.align_destroy(res)

        return (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)

    def align(self, query):
        align_flag = 2
        query_len = len(query)
        qNum = self.to_int(query)
        qProfile = self.ssw.ssw_init(qNum, ct.c_int32(query_len), self.mat, len(self.lEle), 2)
        mask_len = 15 if query_len <= 30 else query_len
        res = self.align_one(qProfile, self.rNum, self.reference_len, self.gap_open, self.gap_extend, align_flag, mask_len)
        align_score, ref_position_start, ref_position_end, query_position_start, query_position_end, hexadecimal_cigar = res[0], res[2], res[3], res[4], res[5], res[8]

        cigar = self.get_cigar(hexadecimal_cigar, ref_position_start, ref_position_end, query_position_start, query_position_end, query)
        return align_score, cigar, ref_position_start


