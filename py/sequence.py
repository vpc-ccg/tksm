#!/usr/bin/env python3

import json
import os
import pathlib
import re
import math
import statistics
import sys
import uuid
import random
import argparse
import functools
from multiprocessing import Pool
import gzip

import numpy as np
from joblib import Parallel, delayed
from sklearn.neighbors import KernelDensity
import scipy.special

import edlib
from tqdm import tqdm

### START OF BADREAD CODE ###
"""

Code in this block (i.e BADREAD CODE block) is part of Badread. It has been slightly modified from the original.
The following is the copyright notice and warranty use from the original code:

Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread
Badread is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Badread is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Badread.
If not, see <http://www.gnu.org/licenses/>.
"""

### START OF SETTINGS_PY ###
# When adding errors to a sequence (sequence_fragment function in simulate.py), it is necessary to
# keep a running estimate of the sequence identity. This estimate can be greatly improved by
# actually aligning the error-added sequence to the original sequence, but this is computationally
# expensive. These settings control how often such alignment takes place and how much sequence is
# used in the alignments. Set ALIGNMENT_INTERVAL to a very large value to turn alignments off, so
# the identity estimates will be based on error count alone.
class SETTINGS_PY:
    ALIGNMENT_INTERVAL = 25
    ALIGNMENT_SIZE = 1000
    # I don't let users set a very small minimum mean read length (e.g. 2) or very low minimum read
    # identity (e.g. 50%) as that might break some things. These settings control how low they can go.
    MIN_MEAN_READ_LENGTH = 100
    MIN_MEAN_READ_IDENTITY = 50
    # If a random qscore model is used, this controls the range of the qualities.
    RANDOM_QSCORE_MIN = 1
    RANDOM_QSCORE_MAX = 20
    # If an ideal qscore model is used, these settings control the range of qualities.
    IDEAL_QSCORE_RANK_1_MIN, IDEAL_QSCORE_RANK_1_MAX = 1, 3
    IDEAL_QSCORE_RANK_2_MIN, IDEAL_QSCORE_RANK_2_MAX = 4, 7
    IDEAL_QSCORE_RANK_3_MIN, IDEAL_QSCORE_RANK_3_MAX = 8, 20
    IDEAL_QSCORE_RANK_4_MIN, IDEAL_QSCORE_RANK_4_MAX = 21, 30
    IDEAL_QSCORE_RANK_5_MIN, IDEAL_QSCORE_RANK_5_MAX = 31, 40
    IDEAL_QSCORE_RANK_6_MIN, IDEAL_QSCORE_RANK_6_MAX = 41, 50
    # Chimeric reads may or may not get adapters in the middle.
    CHIMERA_START_ADAPTER_CHANCE = 0.25
    CHIMERA_END_ADAPTER_CHANCE = 0.25
### END OF SETTINGS_PY ###

### START OF ERROR_MODEL_PY ###
class ERROR_MODEL_PY:
    error_model_names = {'random', 'nanopore2018', 'nanopore2020', 'pacbio2016'}
    class ErrorModel(object):
        def __init__(self, model_path, model_name, output=sys.stderr):
            self.kmer_size = None
            self.alternatives = {}
            self.probabilities = {}
            if model_name == 'random':
                print('\nUsing a random error model', file=output)
                self.type = 'random'
                self.kmer_size = 1
            elif model_name in ERROR_MODEL_PY.error_model_names:
                self.load_from_file(f'{model_path}/{model_name}.error.gz', output)
            else:
                self.load_from_file(model_name, output)

        def load_from_file(self, filename, output):
            print('\nLoading error model from {}'.format(filename), file=output)
            self.type = 'model'
            count = 0
            with MISC_PY.get_open_func(filename)(filename, 'rt') as model_file:
                for line in model_file:
                    kmer = line.split(',', 1)[0]
                    print('\r  ' + kmer, file=output, end='')

                    # All k-mers in the model must be the same size.
                    if self.kmer_size is None:
                        self.kmer_size = len(kmer)
                    else:
                        assert self.kmer_size == len(kmer)

                    alternatives = [x.split(',') for x in line.strip().split(';') if x]
                    assert alternatives[0][0] == kmer

                    self.alternatives[kmer] = [ERROR_MODEL_PY.align_kmers(kmer, x[0]) for x in alternatives]
                    self.probabilities[kmer] = [float(x[1]) for x in alternatives]
                    count += 1
            print(f'\r  done: loaded error distributions for {count} {self.kmer_size}-mers',
                file=output)

        def add_errors_to_kmer(self, kmer):
            """
            Takes a k-mer and returns a (possibly) mutated version of the k-mer, along with the edit
            distance.
            """
            if self.type == 'random':
                return ERROR_MODEL_PY.add_one_random_change(kmer)

            if kmer not in self.alternatives:
                return ERROR_MODEL_PY.add_one_random_change(kmer)

            alts = self.alternatives[kmer]
            probs = self.probabilities[kmer]

            # The model probabilities for alternate k-mers should total to 1 or a bit less than 1.
            # If less, then the remaining probability is given to random change.
            random_change_prob = 1.0 - sum(probs)
            if random_change_prob > 0.0:
                alts.append(None)
                probs.append(random_change_prob)

            alt = random.choices(alts, weights=probs)[0]
            if alt is None:
                return ERROR_MODEL_PY.add_one_random_change(kmer)
            else:
                return alt
    
    @classmethod
    def align_kmers(cls, kmer, alt):
        """
        This function aligns a k-mer to its alternative, returning a list of strings that join to make
        the alternative, but positioned according to the k-mer.
        Examples:
            align_kmers('ACGT', 'ACGTT') -> ['A', 'C', 'GT', 'T']
            align_kmers('ACGT', 'ACT') -> ['A', 'C', '', 'T']
        """
        # The reference k-mer must be at least a 3-mer (but the alt can be a 2-mer).
        assert len(kmer) > 2
        assert len(alt) > 1

        result = [kmer[0]] + [None] * (len(kmer) - 2) + [kmer[-1]]

        # The alternative k-mer should match at the start and end, so we can trim the sequences down by
        # one before aligning them.
        assert kmer[0] == alt[0] and kmer[-1] == alt[-1]
        kmer, alt = kmer[1:-1], alt[1:-1]

        # Edlib doesn't seem to like it when the first sequence is empty.
        if len(alt) == 0:
            cigar = '{}D'.format(len(kmer))
        else:
            cigar = edlib.align(alt, kmer, task='path')['cigar']

        cigar_parts = re.findall(r'\d+[IDX=]', cigar)
        kmer_pos, alt_pos = 0, 0
        for c in cigar_parts:
            cigar_type = c[-1]
            cigar_size = int(c[:-1])
            if cigar_type == '=' or cigar_type == 'X':
                for i in range(cigar_size):
                    result[kmer_pos+1] = alt[alt_pos]
                    alt_pos += 1
                    kmer_pos += 1
            elif cigar_type == 'D':
                for i in range(cigar_size):
                    result[kmer_pos+1] = ''
                    kmer_pos += 1
            else:
                assert cigar_type == 'I'
                result[kmer_pos] += alt[alt_pos:alt_pos+cigar_size]
                alt_pos += cigar_size

        # If the insertion landed on the first base, shift it over to the second base (ensures that
        # the first and last base are always the same).
        if len(result[0]) == 2:
            first_base, inserted_base = result[0]
            result[0] = first_base
            result[1] = inserted_base + result[1]
        return result

    @classmethod
    def add_one_random_change(cls, kmer):
        result = [x for x in kmer]  # Change 'ACGT' to ['A', 'C', 'G', 'T']
        error_type = random.choice(['s', 'i', 'd'])
        error_pos = random.randint(0, len(kmer) - 1)
        if error_type == 's':  # substitution
            result[error_pos] = MISC_PY.get_random_different_base(result[error_pos])
        elif error_type == 'i':  # insertion
            if random.random() < 0.5:
                result[error_pos] = result[error_pos] + MISC_PY.get_random_base()
            else:
                result[error_pos] = MISC_PY.get_random_base() + result[error_pos]
        else:  # deletion
            result[error_pos] = ''
        return result
### END OF ERROR_MODEL_PY ###

### START OF MISC.PY ###
class MISC_PY:
    RANDOM_SEQ_DICT = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    @classmethod
    def get_random_base(cls):
        """
        Returns a random base with 25% probability of each.
        """
        return cls.RANDOM_SEQ_DICT[random.randint(0, 3)]
    @classmethod
    def get_random_different_base(cls, b):
        random_base = cls.get_random_base()
        while b == random_base:
            random_base = cls.get_random_base()
        return random_base
    @classmethod
    def get_random_sequence(cls, length):
        """
        Returns a random sequence of the given length.
        """
        return ''.join([cls.get_random_base() for _ in range(length)])
    @classmethod
    def identity_from_edlib_cigar(cls, cigar):
        matches, alignment_length = 0, 0
        cigar_parts = re.findall(r'\d+[IDX=]', cigar)
        for c in cigar_parts:
            cigar_type = c[-1]
            cigar_size = int(c[:-1])
            alignment_length += cigar_size
            if cigar_type == '=':
                matches += cigar_size
        try:
            return matches / alignment_length
        except ZeroDivisionError:
            return 0.0
    @classmethod
    def get_open_func(cls, filename):
        if cls.get_compression_type(filename) == 'gz':
            return gzip.open
        else:  # plain text
            return open
    @classmethod
    def get_compression_type(cls,filename):
        """
        Attempts to guess the compression (if any) on a file using the first few bytes.
        http://stackoverflow.com/questions/13044562
        """
        magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                    'bz2': (b'\x42', b'\x5a', b'\x68'),
                    'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
        max_len = max(len(x) for x in magic_dict)

        unknown_file = open(filename, 'rb')
        file_start = unknown_file.read(max_len)
        unknown_file.close()
        compression_type = 'plain'
        for file_type, magic_bytes in magic_dict.items():
            if file_start.startswith(magic_bytes):
                compression_type = file_type
        if compression_type == 'bz2':
            sys.exit('Error: cannot use bzip2 format - use gzip instead')
        if compression_type == 'zip':
            sys.exit('Error: cannot use zip format - use gzip instead')
        return compression_type
    @classmethod
    def print_in_two_columns(cls, l1p1, l2p1, l3p1, l1p2, l2p2, l3p2, output, space_between=6):
        part_1_len = max(len(l1p1), len(l2p1), len(l3p1)) + space_between
        format_str = '{:<' + str(part_1_len) + '}'
        l1p1 = format_str.format(l1p1)
        l2p1 = format_str.format(l2p1)
        l3p1 = format_str.format(l3p1)
        print(l1p1 + l1p2, file=output)
        print(l2p1 + l2p2, file=output)
        print(l3p1 + l3p2, file=output)
    @classmethod
    def float_to_str(cls, v, decimals=1, trim_zeros=False):
        if float(int(v)) == v:
            return str(int(v))
        else:
            formatter = '%.' + str(decimals) + 'f'
            result = formatter % v
            if trim_zeros:
                while result.endswith('0'):
                    result = result[:-1]
            return result
### END MISC.PY ###

### START OF SIMULATE_PY ###
class SIMULATE_PY:
    @classmethod
    def sequence_fragment(cls, fragment, target_identity, error_model, qscore_model, tail_model, compute_qscores=True):

        # Buffer the fragment a bit so errors can be added to the first and last bases.
        k_size = error_model.kmer_size
        tail_noise_seq = tail_model.noise_seq(len(fragment))
        fragment = MISC_PY.get_random_sequence(k_size) + fragment +  tail_noise_seq +  MISC_PY.get_random_sequence(k_size)
        frag_len = len(fragment)

        # A list to hold the bases for the errors-added fragment. Note that these values can be ''
        # (meaning the base was deleted) or more than one base (meaning there was an insertion).
        new_fragment_bases = [x for x in fragment]

        errors = 0.0
        change_count, loop_count = 0, 0
        max_kmer_index = len(new_fragment_bases) - 1 - k_size
        while True:
            # A precaution to make sure we don't get caught in an infinite loop.
            loop_count += 1
            if loop_count > 100 * frag_len:
                break

            # If we have changed almost every base in the fragment, then we can give up (the identity
            # is about as low as we can make it). This is likely to only happen when the target
            # identity is very low (below 60%).
            if change_count > 0.9 * frag_len:
                break

            # To gauge the identity, we first use the number of changes we've added to the fragment,
            # which will probably under-estimate the identity, but it's fast.
            estimated_identity = 1.0 - (errors / frag_len)
            if estimated_identity <= target_identity:
                break

            i = random.randint(0, max_kmer_index)
            kmer = fragment[i:i+k_size]
            new_kmer = error_model.add_errors_to_kmer(kmer)

            # If the error model didn't make any changes (quite common with a non-random error model),
            # we just try again at a different position.
            if kmer == ''.join(new_kmer):
                continue

            for j in range(k_size):
                fragment_base = fragment[i+j]
                new_base = new_kmer[j]  # can actually be more than one base, in cases of insertion

                # If this base is changed in the k-mer and hasn't already been changed, then we apply
                # the change.
                if new_base != fragment_base and fragment_base == new_fragment_bases[i+j]:
                    new_fragment_bases[i+j] = new_base
                    change_count += 1
                    if len(new_base) < 2:  # deletion or substitution
                        new_errors = 1
                    else:  # insertion
                        new_errors = len(new_base) - 1

                    # As the identity gets lower, adding errors has less effect (presumably because
                    # adding an error can shift the alignment in a way that makes the overall identity
                    # no worse or even better). So we scale our new error count down a bit using our
                    # current estimate of the identity.
                    errors += new_errors * (estimated_identity ** 1.5)

                    # Every now and then we actually align a piece of the new sequence to its original
                    # to improve our estimate of the read's identity.
                    if change_count % SETTINGS_PY.ALIGNMENT_INTERVAL == 0:

                        # If the sequence is short enough, we align the whole thing and get an exact
                        # identity.
                        if frag_len <= SETTINGS_PY.ALIGNMENT_SIZE:
                            cigar = edlib.align(fragment, ''.join(new_fragment_bases),
                                                task='path')['cigar']
                            actual_identity = MISC_PY.identity_from_edlib_cigar(cigar)
                            errors = (1.0 - actual_identity) * frag_len

                        # If the sequence is longer, we align a random part of the sequence and use
                        # the result to update the error estimate.
                        else:
                            pos = random.randint(0, frag_len - SETTINGS_PY.ALIGNMENT_SIZE)
                            pos2 = pos+SETTINGS_PY.ALIGNMENT_SIZE
                            cigar = edlib.align(fragment[pos:pos2],
                                                ''.join(new_fragment_bases[pos:pos2]),
                                                task='path')['cigar']
                            actual_identity = MISC_PY.identity_from_edlib_cigar(cigar)
                            estimated_errors = (1.0 - actual_identity) * frag_len
                            weight = SETTINGS_PY.ALIGNMENT_SIZE / frag_len
                            errors = (estimated_errors * weight) + (errors * (1-weight))

        start_trim = len(''.join(new_fragment_bases[:k_size]))
        end_trim = len(''.join(new_fragment_bases[-k_size:]))

        seq = ''.join(new_fragment_bases)
        if compute_qscores:
            qual, actual_identity, identity_by_qscores = QSCOREMODEL_PY.get_qscores(seq, fragment, qscore_model)
        else:
            qual = 'K'*len(seq)
            actual_identity = 1.0 - (errors / frag_len)
            identity_by_qscores = None
        assert(len(seq) == len(qual))

        seq = seq[start_trim:-end_trim]
        qual = qual[start_trim:-end_trim]

        return seq, qual, actual_identity, identity_by_qscores
### END OF SIMULATE_PY ###

### START OF QSCORE_MODEL.PY ###
class QSCOREMODEL_PY:
    qscore_model_names = {'random', 'ideal', 'nanopore2018', 'nanopore2020', 'pacbio2016'}
    class QScoreModel(object):
        def __init__(self, model_path, model_name, output=sys.stderr):
            self.scores, self.probabilities = {}, {}
            self.kmer_size = 1
            self.type = None

            if model_name == 'random':
                self.set_up_random_model(output)
            elif model_name == 'ideal':
                self.set_up_ideal_model(output)
            elif model_name in QSCOREMODEL_PY.qscore_model_names:
                self.load_from_file(f'{model_path}/{model_name}.qscore.gz', output)
            else:
                self.load_from_file(model_name, output)

            # These three cigars must be in the model, as they are the simplest 1-mer cigars and
            # without them the get_qscore method can fail.
            assert '=' in self.scores
            assert 'X' in self.scores
            assert 'I' in self.scores

        def set_up_random_model(self, output):
            print('\nUsing a random qscore model', file=output)
            self.type = 'random'
            self.kmer_size = 1
            for c in ['=', 'X', 'I']:
                self.scores[c], self.probabilities[c] = \
                    QSCOREMODEL_PY.uniform_dist_scores_and_probs(SETTINGS_PY.RANDOM_QSCORE_MIN,
                                                SETTINGS_PY.RANDOM_QSCORE_MAX)

        def set_up_ideal_model(self, output):
            print('\nUsing an ideal qscore model', file=output)
            self.type = 'ideal'
            self.kmer_size = 9

            # Lowest quality: mismatches and insertions.
            for c in ['X', 'I']:
                self.scores[c], self.probabilities[c] = \
                    QSCOREMODEL_PY.uniform_dist_scores_and_probs(SETTINGS_PY.IDEAL_QSCORE_RANK_1_MIN,
                                                SETTINGS_PY.IDEAL_QSCORE_RANK_1_MAX)

            # Increasing quality with length of match run
            self.scores['='], self.probabilities['='] = \
                QSCOREMODEL_PY.uniform_dist_scores_and_probs(SETTINGS_PY.IDEAL_QSCORE_RANK_2_MIN,
                                            SETTINGS_PY.IDEAL_QSCORE_RANK_2_MAX)
            self.scores['==='], self.probabilities['==='] = \
                QSCOREMODEL_PY.uniform_dist_scores_and_probs(SETTINGS_PY.IDEAL_QSCORE_RANK_3_MIN,
                                            SETTINGS_PY.IDEAL_QSCORE_RANK_3_MAX)
            self.scores['====='], self.probabilities['====='] = \
                QSCOREMODEL_PY.uniform_dist_scores_and_probs(SETTINGS_PY.IDEAL_QSCORE_RANK_4_MIN,
                                            SETTINGS_PY.IDEAL_QSCORE_RANK_4_MAX)
            self.scores['======='], self.probabilities['======='] = \
                QSCOREMODEL_PY.uniform_dist_scores_and_probs(SETTINGS_PY.IDEAL_QSCORE_RANK_5_MIN,
                                            SETTINGS_PY.IDEAL_QSCORE_RANK_5_MAX)
            self.scores['========='], self.probabilities['========='] = \
                QSCOREMODEL_PY.uniform_dist_scores_and_probs(SETTINGS_PY.IDEAL_QSCORE_RANK_6_MIN,
                                            SETTINGS_PY.IDEAL_QSCORE_RANK_6_MAX)

        def load_from_file(self, filename, output):
            print('\nLoading qscore model from {}'.format(filename), file=output)
            self.type = 'model'
            last_cigar_len = 0
            count = 0
            with MISC_PY.get_open_func(filename)(filename, 'rt') as model_file:
                for line in model_file:
                    parts = line.strip().split(';')
                    try:
                        if parts[0] == 'overall':
                            continue
                        cigar = parts[0]
                        k = len(cigar.replace('D', ''))
                        if k > self.kmer_size:
                            self.kmer_size = k
                        print('\r  ' + cigar + (' ' * (last_cigar_len - len(cigar))),
                            file=output, end='')
                        last_cigar_len = len(cigar)
                        scores_and_probs = [x.split(':') for x in parts[2].split(',') if x]
                        self.scores[cigar] = [int(x[0]) for x in scores_and_probs]
                        self.probabilities[cigar] = [float(x[1]) for x in scores_and_probs]
                        count += 1
                    except (IndexError, ValueError):
                        sys.exit(f'Error: {filename} does not seem to be a valid qscore model file')
                print(f'\r  done: loaded qscore distributions for {count} alignments',
                    file=output)

        def get_qscore(self, cigar):
            """
            If the cigar is in the model, then we use it to choose a qscore. If not, then we trim the
            cigar down by 2 (1 off each end) and try again with the simpler cigar.
            """
            while True:
                assert len(cigar.replace('D', '')) % 2 == 1
                if cigar in self.scores:
                    scores = self.scores[cigar]
                    probs = self.probabilities[cigar]
                    qscore = random.choices(scores, weights=probs)[0]
                    break
                else:
                    cigar = cigar[1:-1].strip('D')
            return QSCOREMODEL_PY.qscore_val_to_char(qscore)

    @classmethod
    def uniform_dist_scores_and_probs(cls, bottom_q, top_q):
        count = top_q - bottom_q + 1
        scores = list(range(bottom_q, top_q + 1))
        probabilities = [1 / count] * count
        return scores, probabilities

    @classmethod
    def get_qscores(cls, seq, frag, qscore_model):
        assert len(seq) > 0

        # TODO: I fear this full sequence alignment will be slow for long and inaccurate sequences.
        #       Can I break it into chunks for better performance?
        cigar = edlib.align(seq, frag, task='path')['cigar']
        actual_identity = MISC_PY.identity_from_edlib_cigar(cigar)

        aligned_seq, aligned_frag, full_cigar = cls.align_sequences_from_edlib_cigar(seq, frag, cigar)
        unaligned_len = len(seq)
        margins = (qscore_model.kmer_size - 1) // 2

        qscores, error_probs = [], []

        seq_pos_to_alignment_pos = {}
        i, j = 0, 0
        for c in full_cigar:
            if c != 'D':
                seq_pos_to_alignment_pos[i] = j
                i += 1
            j += 1

        for i in range(unaligned_len):
            start = i - margins
            end = i + margins
            while start < 0 or end >= unaligned_len:  # pull back to a smaller k-mer near the seq ends
                start += 1
                end -= 1
            start = seq_pos_to_alignment_pos[start]
            end = seq_pos_to_alignment_pos[end]
            partial_cigar = full_cigar[start:end + 1]
            assert not partial_cigar.startswith('D')
            assert not partial_cigar.endswith('D')
            k_size = len(partial_cigar.replace('D', ''))
            assert k_size <= qscore_model.kmer_size
            assert k_size % 2 == 1  # should be an odd length k-mer
            q = qscore_model.get_qscore(partial_cigar)

            qscores.append(q)
            error_probs.append(cls.qscore_char_to_error_prob(q))

        identity_by_qscores = 1.0 - statistics.mean(error_probs)

        return ''.join(qscores), actual_identity, identity_by_qscores

    @classmethod
    def align_sequences_from_edlib_cigar(cls, seq, frag, cigar, gap_char='-'):
        aligned_seq, aligned_frag, full_cigar = [], [], []
        seq_pos, frag_pos = 0, 0
        cigar_parts = re.findall(r'\d+[IDX=]', cigar)
        for c in cigar_parts:
            cigar_type = c[-1]
            cigar_size = int(c[:-1])
            if cigar_type == '=' or cigar_type == 'X':
                aligned_seq.append(seq[seq_pos:seq_pos+cigar_size])
                aligned_frag.append(frag[frag_pos:frag_pos+cigar_size])
                seq_pos += cigar_size
                frag_pos += cigar_size
            elif cigar_type == 'I':
                aligned_seq.append(seq[seq_pos:seq_pos+cigar_size])
                aligned_frag.append(gap_char * cigar_size)
                seq_pos += cigar_size
            elif cigar_type == 'D':
                aligned_seq.append(gap_char * cigar_size)
                aligned_frag.append(frag[frag_pos:frag_pos+cigar_size])
                frag_pos += cigar_size
            full_cigar.append(cigar_type * cigar_size)
        return ''.join(aligned_seq), ''.join(aligned_frag), ''.join(full_cigar)
    @classmethod
    def qscore_char_to_error_prob(cls, q):
        return cls.qscore_val_to_error_prob(cls.qscore_char_to_val(q))
    @classmethod
    def qscore_char_to_val(cls, q):
        return ord(q) - 33

    @classmethod
    def qscore_val_to_char(cls, q):
        return chr(q + 33)

    @classmethod
    def qscore_val_to_error_prob(cls, q):
        return 10.0 ** (-q/10.0)
### END OF QSCORE_MODEL_PY ###

### START OF IDENTITIES_PY ###
class IDENTITIES_PY:
    class Identities(object):
        def __init__(self, mean, stdev, max_identity, output=sys.stderr):
            # Divide by 100 to convert from percentage to fraction
            self.mean = mean / 100.0
            self.stdev = stdev / 100.0
            self.max_identity = max_identity / 100.0
            print('', file=output)

            if self.mean == self.max_identity:
                self.beta_a, self.beta_b = None, None
                print(f'Using a constant read identity of {self.mean * 100}%', file=output)
            elif self.stdev == 0.0:
                self.max_identity = self.mean
                print(f'Using a constant read identity of {self.mean * 100}%', file=output)
            else:  # beta distribution
                print('Generating read identities from a beta distribution:', file=output)
                self.beta_a, self.beta_b = IDENTITIES_PY.beta_parameters(mean, stdev, max_identity)
                MISC_PY.print_in_two_columns(f'  mean  = {MISC_PY.float_to_str(self.mean * 100):>3}%',
                                    f'  max   = {MISC_PY.float_to_str(self.max_identity * 100):>3}%',
                                    f'  stdev = {MISC_PY.float_to_str(self.stdev * 100):>3}%',
                                    'shape parameters:',
                                    f'  alpha = {self.beta_a:.4e}',
                                    f'  beta  = {self.beta_b:.4e}',
                                    output=output)
                QUICKHIST_PY.quickhist_beta(self.beta_a, self.beta_b, self.max_identity, 8, output=output)

        def get_identity(self):
            if self.mean == self.max_identity:
                return self.mean
            else:  # beta distribution
                return self.max_identity * np.random.beta(self.beta_a, self.beta_b)

    @classmethod
    def beta_parameters(cls, beta_mean, beta_stdev, beta_max):
        u, s, m = beta_mean, beta_stdev, beta_max
        beta_a = (((1-(u/m)) / ((s/m)**2)) - (m/u)) * ((u/m)**2)
        beta_b = beta_a * ((m/u) - 1)
        if beta_a < 0.0 or beta_b < 0.0:
            sys.exit('Error: invalid beta parameters for identity distribution - trying increasing '
                    'the maximum identity or reducing the standard deviation')
        return beta_a, beta_b
### END OF IDENTITIES_PY ###

### START OF QUICKHIST_PY ###
class QUICKHIST_PY:
    @classmethod
    def quickhist_beta(cls, a, b, max_identity, height, output=sys.stderr):
        hist_min, hist_max = 50, 100
        tick_interval = 10
        if cls.get_max_width() > 120:
            bin_size = 0.5
        else:
            bin_size = 1
        bins = (np.arange(hist_min, hist_max, bin_size) + bin_size) / 100 / max_identity
        y = []
        for b_val in bins:
            x = b_val - (bin_size / 200)  # take the value from the middle of the bin
            if x < 1:
                # This is the original function to get the density:
                # beta = x**(a-1) * (1-x)**(b-1) / scipy.special.beta(a, b)
                # But to avoid overflows, I had to log-ify it:
                beta = np.exp((a-1)*np.log(x) + (b-1)*np.log(1-x) - scipy.special.betaln(a, b))
            else:
                beta = 0.0
            y.append(beta)
        bins *= max_identity
        shape = (hist_min, hist_max)
        cls.draw_hist(y, shape, len(bins), height, tick_interval, output=output)
    @classmethod
    def get_max_width(cls):
        col_count = cls.get_terminal_size_stderr().columns
        if col_count < 80:
            return 80
        if col_count > 160:
            return 160
        return col_count
    @classmethod
    def get_terminal_size_stderr(cls, fallback=(80, 24)):
        """
        Unlike shutil.get_terminal_size, which looks at stdout, this looks at stderr.
        """
        try:
            size = os.get_terminal_size(sys.__stderr__.fileno())
        except (AttributeError, ValueError, OSError):
            size = os.terminal_size(fallback)
        return size
    @classmethod
    def draw_hist(cls, y, shape, bins, height, x_tick_interval, y_label='', y_label_space=0,
                print_labels=True, output=sys.stderr):
        max_count = max(y)
        normed_hist_list = [float(x) * height / max_count for x in y]

        # Build plot from top level
        i = 0
        for depth in range(height-1, -1, -1):
            if 0 <= i-y_label_space < len(y_label):
                print(y_label[i-2], end='', file=output)
            else:
                print(' ', end='', file=output)
            print(' │', end='', file=output)

            for item in normed_hist_list:
                floored_item = math.floor(item)
                if floored_item >= depth:
                    if floored_item == depth and 0.75 > item % 1 > 0.25:
                        print('\u2596', end='', file=output)  # half bar
                    elif floored_item == depth and item % 1 > 0.75:
                        print('\u258c', end='', file=output)  # full bar
                    elif floored_item > depth:
                        print('\u258c', end='', file=output)  # full bar
                    else:
                        print(' ', end='', file=output)
                else:
                    print(' ', end='', file=output)
            print('', file=output)
            i += 1

        # Draw X axis with labels
        line, labels = '  ', '  '
        label = shape[0]
        bin_size = (shape[1] - shape[0]) / bins
        label_step = int(x_tick_interval * bin_size)
        for i in range(bins+1):
            if i == 0:
                line += '├'
                labels += str(label)
            elif i % x_tick_interval == 0:
                line += '┐' if i == bins else '┬'
                label += label_step
                labels += str(label)
            else:
                line += '─'
                labels += ' ' * (len(line) - len(labels))
        print(line, file=output)
        if print_labels:
            print(labels, file=output)
### END OF QUICKHIST_PY ###

### END OF BADREAD CODE ###


### START OF TAIL_NOISE_MODEL_PY ###
class TAIL_NOISE_MODEL_PY:
    tail_model_names = {'nanopore', 'pacbio', 'no-noise'}
    class KDE_noise_generator:
        @classmethod
        def from_data(cls, mapped, unmapped, lx, ly, transition_matrix, bw, threads=64):
            kde_vals = KernelDensity(bandwidth=bw).fit(list(zip(mapped,unmapped)))
            
            full_grid = np.array( [[x, y] for x in lx for y in ly]).reshape(lx.shape[0],ly.shape[0],2)
            results=Parallel(n_jobs=threads)(delayed(kde_vals.score_samples)(row) for row in full_grid)

            grid = np.exp(np.stack(results))
            ratio = np.sum(np.array(unmapped) > 0) / len(unmapped)

            return cls(lx, ly, grid, transition_matrix, ratio)

        def __init__(self, lx, ly, grid, transition_matrix, ratio, bases=list("AGTC")):
            self.ly = np.array(ly)
            self.lx = np.array(lx)
            self.grid = np.array(grid)
            self.ratio = ratio
            self.begin_W = np.array(transition_matrix[0])
            self.transition_matrix = np.array(transition_matrix[1])
            
            self.disto = TAIL_NOISE_MODEL_PY.Custom2Dist(self.lx, self.ly, self.grid)
            self.BASES = bases
        

        def sample(self, x):
            return self.disto(x)
        

        def noise_seq(self, i):
            if random.uniform(0,1) > self.ratio:
                return ""
            x = self.sample(i)
            if( x == 0):
                return ""
            if (isinstance(x, list) and len(x) == 1):
                x = x[0]
            index = random.choices(np.arange(0,4,dtype=np.int32))[0]
            gen_seq = ["A"]*x
            for i in range(x):
                w = self.transition_matrix[index,:]
                index = random.choices(np.arange(0,4,dtype=np.int32),w)[0]
                gen_seq[i] = self.BASES[index]
            return "".join(gen_seq)
        
        def save(self,file):
            dc = self.disto.serialize()
            dc["bases"] = self.BASES
            dc["ratio"] = self.ratio
            dc["trans"] = self.transition_matrix.tolist()
            dc["begin"] = self.begin_W.tolist()
            json.dump(dc, file)
        @classmethod
        def load(cls, model_path, model_name):
            if model_name in ["no_noise", "pacbio"]:
                return TAIL_NOISE_MODEL_PY.Mock_noise_generator()
            elif model_name in TAIL_NOISE_MODEL_PY.tail_model_names:
                model_name = f'{model_path}/{model_name}.tail.gz'
            dc = json.load(
                gzip.open(model_name,'rt') 
                if model_name.endswith(".gz")
                else open(model_name,'r')
            )
            return cls(dc["lx"], dc["ly"], dc["grid"], (dc["begin"], dc["trans"]), dc["ratio"], dc["bases"])

    class Mock_noise_generator:
        def __init__(self):
            pass
        def sample(self, x):
            return 0
        def noise_seq(self, i):
            return ""

    class CustomDist:
        def __init__(self, pdf, indices, ):
            self.pdf = pdf
            self.indices = indices
            
            sum_pdf = np.sum(self.pdf)
            self.cdf = np.zeros(len(self.pdf))
            for i, v in enumerate(self.pdf):
                self.cdf[i] = v / sum_pdf + self.cdf[i-1]
            
        def __call__(self):
            val = random.uniform(0,1)
            pos = np.searchsorted(self.cdf, val)
            return self.indices[pos]

    class Custom2Dist:
        @classmethod
        def from_file(cls, path):
            lx = np.load(f'{path}/labelsx.npy')
            ly = np.load(f'{path}/labelsy.npy')
            grid = np.load(f'{path}/grid.npy').T
            return cls(lx,ly,grid)
        @classmethod
        def from_dict(cls, dc):
            lx = np.array(dc["lx"])
            ly = np.array(dc["ly"])
            grid = np.array(dc["grid"])
            return cls(lx,ly,grid)
        def __init__(self, lx, ly, grid):
            self.lx = lx
            self.ly = ly
            self.grid = grid
            self.distributions = []
            for i, y in enumerate(self.ly):
                disto = TAIL_NOISE_MODEL_PY.CustomDist(self.grid[i,:], self.lx)
                self.distributions.append(disto)
        def serialize(self):
            return {"lx" : self.lx.tolist(), 
                    "ly" : self.ly.tolist(),
                    "grid" : self.grid.tolist()}
        def save(self, file):
            json.dump( self.serialize(), file)
        @classmethod
        def load(cls, file):
            return cls.from_dict(json.load(file))
        def __call__(self, y):
            pos = np.searchsorted(self.ly,y)
            if pos < len(self.ly) - 1:
                if np.abs(self.ly[pos] - y) > np.abs(self.ly[pos+1] - y):
                    pos+=1
            if pos >= len(self.ly):
                mult = pos / self.ly[-1]
                pos = len(self.ly) - 1
            else:
                mult = 1
            return int(self.distributions[pos]() * mult)
        def __getitem__(self, y): #returns the subprobability dist without calling
            pos = np.searchsorted(self.ly,y)
            if pos < len(self.ly) - 1:
                if np.abs(self.ly[pos] - y) > np.abs(self.ly[pos+1] - y):
                    pos+=1
            if pos >= len(self.ly):
                mult = y / self.ly[-1]
                pos = len(self.ly) - 1
            else:
                mult = 1
            return lambda : int(mult * self.distributions[pos]())
### END OF TAIL_NOISE_MODEL_PY ###

def parse_args():
    parser = argparse.ArgumentParser(
        description="Sequencer module of RNAInFuser. Generates FASTQ file from MDF (molecule description file) and reference FASTA files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i",
                        "--input",
                        type=str,
                        required=True,
                        help="MDF file.")
    parser.add_argument("-r",
                        "--references",
                        type=str,
                        nargs="+",
                        required=True,
                        help="Reference FASTA files, space separated.")
    parser.add_argument("-o",
                        "--badread",
                        type=str,
                        help="Badread reads output file. Default: no output.")
    parser.add_argument("--perfect",
                        type=str,
                        help="Perfect reads output file. Default: no output.")
    parser.add_argument("--skip-qual-compute",
                        default=False,
                        action="store_true",
                        help="Use all 'K' for q-scores. Default: compute quality score using Edlib alignment to input sequences.")
    parser.add_argument("-O",
                        "--output-format",
                        type=str,
                        default=None,
                        choices=['fastq', 'fasta'],
                        help="Output format. Default: infer from output file extension and fallback to FASTA. Accepts .gz files too.")
    parser.add_argument('-t',
                        '--threads',
                        type=int,
                        default=1,
                        help="Number of threads.",
                        )
    parser.add_argument('--badread-model-path',
                        type=str,
                        default='.',
                        help='Path to the directory containing the Badread models.'
                        )
    parser.add_argument('--badread-identity',
                        type=str,default='87.5,97.5,5',
                        help='Sequencing identity distribution (mean, max and stdev)')
    parser.add_argument('--badread-error-model',
                        type=str,
                        default='nanopore2020',
                        choices=ERROR_MODEL_PY.error_model_names,
                        # help='Can be "nanopore2018", "nanopore2020", "pacbio2016", "random" or a model filename'
                        )
    parser.add_argument('--badread-qscore-model',
                        type=str,
                        default='nanopore2020',
                        choices=QSCOREMODEL_PY.qscore_model_names,
                        # help='Can be "nanopore2018", "nanopore2020", "pacbio2016", "random", "ideal" or a model filename'
                        )
    parser.add_argument('--badread-tail-model',
                        type=str,
                        default='no_noise',
                        choices=TAIL_NOISE_MODEL_PY.tail_model_names,
                        # help='Can be "nanopore", "pacbio", "no_noise" or a model filename')
                        )
    args = parser.parse_args()

    # Process arguments and check for errors
    try:
        identity_parameters = [float(x) for x in args.badread_identity.split(',')]
        assert len(identity_parameters) == 3, "Must specify 3 values for --badread-identity"
        args.badread_mean_identity = identity_parameters[0]
        args.badread_max_identity = identity_parameters[1]
        args.badread_identity_stdev = identity_parameters[2]
    except (ValueError, IndexError):
        sys.exit('Error: could not parse --identity values')
    if args.badread_mean_identity > 100.0:
        sys.exit('Error: mean read identity cannot be more than 100')
    if args.badread_max_identity > 100.0:
        sys.exit('Error: max read identity cannot be more than 100')
    if args.badread_mean_identity <= SETTINGS_PY.MIN_MEAN_READ_IDENTITY:
        sys.exit(f'Error: mean read identity must be at least {SETTINGS_PY.MIN_MEAN_READ_IDENTITY}')
    if args.badread_max_identity <= SETTINGS_PY.MIN_MEAN_READ_IDENTITY:
        sys.exit(f'Error: max read identity must be at least {SETTINGS_PY.MIN_MEAN_READ_IDENTITY}')
    if args.badread_mean_identity > args.badread_max_identity:
        sys.exit(f'Error: mean identity ({args.badread_mean_identity}) cannot be larger than max '
                 f'identity ({args.badread_max_identity})')
    if args.badread_identity_stdev < 0.0:
        sys.exit('Error: read identity stdev cannot be negative')

    if not args.badread and not args.perfect:
        parser.error("Must specify either --output or --perfect-reads.")
    return args

def generate_fasta(fasta):
    name = ""
    seq = list()
    if fasta.endswith('.gz'):
        f = gzip.open(fasta, 'rt')
    else:
        f = open(fasta, 'r')
    for _, l in enumerate(f):
        l = l.rstrip('\n')
        if l[0]=='>':
            if len(seq) == 0:
                name = l[1:].split(' ')[0]
                continue
            yield (name, ''.join(seq))
            seq = list()
            name = l[1:].split(' ')[0]
        else:
            seq.append(l)
    yield (name, ''.join(seq))

def get_reference_seqs(reference):
    reference_seqs = dict()
    for ref in reference:
        print(f"Loading reference {ref}...")
        reference_seqs.update({
            name: seq for name, seq in generate_fasta(ref)
        })
    return reference_seqs


def mdf_generator(f):
    # MDF format:
    # Each entry starts with a header: +<mol_id>\t<depth>\t<comments>; 
    #   where <mol_id> is the molecule ID, <depth> is the depth of the molecule, and each <comments> entry is key=value pair terminated by a semicolon.
    # Each line after the header, and until the next header or end of file, is an interval description: <ref_id>\t<start:int>\t<end:int>\t<strand>\t<modifications>
    read_id = None
    intervals = list()

    for line in f:
        line = line.strip('\n').split('\t')
        if line[0][0] == '+':
            if read_id is not None:
                for _ in range(depth):
                    yield read_id, intervals
            read_id = line[0][1:]
            depth = int(line[1])
            intervals = list()
        else:
            chrom, start, end, strand, modifications = line

            start, end = int(start), int(end)
            intervals.append(
                (chrom, start, end, strand, modifications)
            )
    if read_id is not None:
        for _ in range(depth):
            yield read_id, intervals

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def apply_modifications(seq, modifications):
    if modifications == '':
        return seq
    else:
        modifications = modifications.split(',')
        seq = list(seq)
        for mod in modifications:
            char = mod[-1]
            pos = int(mod[:-1])
            seq[pos] = char
        return ''.join(seq)


def badread(molecule_id, raw_seq, compute_qual=False):
    target_identity = identities.get_identity()
    (
        seq, 
        quals, 
        actual_identity, 
        _
    ) = SIMULATE_PY.sequence_fragment(
        raw_seq,
        target_identity,
        error_model,
        qscore_model,
        tail_model,
        compute_qual,
    )
    info = [
        f"length={len(seq)}",
        f"error_free_length={len(raw_seq)}",
        f"read_identity={actual_identity * 100.0:.2f}%",
        f"molecule_id={molecule_id}",
    ]
    return (seq, quals, info)

def perfect(molecule_id, seq):
    quals = 'K' * len(seq)
    info = f"length={len(seq)} molecule_id={molecule_id}"
    info = [
        f"length={len(seq)}",
        f"error_free_length={len(seq)}",
        f"read_identity={1 * 100.0:.2f}%",
        f"molecule_id={molecule_id}",
    ]
    return (seq, quals, info)

def fastq_formatter(read_id, seq, quals, info):
    result = [
        f"@{read_id} {' '.join(info)}",
        seq,
        '+',
        quals,
    ]
    return '\n'.join(result)

def fasta_formatter(read_id, seq, _, info):
    result = [
        f">{read_id} {' '.join(info)}",
        seq,
    ]
    return '\n'.join(result)

def get_output_file(outfile):
    if outfile.endswith('.gz'):
        f = gzip.open(outfile, 'wt+')
        outfile = outfile[:-3]
    else:
        f = open(outfile, 'w+')
    if outfile.endswith('.fastq') or outfile.endswith('.fq'):
        return f, fastq_formatter
    else:
        return f, fasta_formatter

def mdf_to_seq(mdf, targets=dict()):
    molecule_id, intervals = mdf
    seq = list()
    for chrom, start, end, strand, modifications in intervals:
        segment = reference_seqs[chrom][start:end].upper()
        segment = apply_modifications(segment, modifications)
        if strand == '+':
            seq.append(segment)
        else:
            seq.append(reverse_complement(segment))
    seq = ''.join(seq)

    results = dict()
    read_id = uuid.UUID(int=random.getrandbits(128))
    for k,(seq_producer, read_formatter) in targets.items():
        seq, quals, info = seq_producer(molecule_id, seq)
        results[k] = read_formatter(read_id, seq, quals, info)        
    return results

if __name__ == "__main__":
    args = parse_args()
    reference_seqs = get_reference_seqs(args.references)

    targets = dict()
    target_outfiles = dict()
    if args.badread:
        identities = IDENTITIES_PY.Identities(args.badread_mean_identity, args.badread_identity_stdev, args.badread_max_identity, sys.stderr)
        error_model = ERROR_MODEL_PY.ErrorModel(args.badread_model_path, args.badread_error_model, sys.stderr)
        qscore_model = QSCOREMODEL_PY.QScoreModel(args.badread_model_path, args.badread_error_model, sys.stderr)
        tail_model = TAIL_NOISE_MODEL_PY.KDE_noise_generator.load(args.badread_model_path, args.badread_tail_model)

        output_file, output_file_formatter = get_output_file(args.badread)
        target_outfiles['badread'] = output_file
        partial_badread = functools.partial(badread, compute_qual=(not args.skip_qual_compute) and output_file_formatter == fastq_formatter)
        targets['badread']= [partial_badread, output_file_formatter]
    if args.perfect:
        output_file, output_file_formatter = get_output_file(args.perfect)
        target_outfiles['perfect'] = output_file
        targets['perfect']= [perfect, output_file_formatter]

    mdf_to_seq_targeted = functools.partial(mdf_to_seq, targets=targets)
    with open(args.input, 'r') as f:
        mdg = mdf_generator(f)
        if args.threads > 1:
            mapper = functools.partial(Pool(args.threads).imap_unordered, chunksize=10)
        else:
            mapper = map
        for read_dict in tqdm(mapper(mdf_to_seq_targeted, mdg), desc="Sequencing"):
            for k,v in read_dict.items():
                print(v, file=target_outfiles[k])

        if args.threads > 1:
            Pool(args.threads).close()

    for v in target_outfiles.values():
        v.close()
