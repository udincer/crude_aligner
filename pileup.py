import configuration as c
import utils
import cPickle as pickle
from collections import Counter
import numpy as np
from distance import hamming, levenshtein
import time
import sys
import os

CALL_LOOKUP = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
INV_LOOKUP  = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def consensus(res):
    # possible_indel_locations = []
    donor = ['X'] * c.OUTPUT_SIZE
    for i in xrange(len(res)):
        if max(res[i]) < c.MIN_COVERAGE_FOR_INCLUSION: # ultra low coverage
            donor[i] = 'X'
            continue
        call_index = res[i].argmax()
        donor[i] = CALL_LOOKUP[call_index]
        # if res[i][call_index] > 3 and sum(res[i])/res[i][call_index]-1 > 0.2:
        #     possible_indel_locations.append(i)
    # print 'possible {}'.format(possible_indel_locations)
    return donor

def pileup2(read_pos_tuple_list):
    result = np.zeros((c.OUTPUT_SIZE, 4))
    #result = [[0,0,0,0] for x in xrange(c.OUTPUT_SIZE)]
    count = 0.0
    for read_pos_tuple in read_pos_tuple_list:
        read = utils.integer_to_key(read_pos_tuple[1],c.KEY_SIZE)
        position = read_pos_tuple[0]

        for i in xrange(position,position+len(read)):
            j = i - position
            if read[j]=='A':
                result[i][0] += 1
            elif read[j]=='C':
                result[i][1] += 1
            elif read[j]=='G':
                result[i][2] += 1
            elif read[j]=='T':
                result[i][3] += 1
                #result[position:position+len(read)] = read
        count += 1
        if count % 100000 == 0:
            print 'done: {:.2f}'.format(count/len(read_pos_tuple_list))

    donor = consensus(result)
    return ''.join(donor)

def get_alignment_from_pickle(pickle_fn):
    print 'loading alignment pickle'
    mapping = pickle.load(open(pickle_fn, 'rb'))
    print 'alignment pickle loaded'
    return mapping


def all_snp(donor, ref):
    snp_list = [[ref[i],donor[i],str(i)] if ref[i] != donor[i] else [] for i in xrange(min(len(donor),len(ref)))]
    snp_list_noempty = [x for x in snp_list if x]

    snp_list_filtered = [snp_list_noempty[i] for i in xrange(len(snp_list_noempty)-1) if int(snp_list_noempty[i+1][2])-int(snp_list_noempty[i][2]) > 2]
    # TODO print stuff
    with open("{}_out.txt".format(c.DATASET), "w") as output_file:
        output_file.write('>{}_{}\n'.format(c.DATASET, 'chr_1'))
        output_file.write('>SNP')
        for snp in snp_list_filtered:
            if snp[1]!='X':
                output_file.write('\n{},{},{}'.format(snp[0],snp[1],snp[2]))

    return snp_list_filtered

#
# def generate_donor_pieces_with_sorted_list(stretches, ref, sorted_nprt):
#     STRETCH_LIMIT = 120
#     MARGIN = 2 + c.READ_SIZE/4
#     count = 0.0
#     the_donors = []
#     donor_indexes = []
#     for stretch in stretches[0:100]:
#         stretch_length = stretch[1]-stretch[0]
#         if stretch_length > STRETCH_LIMIT:
#             print '{} is over stretch limit, skipping.'.format(stretch)
#             the_donors.append(None)
#             continue
#         donor = ['.']*(c.READ_SIZE + 2*MARGIN + stretch_length)
#         read_tuples = []
#         (start, end) = (stretch[0]-c.READ_SIZE, stretch[1]+c.READ_SIZE)
#         if start < 0 or end > len(ref):
#             continue
#
#         for i in xrange(start,end):
#             try:
#                 read_tuples.append((i,utils.sorted_get_index(sorted_nprt,i)))
#             except ValueError:
#                 pass
#         # print 'reads: {}'.format(read_tuples)
#         distances = []
#         #for read_tuple in read_tuples:
#             #ref_piece = ref[read_tuple[0]-MARGIN:read_tuple[0]+c.READ_SIZE+MARGIN]
#             # read_str = utils.integer_to_key(read[1], c.READ_SIZE)
#             #read_str = read_tuple
#             #distances.append(utils.sliding_window(read_str, ref_piece))
#             # print 'ref: {}\nrea: {}'.format(ref_piece, read_str)
#             # print 'distances {}'.format(distances)
#
#         # seed generation!
#         #argmin = distances.index(min(distances))
#
#         argmin = 0 # first one always behaves well!
#         #pos = read_tuples[argmin][0]
#         str = read_tuples[argmin][1]
#         donor[argmin:argmin+c.READ_SIZE] = list(utils.integer_to_key(str,c.READ_SIZE))
#
#         iteration_count = xrange(5)
#         to_be_removed = []
#         threshold = -40
#         for _ in iteration_count:
#             threshold += 2
#             for item in to_be_removed:
#                 try:
#                     read_tuples.remove(item)
#                 except ValueError:
#                     print 'Value not in list problem.'
#                     break
#             to_be_removed = []
#
#             chosen_ones = []
#
#             #print '\n{} -> {}'.format(''.join(donor), stretch)
#             for read_tuple in read_tuples:
#                 read_num = read_tuple[1]
#                 read = utils.integer_to_key(read_num, c.READ_SIZE)
#
#                 for offset in xrange(len(donor) - len(read) + 1):
#                     j = len(donor) - len(read) - offset
#                     pre =  '.'*offset
#                     post = '.'*j
#                     padded = pre + read + post
#                     ham = hamming_ignore_dots(donor, padded)
#                     if ham < threshold:
#                         #print '{} -> {}'.format(padded, ham)
#                         chosen_ones.append(padded)
#                         to_be_removed.append(read_tuple)
#
#             piece_of_donor = pileup_ignore_dots(chosen_ones, donor)
#             donor = piece_of_donor # new seed!
#             #print donor
#         the_donors.append(donor)
#         donor_indexes.append(start)
#         count += 1
#         print 'Generating donors: {:.2f}% complete'.format(count/len(stretches)*100)
#         print '{}\n{}'.format(ref[start-20:start+100], '                    ' + donor)
#     return (the_donors, donor_indexes)
#

def generate_donor_pieces(stretches, ref, pos_to_read):
    STRETCH_LIMIT = 200
    MARGIN_LEFT = c.READ_SIZE + 1
    MARGIN_RIGHT = 5
    count = 0.0
    the_donors = []
    donor_indexes = []
    for stretch in stretches:
        stretch_length = stretch[1]-stretch[0]
        if stretch_length > STRETCH_LIMIT:
            print '{} is over stretch limit, skipping.'.format(stretch)
            the_donors.append(None)
            continue
        donor = ['.']*(c.READ_SIZE + 1*MARGIN_LEFT + stretch_length)
        read_tuples = []
        (start, end) = (stretch[0]-MARGIN_LEFT, stretch[1]+MARGIN_RIGHT)
        if start < 0 or end > len(ref):
            continue

        for i in xrange(start,end):
            try:
                read_tuples.append((i,utils.integer_to_key(pos_to_read[i],c.READ_SIZE)))
            except KeyError:
                pass
        # print 'reads: {}'.format(read_tuples)
        distances = []
        #for read_tuple in read_tuples:
            #ref_piece = ref[read_tuple[0]-MARGIN:read_tuple[0]+c.READ_SIZE+MARGIN]
            # read_str = utils.integer_to_key(read[1], c.READ_SIZE)
            #read_str = read_tuple
            #distances.append(utils.sliding_window(read_str, ref_piece))
            # print 'ref: {}\nrea: {}'.format(ref_piece, read_str)
            # print 'distances {}'.format(distances)

        # seed generation!
        #argmin = distances.index(min(distances))

        argmin = 0 # first one always behaves well!
        pos = read_tuples[argmin][0]
        str = read_tuples[argmin][1]
        donor[argmin:argmin+c.READ_SIZE] = list(str)

        iteration_count = xrange(10)
        to_be_removed = []
        threshold = -30
        for _ in iteration_count:
            threshold += 2
            for item in to_be_removed:
                try:
                    read_tuples.remove(item)
                except ValueError:
                    print 'Value not in list problem.'
                    break
            to_be_removed = []

            chosen_ones = []

            #print '\n{} -> {}'.format(''.join(donor), stretch)
            for read_tuple in read_tuples:
                read = read_tuple[1]

                for offset in xrange(len(donor) - len(read) + 1):
                    j = len(donor) - len(read) - offset
                    pre =  '.'*offset
                    post = '.'*j
                    padded = pre + read + post
                    ham = hamming_ignore_dots(donor, padded)
                    if ham < threshold:
                        #print '{} -> {}'.format(padded, ham)
                        chosen_ones.append(padded)
                        to_be_removed.append(read_tuple)

            piece_of_donor = pileup_ignore_dots(chosen_ones, donor)
            donor = piece_of_donor # new seed!
            #print donor
        the_donors.append(donor.strip('.'))
        donor_indexes.append(pos)
        count += 1
        #matching_positions = ['|' if ref[pos+i] == donor[i] else ' ' for i in xrange(len(donor))]
        pipes = ''.join(['|' if ref[pos+i] == donor[i] else ' ' for i in xrange(len(donor))])
        if pipes.count('|') < len(donor.strip('.')) - 5:
            print 'Generating donors: {:.2f}% complete'.format(count/len(stretches)*100)
            print '{}\n{}\n{}'.format(ref[pos:pos+len(donor)], pipes, donor)
    return (the_donors, donor_indexes)



        # TODO

def pileup_ignore_dots(list_good_reads, seed):
    SEED_ADVANTAGE = 4
    donor = seed
    result = np.zeros((len(seed), 4))
    for k in xrange(len(seed)):
        if seed[k] == '.':
            continue
        result[k][INV_LOOKUP[seed[k]]] += SEED_ADVANTAGE

    for good_read in list_good_reads:
        for i in xrange(len(donor)):
            if good_read[i] == '.':
                continue
            result[i][INV_LOOKUP[good_read[i]]] += 1
    donor = consensus_for_donor(result)
    return ''.join(donor)


def consensus_for_donor(result):
    INCLUSION_THRESH = 2

    donor = ['.'] * len(result)
    for i in xrange(len(result)):
        if max(result[i]) < INCLUSION_THRESH:
            donor[i] = '.'
            continue
        call_index = result[i].argmax()
        donor[i] = CALL_LOOKUP[call_index]
    return donor


def hamming_ignore_dots(s1,s2):
    assert len(s1)==len(s2)
    return sum([0 if (s1[i] == '.' or s2[i] == '.') else -1 if s1[i] == s2[i] else 1 for i in xrange(len(s1))])


### MIKE ###

def edit_distance_matrix(ref, donor):
    """
    Computes the edit distance matrix between the donor and reference
    This algorithm makes substitutions, insertions, and deletions all equal.
    Does that strike you as making biological sense? You might try changing the cost of
    deletions and insertions vs snps.
    :param ref: reference genome (as an ACTG string)
    :param donor: donor genome guess (as an ACTG string)
    :return: complete (len(ref) + 1) x (len(donor) + 1) matrix computing all changes
    """

    output_matrix = np.zeros((len(ref), len(donor)))
    # print len(ref), len(donor)
    # print output_matrix
    # This is a very fast and memory-efficient way to allocate a matrix
    for i in range(len(ref)):
        output_matrix[i, 0] = i
    for j in range(len(donor)):
        output_matrix[0, j] = j
    for j in range(1, len(donor)):
        for i in range(1, len(ref)):  # Big opportunities for improvement right here.
            deletion = output_matrix[i - 1, j] + 1
            insertion = output_matrix[i, j - 1] + 1
            identity = output_matrix[i - 1, j - 1] if ref[i] == donor[j] else np.inf
            substitution = output_matrix[i - 1, j - 1] + 1 if ref[i] != donor[j] else np.inf
            output_matrix[i, j] = min(insertion, deletion, identity, substitution)
    return output_matrix


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) == 2:
        FILE_INDEX_BEGIN = int(args[0])
        FILE_INDEX_END = int(args[1])
    else:
        FILE_INDEX_BEGIN = 0
        FILE_INDEX_END = len(os.listdir('data/{}/reads_split/'.format(c.DATASET)))
    print 'Processing files {} through {}'.format(FILE_INDEX_BEGIN, FILE_INDEX_END)

    print 'start pileup'

    ref = utils.read_reference('data/{}/{}'.format(c.DATASET, c.REF_FILE))
    #mapping_pickles_folder = 'data/{}/pickled_mappings'.format(c.DATASET)
    #pickle_files = ['{}/mappings_part_{}.txt.pkl'.format(mapping_pickles_folder, i) for i in xrange(FILE_INDEX_BEGIN, FILE_INDEX_END)]

    alignment = get_alignment_from_pickle('compiled_list_{}.pkl'.format(c.DATASET))
    # alignment_sorted = sorted(alignment)
    # alignment_flat = [item for sublist in alignment for item in sublist]
    #
    stretches = pickle.load(file('stretches_{}.pkl'.format(c.DATASET), 'rb'))
    # nprt_sl = pickle.load(file('sorted_nprt_{}.pkl'.format(c.DATASET), 'rb'))
    #
    # print 'nonpeft len: {}'.format(len(nprt_sl))
    # print 'align len:   {}'.format(len(alignment))
    # pos_to_read = pickle.load(file('compiled_{}.pkl'.format(c.DATASET)))
    #
    # pickle.dump(bad_stretches, file('bad_stretches.pkl','wb'))
    # pickle.dump(pos_to_read, file('pos_to_read.pkl','wb'))
    #
    # print 'pickled!'

    # bad_stretches = pickle.load(file('bad_stretches.pkl','rb'))
    pos_to_read = pickle.load(file('pos_to_read.pkl','rb'))

    # for key in pos_to_read:
    #     print 'key: {}\tvalue: {}'.format(key, pos_to_read[key])

    (the_donors, donor_indexes) = generate_donor_pieces(stretches, ref, pos_to_read)
    #(the_donors, donor_indexes) = generate_donor_pieces_with_sorted_list(stretches, ref, alignment_sorted)

    pickle.dump(the_donors, file('bad_stretches.pkl','wb'))
    pickle.dump(donor_indexes, file('pos_to_read.pkl','wb'))

    # snp_list = []
    # ins_list = []
    # del_list = []
    # for i in xrange(len(donor_indexes)):
    #     stretch = stretches[i]
    #     (start, end) = (donor_indexes[i], donor_indexes[i] + 180)
    #     local_ref = ref[start:end]
    #     if the_donors[i] == None:
    #         continue
    #     donor = the_donors[i].strip('.')
    #     changes = identify_changes(local_ref,donor,0)
    #     for change in changes:
    #         if change[0] == 'SNP':
    #             f = change[1]
    #             t = change[2]
    #             loc = change[3] + donor_indexes[i]
    #             snp_list.append((f,t,loc))
    #         elif change[0] == 'INS':
    #             ins = change[1]
    #             loc = change[2] + donor_indexes[i]
    #             ins_list.append((ins,loc))
    #         if change[0] == 'DEL':
    #             d = change[1]
    #             loc = change[2] + donor_indexes[i]
    #             del_list.append((d,loc))
    #
    # print snp_list
    # print del_list
    # print ins_list


    print 'done?'
    # print 'piling up...'
    # donor = pileup2(alignment_flat)
    #
    # co = Counter()
    # [co.update(d) for d in donor]
    # print co
    #
    # snps = all_snp(donor, ref)
    # print 'done'
