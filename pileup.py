import configuration as c
import utils
import cPickle as pickle
from collections import Counter
import numpy as np
import time
import sys
import os
import msgpack
import bisect
#import cProfile
import dependencies

CALL_LOOKUP = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
INV_LOOKUP  = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

SCORE_THRESHOLD = 40


def call_indel_zones(donor, ref):
    BIN_SIZE = 25
    indel_zones = []
    mismatches = [i for i in xrange(len(donor)) if donor[i] != ref[i]]
    for s in xrange(0,len(donor),BIN_SIZE):
        f = bisect.bisect_left(mismatches, s)
        l = bisect.bisect_left(mismatches, s + BIN_SIZE)
        num_mismatches_in_interval = l - f
        if 2 < num_mismatches_in_interval < 8:
            if l == len(mismatches):
                continue
            indel_zones.append((mismatches[f],mismatches[l]))
    return indel_zones




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
        read = utils.integer_to_key(read_pos_tuple[1],c.READ_SIZE)
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
    count = 0.0
    the_donors = []
    for stretch in stretches:
        the_donors.append(get_donor_for_stretch(stretch,ref ,pos_to_read))
        count += 1
        if count % 100 == 0:
            print 'progress: {:.2f}'.format(count / len(stretches))
    return the_donors


def get_donor_for_stretch(stretch, ref, pos_to_read):
    STRETCH_LIMIT = 20
    MARGIN_LEFT = c.READ_SIZE
    MARGIN_RIGHT = stretch[1] - stretch[0] + 8

    stretch_length = stretch[1] - stretch[0]
    if stretch_length > STRETCH_LIMIT:
        print '{} is over stretch limit, skipping.'.format(stretch)
        return
    donor = ['.'] * (MARGIN_RIGHT + MARGIN_LEFT + stretch_length + 1)
    read_tuples = []
    (start, end) = (stretch[0] - MARGIN_LEFT, stretch[1] + MARGIN_RIGHT)
    if start < 0 or end > len(ref):
        return

    for i in xrange(start, end-c.READ_SIZE): # we don't want the extras on the right
        try:
            read_tuples.append((i, utils.integer_to_key(pos_to_read[i], c.READ_SIZE)))
        except KeyError:
            pass
            # print 'reads: {}'.format(read_tuples)
            # distances = []
            # for read_tuple in read_tuples:
            # ref_piece = ref[read_tuple[0]-MARGIN:read_tuple[0]+c.READ_SIZE+MARGIN]
            # read_str = utils.integer_to_key(read[1], c.READ_SIZE)
            # read_str = read_tuple
            # distances.append(utils.sliding_window(read_str, ref_piece))
            # print 'ref: {}\nrea: {}'.format(ref_piece, read_str)
            # print 'distances {}'.format(distances)

    # seed generation!
    # argmin = distances.index(min(distances))
    if len(read_tuples) < stretch_length:
        print 'skipping {} low read tuple count'.format(stretch)
        return
    elif len(read_tuples) > 250:
        print 'skipping {} HIGH read tuple count'.format(stretch)
        return

    # ham = []
    # for s in read_tuples:
    #     rr = s[1]
    #     po = s[0]
    #     ham.append(hamming_ignore_dots_list_of_char(ref[po:po + c.READ_SIZE], rr))
    # argmin = ham.index(min(ham))
    # print 'ARGMIN:{}.'.format(argmin)

    # SEED NUMBER 1
    #argmin = 0  # first one always behaves well!
    # try:
    # #    pos = read_tuples[argmin][0]
    # #    str = read_tuples[argmin][1]
    #     #if hamming_ignore_dots_list_of_char(ref[pos:pos + c.READ_SIZE], str) > -1 * c.READ_SIZE + 1:
    #     if sum([ref[pos+i] == str[i] for i in xrange(len(str))])<49:
    #         print 'skipping due to bad initial read'
    #         #print pos
    #         return
    donor[0:0 + c.READ_SIZE] = list(ref[start:start + c.READ_SIZE])
    #donor[-1*c.READ_SIZE-1:-1] = list(ref[end - c.READ_SIZE:end])
    # except IndexError:
    #     return


    # argmax = -1
    # try:
    #     pos = read_tuples[argmax][0]
    #     str = read_tuples[argmax][1]
    #     if sum([ref[pos+i] == str[i] for i in xrange(len(str))])<50:
    #         print 'skipping due to bad initial read'
    #         return
    #     donor[argmax:argmax + c.READ_SIZE] = list(str)
    # except IndexError:
    #     return

    #print 'initial state of donor:\n{}'.format(''.join(donor))

    iteration_count = xrange(6)
    to_be_removed = []
    still_unused = []
    threshold = -40
    for _ in iteration_count:
        threshold += 3
        for item in to_be_removed:
            try:
                read_tuples.remove(item)
            except ValueError:
                print 'Value not in list problem. repetitive region'
                return
        to_be_removed = []

        chosen_ones = []
        for read_tuple in read_tuples:
            if read_tuple == None:
                continue
            read = read_tuple[1]
            hams = []
            for offset in xrange(0, len(donor) - len(read)):
                j = len(donor) - offset - len(read)
                pre = ['.'] * offset
                post = ['.'] * j
                padded = pre + list(read) + post
                ham = hamming_ignore_dots_list_of_char(donor, padded)
                # if ham < -49:
                #     print 'repetitive region! skipping...'
                #     to_be_removed.extend(read_tuples)
                #     break

                hams.append(ham)
            if min(hams) < threshold:
                offset = hams.index(min(hams))
                j = len(donor) - offset - len(read)
                pre = ['.'] * offset
                post = ['.'] * j
                padded = pre + list(read) + post
                chosen_ones.append(padded)
                to_be_removed.append(read_tuple)
                #print '{} -> {}'.format(''.join(padded), min(hams))

        piece_of_donor = pileup_ignore_dots(chosen_ones, donor)
        donor = piece_of_donor  # new seed!


    #print '\n{} -> {}'.format(''.join(donor), stretch)
    return (start, donor.strip('.'))


def visualize_lines(donor, pos, ref, stretch = ''):
    try:
        pipes = ''.join(['|' if ref[pos+i] == donor[i] else ' ' for i in xrange(len(donor))])
        print '\n{}\n{} -> {}\n{}'.format(ref[pos:pos+len(donor)], pipes, stretch, donor)
    except TypeError:
        pass

def pileup_ignore_dots(list_good_reads, seed):
    SEED_ADVANTAGE = 3
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


def hamming_ignore_dots(s1,s2): # BUG
    assert len(s1)==len(s2)
    r1 = min(s1.find('.'), s2.find('.')) # BUG w find
    r2 = max(s1.rfind('.'), s2.rfind('.'))
    return sum([s1[i]!=s2[i] for i in xrange(r1,r2)])

    #return sum([0 if (s1[i] == '.' or s2[i] == '.') else -1 if s1[i] == s2[i] else 1 for i in xrange(len(s1))])

def hamming_ignore_dots_list_of_char(s1,s2):
    assert len(s1)==len(s2)
    total = 0
    r1 = -1
    r2 = -1
    for i in xrange(len(s1)):
        if s1[i] != '.' and s2[i] != '.':
            r1 = i
            break
    for i in xrange(len(s1)-1,-1,-1):
        if s1[i] != '.' and s2[i] != '.':
            r2 = i
            break
    for i in xrange(r1,r2):
        total += s1[i]!=s2[i]
    return total - r2 + r1

def last_part():
    the_donors = msgpack.load(file('donors_{}'.format(c.DATASET),'rb'))

    #clean up
    the_donors = [donors for donors in the_donors if donors]
    good_changes = []
    for donor_tuple in the_donors:
        donor = donor_tuple[1]
        pos = donor_tuple[0]
        ref_piece = ref[pos:pos+len(donor)]
        changes, score = dependencies.identify_changes(ref_piece,donor,0)
        if score < 10:
            #print visualize_lines(donor,pos,ref_piece)
            for cc in changes:
                try:
                    cc[3] += pos
                except:
                    cc[2] += pos

            good_changes.extend(changes)

    utils.write_indels(good_changes, DO_NUMBER)


if __name__ == '__main__':




    args = sys.argv[1:]
    if len(args) == 2:
        DIVIDE_INTO = int(args[0])
        DO_NUMBER = int(args[1])
    else:
        DIVIDE_INTO = 1
        DO_NUMBER = 1
    print 'Doing part {} of {}'.format(DO_NUMBER, DIVIDE_INTO)

    print 'start pileup'

    ref = utils.read_reference('data/{}/{}'.format(c.DATASET, c.REF_FILE))

    #mapping_pickles_folder = 'data/{}/pickled_mappings'.format(c.DATASET)
    #pickle_files = ['{}/mappings_part_{}.txt.pkl'.format(mapping_pickles_folder, i) for i in xrange(FILE_INDEX_BEGIN, FILE_INDEX_END)]

    alignment = get_alignment_from_pickle('compiled_list_{}.pkl'.format(c.DATASET))
    # alignment_sorted = sorted(alignment)
    # alignment_flat = [item for sublist in alignment for item in sublist]
    #

    if DO_NUMBER == 0:
        print 'all snps'
        d = pileup2(alignment)
        all_snp(d, ref)
        exit()

    #d = pileup2(alignment)
    #msgpack.dump(d, file('pileup2.msg','wb'))
    d = msgpack.load(file('pileup2.msg','rb'))
    print 'pileup complete'
    if DO_NUMBER == 0:
        print 'all snps'
        all_snp(d, ref)
        exit()

    LENGTH_FILTER = 10
    print 'calling indel zones'
    #iz = call_indel_zones(d,ref)
    print 'called indel zones'
    iz = msgpack.load(file('stretches_{}.msg'.format(c.DATASET), 'rb'))
    print 'before length filtering: {}'.format(len(iz))
    iz = [x for x in iz if x[1] - x[0] <= LENGTH_FILTER]
    print '{} spots will be checked for indels'.format(len(iz))


    pos_to_read = pickle.load(file('compiled_{}.pkl'.format(c.DATASET),'rb'))
    #pos_to_read = pickle.load(file('pos_to_read.pkl','rb'))
    # msgpack.dump(d, file('pos_to_read.msg','wb'))
    # pos_to_read = msgpack.load(file('pos_to_read.msg','rb'))

    # cProfile.run('generate_donor_pieces(iz[0:100],ref,pos_to_read)')

    iz_piece = iz[len(iz)//DIVIDE_INTO*(DO_NUMBER-1):len(iz)//DIVIDE_INTO*DO_NUMBER]
    print 'now checking {} positions'.format(len(iz_piece))

    the_donors = generate_donor_pieces(iz_piece, ref, pos_to_read)

    #the_donors = msgpack.load(file('donors_{}'.format(c.DATASET),'rb'))
    #clean up
    the_donors = [donors for donors in the_donors if donors]


    # all_snp(d, ref)
    #msgpack.dump(the_donors,file('donors_{}'.format(c.DATASET),'wb'))
    #exit()

    #

    good_changes = []
    for donor_tuple in the_donors:
        donor = donor_tuple[1]
        pos = donor_tuple[0]
        ref_piece = ref[pos:pos+len(donor)]
        changes, score = dependencies.identify_changes(ref_piece,donor,-1)
        if score < len(donor) * 0.4: # thanks leah
            #visualize_lines(donor,pos,ref, pos)
            #print 'C.\n{}'.format(changes)
            for cc in changes:
                try:
                    if cc[2] == '.':
                        changes.remove(cc)
                        continue
                    cc[3] += pos
                except:
                    cc[2] += pos

            good_changes.extend(changes)

    for _ in good_changes:
        print _

    utils.write_indels(good_changes, DO_NUMBER)

    msgpack.dump(the_donors,file('donors_{}'.format(c.DATASET),'wb'))

    print 'done?'