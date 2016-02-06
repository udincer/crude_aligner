import configuration as c
import utils
import cPickle as pickle
from collections import Counter
import numpy as np
from distance import hamming, levenshtein
import time

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


def get_nonperfect_stretches(pos_read_tuples, ref):
    # identify bad regions

    good_areas = [False] * len(ref)
    start = time.time()
    count = 0.0
    for t in pos_read_tuples:
        read = utils.integer_to_key(t[1], c.READ_SIZE)
        pos = t[0]
        ref_piece = ref[pos:pos+c.READ_SIZE]
        if hamming(read, ref_piece) == 0:
            good_areas[pos:pos+c.READ_SIZE] = [True] * c.READ_SIZE
        count += 1
        if count % 10000 == 0:
            print '{}'.format(count/len(pos_read_tuples))
            print '{}'.format(time.time() - start)

        #print 'read: {}\nreff: {}\ndist:{}'.format(read, ref_piece, dist)
    print 'done in {}'.format((time.time() - start))

    # str = ''.join(['O' if good_areas[i] else 'X' for i in xrange(len(good_areas))])
    # for i in xrange(0,12000,120):
    #     print str[i:i+120]

    stretches = []
    start = 0
    for i in xrange(1,len(good_areas)):
        if good_areas[i-1] and not good_areas[i]:
            start = i
        elif not good_areas[i-1] and good_areas[i]:
            end = i
            stretches.append((start, end))

    return stretches


def generate_donor_pieces(stretches, ref, pos_to_read):
    STRETCH_LIMIT = 200
    MARGIN = 8 + c.READ_SIZE/2
    count = 0.0
    the_donors = []
    for stretch in stretches:
        stretch_length = stretch[1]-stretch[0]
        if stretch_length > STRETCH_LIMIT:
            print '{} is over stretch limit, skipping.'.format(stretch)
            the_donors.append(None)
            continue
        donor = ['.']*(c.READ_SIZE + 2*MARGIN + stretch_length)
        read_tuples = []
        (start, end) = (stretch[0]-c.READ_SIZE, stretch[1]+c.READ_SIZE)
        if start < 0 or end > len(ref):
            continue

        for i in xrange(start,end):
            try:
                read_tuples.append((i,utils.integer_to_key(pos_to_read[i],c.READ_SIZE)))
            except KeyError:
                pass
        # print 'reads: {}'.format(read_tuples)
        distances = []
        for read_tuple in read_tuples:
            ref_piece = ref[read_tuple[0]-MARGIN:read_tuple[0]+c.READ_SIZE+MARGIN]
            # read_str = utils.integer_to_key(read[1], c.READ_SIZE)
            read_str = read_tuple
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
        threshold = -40
        for _ in iteration_count:
            threshold += 2
            for item in to_be_removed:
                read_tuples.remove(item)
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
        the_donors.append(donor)
        count += 1
        print 'Generating donors: {}% complete'.format(count/len(stretches))
    return the_donors



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


if __name__ == '__main__':
    print 'start pileup'

    ref = utils.read_reference('data/{}/{}'.format(c.DATASET, c.REF_FILE))
    #mapping_pickles_folder = 'data/{}/pickled_mappings'.format(c.DATASET)
    #pickle_files = ['{}/mappings_part_{}.txt.pkl'.format(mapping_pickles_folder, i) for i in xrange(FILE_INDEX_BEGIN, FILE_INDEX_END)]

    alignment = get_alignment_from_pickle('compiled_list_{}.pkl'.format(c.DATASET))
    # alignment_flat = [item for sublist in alignment for item in sublist]
    #
    # bad_stretches = get_nonperfect_stretches(alignment, ref)
    # pos_to_read = pickle.load(file('compiled_{}.pkl'.format(c.DATASET)))
    #
    # pickle.dump(bad_stretches, file('bad_stretches.pkl','wb'))
    # pickle.dump(pos_to_read, file('pos_to_read.pkl','wb'))
    #
    # print 'pickled!'

    bad_stretches = pickle.load(file('bad_stretches.pkl','rb'))
    pos_to_read = pickle.load(file('pos_to_read.pkl','rb'))

    # for key in pos_to_read:
    #     print 'key: {}\tvalue: {}'.format(key, pos_to_read[key])
    generate_donor_pieces(bad_stretches, ref, pos_to_read)
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
