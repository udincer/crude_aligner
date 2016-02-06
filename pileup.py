import configuration as c
import utils
import cPickle as pickle
from collections import Counter
import numpy as np
from distance import hamming, levenshtein
import time

CALL_LOOKUP = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}


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
    for stretch in stretches:
        (start, end) = (stretches[0], stretches[1])
        reads = pos_to_read[start:end]


if __name__ == '__main__':
    print 'start pileup'

    ref = utils.read_reference('data/{}/{}'.format(c.DATASET,c.REF_FILE))
    #mapping_pickles_folder = 'data/{}/pickled_mappings'.format(c.DATASET)
    #pickle_files = ['{}/mappings_part_{}.txt.pkl'.format(mapping_pickles_folder, i) for i in xrange(FILE_INDEX_BEGIN, FILE_INDEX_END)]

    alignment = get_alignment_from_pickle('alignments_{}_{}.pkl'.format(c.DATASET, c.KEY_SIZE))
    alignment_flat = [item for sublist in alignment for item in sublist]

    bad_stretches = get_nonperfect_stretches(alignment_flat, ref)
    print 'piling up...'
    donor = pileup2(alignment_flat)

    co = Counter()
    [co.update(d) for d in donor]
    print co

    snps = all_snp(donor, ref)
    print 'done'
