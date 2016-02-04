import configuration as c
import utils
import cPickle as pickle
from collections import Counter

def consensus(res):
    donor = ['X'] * c.OUTPUT_SIZE
    for i in xrange(len(res)):
        if res[i] == [0,0,0,0]:
            call = 'X'
            continue
        elif max(res[i]) < c.MIN_COVERAGE_FOR_INCLUSION:
            # ultra low coverage
            call = 'X'
            continue
        call_index = res[i].index(max(res[i]))
        if call_index == 0:
            call = 'A'
        elif call_index == 1:
            call = 'C'
        elif call_index == 2:
            call = 'G'
        elif call_index == 3:
            call = 'T'
        donor[i] = call
    return donor

def pileup2(read_pos_tuple_list):
    result = [[0,0,0,0] for x in xrange(c.OUTPUT_SIZE)]
    count = 0.0
    for read_pos_tuple in read_pos_tuple_list:
        read = read_pos_tuple[1]
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
            print 'done: {}'.format(count/len(read_pos_tuple_list))

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

    snp_list_filtered = [snp_list_noempty[i] for i in xrange(len(snp_list_noempty)-1) if int(snp_list_noempty[i+1][2])-int(snp_list_noempty[i][2]) > 3]
    # TODO print stuff
    with open("{}_out.txt".format(c.DATASET), "w") as output_file:
        output_file.write('>{}_{}\n'.format(c.DATASET, 'chr1'))
        output_file.write('>SNP')
        for snp in snp_list_filtered:
            if snp[1]!='X':
                output_file.write('\n{},{},{}'.format(snp[0],snp[1],snp[2]))

    return snp_list_filtered


if __name__ == '__main__':
    print 'start pileup'

    ref = utils.read_reference('data/{}/{}'.format(c.DATASET,c.REF_FILE))
    #mapping_pickles_folder = 'data/{}/pickled_mappings'.format(c.DATASET)
    #pickle_files = ['{}/mappings_part_{}.txt.pkl'.format(mapping_pickles_folder, i) for i in xrange(FILE_INDEX_BEGIN, FILE_INDEX_END)]

    alignment = get_alignment_from_pickle('alignments_{}_{}.pkl'.format(c.DATASET, c.KEY_SIZE))
    alignment_flat = [item for sublist in alignment for item in sublist]

    print 'piling up...'
    donor = pileup2(alignment_flat)

    co = Counter()
    [co.update(d) for d in donor]
    print co

    snps = all_snp(donor, ref)
    print 'done'
