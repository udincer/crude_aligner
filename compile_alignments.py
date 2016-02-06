import os
import configuration as c
import cPickle as pickle
from collections import defaultdict

def compile_alignments_and_sort():
    directory = '{}/{}/alignments/'.format(c.DATA_PATH, c.DATASET)
    part_max = len(os.listdir(directory))
    fns = ['alignments_{}_{}_part_{}.pkl'.format(c.DATASET, c.KEY_SIZE, file_num) for file_num in xrange(0, part_max)]

    location_to_read_dict = {}
    for fn in fns:
        alignments = pickle.load(open('{}{}'.format(directory,fn), 'rb'))
        alignments = [item for sublist in alignments for item in sublist]
        for t in alignments:
            location_to_read_dict[t[0]] = t[1]
        print 'done {}'.format(fn)

    print location_to_read_dict
    pass

if __name__ == '__main__':
    compile_alignments_and_sort()
