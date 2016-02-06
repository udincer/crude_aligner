import os
import configuration as c
import cPickle as pickle

def compile_alignments_into_dict():
    print 'compiling alignments into dict...'
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

    return location_to_read_dict

def compile_alignments_into_list():
    print 'compiling alignments into list...'
    directory = '{}/{}/alignments/'.format(c.DATA_PATH, c.DATASET)
    part_max = len(os.listdir(directory))
    fns = ['alignments_{}_{}_part_{}.pkl'.format(c.DATASET, c.KEY_SIZE, file_num) for file_num in xrange(0, part_max)]

    l = []
    for fn in fns:
        alignments = pickle.load(open('{}{}'.format(directory,fn), 'rb'))
        alignments = [item for sublist in alignments for item in sublist]
        l.extend(alignments)
        print 'done {}'.format(fn)

    return l


if __name__ == '__main__':
    d = compile_alignments_into_dict()
    pickle.dump(d, file('compiled_{}.pkl'.format(c.DATASET), 'wb'))

    l = compile_alignments_into_list()
    pickle.dump(l, file('compiled_list_{}.pkl'.format(c.DATASET), 'wb'))
