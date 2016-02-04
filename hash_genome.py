from collections import defaultdict
import cPickle as pickle
import utils
import configuration as c


def read_file_hash_and_pickle(file_path, key_size, pickle_filename):
    ref = ''
    with open(file_path, 'r') as f:
        first_line = True
        line_count = 0
        for line in f:
            if first_line:
                first_line = False
                continue
            ref += line.strip()
            if line_count%100000==0:
                print 'lines read: {}'.format(line_count)
            line_count += 1

    print 'initializing dict'
    ref_hash = defaultdict(list)
    for i in xrange(len(ref)-key_size+1):
        key = ref[i:i+key_size]
        intkey = utils.key_to_integer(key)
        ref_hash[intkey].append(i)
        if i%100000==0:
            print 'hashing position {}'.format(i)

    print 'pickling hashed dict'
    pickle.dump(ref_hash, open(pickle_filename, 'wb'))


if __name__ == '__main__':
    print 'hashing the genome...'

    file_path = 'data/{}/{}'.format(c.DATASET, c.REF_FILE)
    pickle_filename = 'ref_hash_{}_{}.pkl'.format(c.DATASET, c.KEY_SIZE)
    read_file_hash_and_pickle(file_path, c.KEY_SIZE, pickle_filename)
