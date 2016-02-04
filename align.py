import cPickle as pickle
import configuration as c
from collections import Counter
import utils
from multiprocessing import Pool, Manager
import time
import os

def load_hash_pickle(pickle_fn):
    ref_hash = pickle.load(open(pickle_fn, 'rb'))
    return ref_hash


def align_single_read(read, ref_hash, key_size,
                      min_num_pieces = c.MINIMUM_NUMBER_OF_PIECES_TO_MATCH):
    read_pieces = [read[i * key_size: (i + 1) * key_size] for i in xrange(len(read) / key_size)]
    read_piece_locations = [ref_hash[utils.key_to_integer(read_pieces[i])]
                            for i in xrange(len(read_pieces))]
    start_positions = [[x - i * key_size for x in read_piece_locations[i]]
                       for i in xrange(len(read_piece_locations))]

    start_position_counter = Counter()

    for start_position in start_positions:
        start_position_counter.update(start_position)

    viable_locations = []
    for location in start_position_counter:
        if start_position_counter[location] >= min_num_pieces:
            viable_locations.append(location)

    return viable_locations


def align_paired_read(paired_read, ref, ref_hash):
    NUM_PIECES_TO_MATCH_FOR_PRIMARY = 3
    NUM_PIECES_TO_MATCH_FOR_SECONDARY = 1

    orientations = [(paired_read[0], paired_read[1][::-1]),
                    (paired_read[1], paired_read[0][::-1])]

    location_read_tuple_list = []

    for oriented_paired_read in orientations:
        primary_read = oriented_paired_read[0]
        secondary_read = oriented_paired_read[1]

        primary_locations = align_single_read(primary_read, ref_hash, c.KEY_SIZE,
                                              NUM_PIECES_TO_MATCH_FOR_PRIMARY)

        secondary_locations = align_single_read(secondary_read, ref_hash, c.KEY_SIZE,
                                                NUM_PIECES_TO_MATCH_FOR_SECONDARY)

        for pl in primary_locations:
            secondary_location_f = [x for x in secondary_locations
                        if  pl + c.PAIRED_GAP_ESTIMATE_INTERVAL[0] < x < pl + c.PAIRED_GAP_ESTIMATE_INTERVAL[1] or
                            pl - c.PAIRED_GAP_ESTIMATE_INTERVAL[1] < x < pl - c.PAIRED_GAP_ESTIMATE_INTERVAL[0]]
            if len(secondary_location_f) > 0:
                lrt_primary = (pl, primary_read)
                lrt_secondary = (secondary_location_f[0], secondary_read) # only take first one because i just can't
                location_read_tuple_list.append(lrt_primary)
                location_read_tuple_list.append(lrt_secondary)

    return location_read_tuple_list


def parallel_align(paired_read_list, ref, ref_hash):
    results = [align_paired_read(paired_read, ref, ref_hash) for paired_read in paired_read_list]


def align_and_pickle_mappings():
    print 'loading ref_hash pickle...'
    pickle_filename = 'ref_hash_{}_{}.pkl'.format(c.DATASET, c.KEY_SIZE)
    ref_hash = load_hash_pickle(pickle_filename)

    print 'loading ref...'
    ref_file = c.REF_FILE
    ref = utils.read_reference('{}/{}/{}'.format(c.DATA_PATH, c.DATASET, ref_file))

    print 'loading reads...'
    reads_file = c.READS_FILE
    reads = utils.read_reads('{}/{}/{}'.format(c.DATA_PATH, c.DATASET, reads_file))
    print 'loaded reads'

    print 'aligning {}'.format(c.DATASET)
    start = time.time()
    alignments = [align_paired_read(pr, ref, ref_hash) for pr in reads]
    alignments = [x for x in alignments if x]
    print 'alignment complete, elapsed: {}'.format(time.time() - start)

    directory = '{}/{}/{}'.format(c.DATA_PATH, c.DATASET, 'alignments/')
    if not os.path.exists(directory):
        print 'creating folder {}'.format(directory)
        os.makedirs(directory)

    start = time.time()
    alignments_pickle_path = 'alignments_{}_{}.pkl'.format(c.DATASET, c.KEY_SIZE)
    pickle.dump(alignments, open(alignments_pickle_path, 'wb'))
    print '{} alignment pickled in {}'.format(reads_file, time.time() - start)


if __name__ == '__main__':
    align_and_pickle_mappings()
