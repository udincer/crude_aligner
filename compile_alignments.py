import os
import configuration as c
import cPickle as pickle
import time
import utils
from bitarray import bitarray
import msgpack

# def compile_alignments_into_dict():
#     print 'compiling alignments into dict...'
#     directory = '{}/{}/alignments/'.format(c.DATA_PATH, c.DATASET)
#     part_max = len(os.listdir(directory))
#     fns = ['alignments_{}_{}_part_{}.pkl'.format(c.DATASET, c.KEY_SIZE, file_num) for file_num in xrange(0, part_max)]
#
#     location_to_read_dict = {}
#     for fn in fns:
#         alignments = pickle.load(open('{}{}'.format(directory,fn), 'rb'))
#         alignments = [item for sublist in alignments for item in sublist]
#         for t in alignments:
#             location_to_read_tdict[t[0]] = t[1]
#         print 'done {}'.format(fn)
#
#     return location_to_read_dict

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


def ident_bad_regions(some_pr_tuples, ref):
    nonperfect_pos_read_tuples = []
    good_areas = bitarray(len(ref))
    good_areas.setall(False)
    start = time.time()
    count = 0.0
    for t in some_pr_tuples:
        read_num = t[1]
        pos = t[0]
        ref_piece_num = utils.key_to_integer(ref[pos:pos+c.READ_SIZE])
        if read_num == ref_piece_num:
            good_areas[pos:pos+c.READ_SIZE] = True
        else:
            pass
        #good_areas[pos:pos+c.READ_SIZE] = True
        #else:
        #    nonperfect_pos_read_tuples.append(t)
        count += 1
        if count % 10000 == 0:
            print '{}'.format(count/len(some_pr_tuples))
            print '{}'.format(time.time() - start)
    print 'piece done in {}'.format((time.time() - start))
    return (good_areas, nonperfect_pos_read_tuples)


def merge_bad_regions(good_areas_list):
    r = bitarray(good_areas_list[0])
    for i in xrange(1,len(good_areas_list)):
        r = r | good_areas_list[i]
    return r


def get_nonperfect_stretches(ref):
    #directory = '{}/{}/alignments/'.format(c.DATA_PATH, c.DATASET)
    directory = 'alignments/'

    part_max = len(os.listdir(directory))
    fns = ['alignments_{}_{}_part_{}.pkl'.format(c.DATASET, c.KEY_SIZE, file_num) for file_num in xrange(0, part_max)]

    ga_overall = None
    first = True
    #nprt_list = []
    for fn in fns:
        print 'getting stretches for: {}'.format(fn)
        alignments = pickle.load(open('{}{}'.format(directory,fn), 'rb'))
        alignments = [item for sublist in alignments for item in sublist]
        (ga, nprt) = ident_bad_regions(alignments, ref)
        if first:
            ga_overall = ga
            first = False
        else:
            ga_overall = merge_bad_regions([ga_overall, ga])
        #nprt_list.extend(nprt)

    stretches = []
    start = 0
    for i in xrange(1,len(ga_overall)):
        if ga_overall[i-1] and not ga_overall[i]:
            start = i
        elif not ga_overall[i-1] and ga_overall[i]:
            end = i
            stretches.append((start, end))

    return stretches


# def compile_alignments_into_sorted_list():
#     print 'compiling alignments into SORTED list...'
#     directory = '{}/{}/alignments/'.format(c.DATA_PATH, c.DATASET)
#     part_max = len(os.listdir(directory))
#     fns = ['alignments_{}_{}_part_{}.pkl'.format(c.DATASET, c.KEY_SIZE, file_num) for file_num in xrange(0, part_max)]
#
#     l = SortedListWithKey(key=lambda val:val[0])
#     for fn in fns:
#         alignments = pickle.load(open('{}{}'.format(directory,fn), 'rb'))
#         alignments = [item for sublist in alignments for item in sublist]
#         l.update(alignments)
#         print 'done {}'.format(fn)
#     return l



if __name__ == '__main__':



    ref = utils.read_reference()
    stretches = get_nonperfect_stretches(ref)
    print stretches

    #print 'sorting...'
    #sl = sorted(nprt_list)
    #print 'sorted: {}'.format(sl)

    pickle.dump(stretches, file('stretches_{}.pkl'.format(c.DATASET), 'wb'))
    #msgpack.dump(stretches, file('stretches_{}.msg'.format(c.DATASET), 'wb'))

    #pickle.dump(sl, file('sorted_nprt_{}.pkl'.format(c.DATASET), 'wb'))

    print 'DONE'

    # start = time.time()
    # print 'compiling sorted list...'
    # sl = compile_alignments_into_sorted_list()
    # pickle.dump(sl, file('compiled_sorted_list_{}.pkl'.format(c.DATASET), 'wb'))
    # print 'elapsed time for sorted list {}'.format(time.time()-start)

    # print 'compiling dict...'
    # d = compile_alignments_into_dict()
    # pickle.dump(d, file('compiled_{}.pkl'.format(c.DATASET), 'wb'))
    # print 'elapsed time for sorted list {}'.format(time.time()-start)


    # start = time.time()
    # print 'compiling regular list...'
    # l = compile_alignments_into_list()
    # pickle.dump(l, file('compiled_list_{}.pkl'.format(c.DATASET), 'wb'))
    # print 'elapsed time for sorted list {}'.format(time.time()-start)
