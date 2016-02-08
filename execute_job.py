import cPickle as pickle
import configuration as c
import pileup
import utils
import dependencies
import msgpack
import sys

if __name__ == '__main__':

    args = sys.argv[1:]
    if len(args) ==  1:
        job_number = int(args[0])
    else:
        job_number = 1

    directory = 'jobs'

    print 'loading {}...'.format(job_number)
    sub, d = pickle.load(file('{}/job_{}_part_{}.pkl'.format(
        directory, c.DATASET, job_number),'rb'))
    print 'loaded {}'.format(job_number)

    ref = utils.read_reference('data/{}/{}'.format(c.DATASET, c.REF_FILE))

    the_donors = pileup.generate_donor_pieces(sub, ref, d)
    the_donors = [donors for donors in the_donors if donors]

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

    utils.write_indels(good_changes, job_number)
    msgpack.dump(the_donors,file('donors_{}_part_{}'.format(c.DATASET, job_number),'wb'))