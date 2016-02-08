import cPickle as pickle
import configuration as c
import os
import msgpack


if __name__ == '__main__':

    TOTAL_JOB_COUNT = 10

    print 'Generating work orders...'
    directory = 'jobs'
    if not os.path.exists(directory):
        os.makedirs(directory)

    print 'Loading stretches'
    stretches = msgpack.load(file('stretches_{}.msg'.format(c.DATASET), 'rb'))

    print 'loading huge dictionary'
    pos_to_read = pickle.load(file('compiled_{}.pkl'.format(c.DATASET),'rb'))

    for job_number in xrange(1,TOTAL_JOB_COUNT+1):
        job_s = len(stretches)//TOTAL_JOB_COUNT*(job_number-1)
        job_e = len(stretches)//TOTAL_JOB_COUNT*(job_number)

        sub = stretches[job_s:job_e]
        d = {k:v for k,v in pos_to_read.iteritems() if sub[0][0] <= k <= sub[-1][1]}

        pickle.dump((sub, d),file('{}/job_{}_part_{}.pkl'.format(
            directory, c.DATASET, job_number),'wb'))
        print 'dumping {}'.format(job_number)