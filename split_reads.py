import os
import utils
import configuration as c

MAX_NUM_LINES = 1000000


def split_reads(read_fn):
    max_num_lines = MAX_NUM_LINES
    file_count = 0
    with open(read_fn, 'r') as f:
        write_buffer = ''
        line_count = 0
        first_line = True
        for line in f:
            if first_line:
                first_line = False
                continue
            write_buffer += line.strip() + '\n'
            line_count += 1
            if line_count >= max_num_lines:
                with open('data/{}/reads_split/part_{}.txt'.format(c.DATASET, file_count),'w') as w:
                    w.write(write_buffer)
                    print 'wrote part {}'.format(file_count)
                file_count += 1
                write_buffer = ''
                line_count = 0
        if len(write_buffer) > 0:
            with open('data/{}/reads_split/part_{}.txt'.format(c.DATASET, file_count),'w') as w:
                w.write(write_buffer)
                print 'wrote part {}'.format(file_count)


if __name__ == '__main__':
    print 'splitting reads...'

    directory = '{}/{}/{}'.format(c.DATA_PATH, c.DATASET, 'reads_split/')
    if not os.path.exists(directory):
        print 'creating folder {}'.format(directory)
        os.makedirs(directory)

    ref = utils.read_reference('data/{}/{}'.format(c.DATASET,c.REF_FILE))
    split_reads('data/{}/{}'.format(c.DATASET, c.READS_FILE))

    print 'reads split'
