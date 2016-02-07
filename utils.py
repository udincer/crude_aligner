import configuration as c
from bisect import bisect_left

nucleotide_to_binary = {'A':'00', 'T':'01', 'C':'10', 'G':'11'}
binary_to_nucleotide = {'00':'A', '01':'T', '10':'C', '11':'G'}

def read_reference(ref_fn = '{}/{}/{}'.format(c.DATA_PATH,c.DATASET,c.REF_FILE)):
    print 'reading ref...'
    f = open(ref_fn, 'r')
    first_line = True
    output_reference = ''
    for line in f:
        if first_line:
            first_line = False
            continue  # We skip the first line, since it
            # only contains the name of the chromosome the reads
            # came from.
        line = line.strip()
        output_reference += line  # We append each line to the output reference string.
    print 'end ref read'
    return output_reference


def read_reads(read_fn):
    print 'reading reads...'
    f = open(read_fn, 'r')
    all_reads = []
    for line in f:
        if line[0] == '>':
            continue
        line = line.strip()
        paired_end_reads = line.split(',')  # The two paired ends are separated by a comma
        if len(paired_end_reads[1])!=c.READ_SIZE:
            continue
        all_reads.append((paired_end_reads[0],paired_end_reads[1]))
    print 'read reads end'
    return all_reads


def key_to_integer(key):
    binary_str = ''
    for char in key:
        binary_str += nucleotide_to_binary[char]
    return int(binary_str,2)


def integer_to_key(integer, key_size):
    key = ''
    str = format(integer,'b')
    if len(str)%2 == 1:
        str = '0' + str
    kk = [str[i:i+2] for i in xrange(0, len(str), 2)]
    for k in kk:
        key += binary_to_nucleotide[k]
    if len(key)< key_size:
        key = 'A'*(key_size-len(key))+key
    return key


def count_mismatches(s1, s2):
    if len(s1)!=len(s2):
        assert False
    return sum([1 if s1[i]!=s2[i] else 0 for i in xrange(len(s1))])


def sliding_window(s1, s2):
    if len(s1) == len(s2):
        return (0, count_mismatches(s1,s2))
    elif len(s2) < len (s1):
        temp = s1
        s1 = s2
        s2 = temp
    assert len(s1) < len(s2)
    corr = [count_mismatches(s1, s2[i:i+len(s1)]) for i in xrange(len(s2) - len(s1)+1)]
    min_distance = min(corr)
    min_index = corr.index(min_distance)
    return (min_index, min_distance)


def sorted_get_index(list, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left([l[0] for l in list], x)
    if i != len(list) and list[i][0] == x:
        return i
    raise ValueError


   # return sum([0 if s1[i]==s2[i] or (s1[i] == '.' or s2[i] == '.') else 1 for i in xrange(len(s1))])

if __name__ == '__main__':

    print hamming_ignore_dots('......AAABBAAAAAA...','...AAAAAAAAA.....A..')
