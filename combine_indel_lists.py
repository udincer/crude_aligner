import configuration as c
import os

if __name__ == '__main__':
    snps = []
    inss = []
    dels = []

    directory = 'indels_list'
    files = ['{}/{}_indels_out_part_{}.txt'.format(directory, c.DATASET, i) for i in xrange(1, len(os.listdir(directory)))]
    files = ['{}/hw2grad_M_1_indels_out_part_{}.txt'.format(directory, i) for i in xrange(1, len(os.listdir(directory)))]

    for fn in files:
        first_line = True
        with open(fn, 'r') as f:
            mode = ''
            for line in f:
                if line[0] == '>':
                    if first_line:
                        first_line = False
                        continue
                    mode = line[1:4]
                elif mode == 'SNP':
                    s1, s2, pos = line.strip().split(',')
                    if s2 != '.' and pos > 100:
                        snps.append((s1,s2,pos))
                elif mode == 'INS':
                    s, pos = line.strip().split(',')
                    if '.' not in s and pos > 100:
                        inss.append((s,pos))
                elif mode == 'DEL':
                    s, pos = line.strip().split(',')
                    if '.' not in s and pos > 100:
                        dels.append((s, pos))

    with open('total_output.txt','w') as out:
        out.write('>{}'.format(c.DATASET))
        out.write('\n>INS')
        for t in inss:
            out.write('\n{},{}'.format(t[0],t[1]))
        out.write('\n>DEL')
        for t in inss:
            out.write('\n{},{}'.format(t[0],t[1]))
        out.write('\n>SNP')
        for t in snps:
            out.write('\n{},{},{}'.format(t[0],t[1],t[2]))


