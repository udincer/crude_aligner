import numpy as np

def edit_distance_matrix(ref, donor):
    """
    Computes the edit distance matrix between the donor and reference
    This algorithm makes substitutions, insertions, and deletions all equal.
    Does that strike you as making biological sense? You might try changing the cost of
    deletions and insertions vs snps.
    :param ref: reference genome (as an ACTG string)
    :param donor: donor genome guess (as an ACTG string)
    :return: complete (len(ref) + 1) x (len(donor) + 1) matrix computing all changes
    """

    output_matrix = np.zeros((len(ref), len(donor)))
    # print len(ref), len(donor)
    # print output_matrix
    # This is a very fast and memory-efficient way to allocate a matrix
    for i in range(len(ref)):
        output_matrix[i, 0] = i
    for j in range(len(donor)):
        output_matrix[0, j] = j
    for j in range(1, len(donor)):
        for i in range(1, len(ref)):  # Big opportunities for improvement right here.
            deletion = output_matrix[i - 1, j] + 1
            insertion = output_matrix[i, j - 1] + 1
            identity = output_matrix[i - 1, j - 1] if ref[i] == donor[j] else np.inf
            substitution = output_matrix[i - 1, j - 1] + 1 if ref[i] != donor[j] else np.inf
            output_matrix[i, j] = min(insertion, deletion, identity, substitution)
    return output_matrix


def identify_changes(ref, donor, offset):
    """
    Performs a backtrace-based re-alignment of the donor to the reference and identifies
    SNPS, Insertions, and Deletions.
    Note that if you let either sequence get too large (more than a few thousand), you will
    run into memory issues.
    :param ref: reference sequence (ATCG string)
    :param donor: donor sequence (ATCG string)
    :param offset: The starting location in the genome.
    :return: SNPs, Inserstions, and Deletions
    """
    # print offset
    ref = '${}'.format(ref)
    donor = '${}'.format(donor)
    edit_matrix = edit_distance_matrix(ref=ref, donor=donor)
    current_row = len(ref) - 1
    current_column = len(donor) - 1
    score = edit_matrix[current_row, current_column]
    changes = []
    while current_row > 0 or current_column > 0:
        if current_row == 0:
            pvs_row = -np.inf
        else:
            pvs_row = current_row - 1

        if current_column == 0:
            pvs_column = -np.inf
        else:
            pvs_column = current_column - 1

        try:
            insertion_dist = edit_matrix[current_row, pvs_column]
        except IndexError:
            insertion_dist = np.inf

        try:
            deletion_dist = edit_matrix[pvs_row, current_column]
        except IndexError:
            deletion_dist = np.inf

        try:
            if ref[current_row] == donor[current_column]:
                identity_dist = edit_matrix[pvs_row, pvs_column]
            else:
                identity_dist = np.inf

            if ref[current_row] != donor[current_column]:
                substitution_dist = edit_matrix[pvs_row, pvs_column]
            else:
                substitution_dist = np.inf
        except (TypeError, IndexError) as e:
            identity_dist = np.inf
            substitution_dist = np.inf

        min_dist = min(insertion_dist, deletion_dist, identity_dist, substitution_dist)

        ref_index = current_row + offset - 1
        if min_dist == identity_dist:
            current_row = pvs_row
            current_column = pvs_column
        elif min_dist == substitution_dist:
            changes.append(['SNP', ref[current_row], donor[current_column], ref_index + 1])
            current_row = pvs_row
            current_column = pvs_column
        elif min_dist == insertion_dist:
            if len(changes) > 0 and changes[-1][0] == 'INS' and changes[-1][-1] == ref_index + 1:
                changes[-1][1] = donor[current_column] + changes[-1][1]
            else:
                changes.append(['INS', donor[current_column], ref_index + 1])
            current_column = pvs_column
        elif min_dist == deletion_dist:
            if len(changes) > 0 and changes[-1][0] == 'DEL' and changes[-1][-1] == ref_index + 1:
                changes[-1] = ['DEL', ref[current_row] + changes[-1][1], ref_index]
            else:
                changes.append(['DEL', ref[current_row], ref_index + 1])
            current_row = pvs_row
        else:
            raise ValueError
    changes = sorted(changes, key=lambda change: change[-1])
    #print str(changes)
    return changes, score

if __name__ == '__main__':
    a ='TCGGATAAGGCACCTGAATTATGAACTTGTCAAGTCGGCACCCGCGGGCTTTGCCACCGATGAACGCCGTTAAAGCCCAGACATTAAAGGACGGTGACGTGT'
    b ='TCGGATAAGTCACCTGAATTATGAACTTGTCAAGTCGGCACCCGCGTGGCTTTGCCACCGATGAACGCCGTTAAAGCCCAGACATTAAAGGACGGTGACTGT'
    changes, score = identify_changes(a,b,0)
    pass