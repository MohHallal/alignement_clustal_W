import alignement_local
import numpy as np
import constants
def align_glob(tree, seqs):
    """ Computing the multiple alignment """
    global align_global
    align_global = {}
    """ Computing the local aligmanment between all the close sequnces """
    traverse_tree_leafs(tree, seqs) #faire les alignement 2 à 2 entre les séquences les plus proches
    """ Computing the multiple alignment one by until finishing all the sequences """
    while len(list(align_global.keys())) != 1:
        clades = list(align_global.keys())
        for clade in clades:
            algnA = clade
            """ getting the name of the clade that is going to be aligned with algnA (the closest according to the guide tree) """
            algnB = traverse_tree(algnA, tree)
            if algnB in clades:
                """ generating the name of the new alignemnt """
                new_algn = (algnA, algnB)
                """ Computing the algnment """
                align_global[new_algn] = calculs_align_glob(align_global[algnA], align_global[algnB])
                del align_global[algnA]
                del align_global[algnB]
    return align_global[tree]
def traverse_tree_leafs(tree,seqs):
    """ The function that computes the local alignment between the pairs that are close according to the guide tree """
    if isinstance(tree[0], tuple):
        for subtree in tree:
            traverse_tree_leafs(subtree,seqs)
    elif len(tree)==2:
        align_global[tree] = alignement_local.needlman_alignement(seqs[tree[0]],seqs[tree[1]])
    else:
        align_global[(tree)] = seqs[tree]

def traverse_tree(algnA, tree):
    """ The function that retrieves the name of the alignement that is close to algnA ( another alignment) from the tree """
    if isinstance(tree[0], tuple):
        for subtree in tree:
            if algnA == subtree:
                result = tree[1]  # Update the result when a match is found
            else:
                result = traverse_tree(algnA, subtree)
            if result is not None and algnA!=result:
                return result
def calculs_align_glob(algnA, algnB):
    """ The function that does the calculs of computing the alignment """
    """ The following instructions needs the names algnA and algnB as tuples so we make sure they are """
    if isinstance(algnA, tuple):
        pass
    else:
        algnA = tuple([algnA])
    if isinstance(algnB, tuple):
        pass
    else:
        algnB = tuple([algnB])
    """ computing the frequenceies of amino acids in each profil """
    freq_algnA = frequencies_calculator(algnA)
    freq_algnB = frequencies_calculator(algnB)
    """ Creating the matrix the has the scores and the cordinates of the cell that came from """
    dt = np.dtype([('float', np.float32), ('list', 'O')])
    matrix_glob_align = np.zeros((len(algnA[0])+1,len(algnB[0])+1),dtype=dt)
    """ The initialisation of the matrix """
    for i in range(1, len(algnA[0])+1):
        matrix_glob_align[i][0][0] = matrix_glob_align[i-1][0][0] + np.dot(constants.SCORE_GAP_vector, freq_algnA[:,i-1])
        matrix_glob_align[i][0][1] = (i-1,0)
    for i in range(1, len(algnB[0])+1):
        matrix_glob_align[0][i][0] = matrix_glob_align[0][i-1][0] + np.dot(constants.SCORE_GAP_vector, freq_algnB[:,i-1])
        matrix_glob_align[0][i][1] = (0,i-1)
    for i in range(1, len(algnA[0])+1):
        for j in range(1, len(algnB[0]) + 1):
            """ Computing the scores and getting the coordinates of the best score """
            score_1 = matrix_glob_align[i-1][j][0] + np.dot(constants.SCORE_GAP_vector,freq_algnA[:,i-1])
            score_2 = matrix_glob_align[i][j-1][0] + np.dot(constants.SCORE_GAP_vector, freq_algnB[:,j-1])
            score_3 = matrix_glob_align[i-1][j-1][0] + np.dot(freq_algnA[:,i-1],freq_algnB[:,j-1])
            if max(score_1 , score_2 , score_3) == score_1:
                matrix_glob_align[i][j][0] = score_1
                matrix_glob_align[i][j][1] = (i-1,j)
            if max(score_1 , score_2 , score_3) == score_2:
                matrix_glob_align[i][j][0] = score_2
                matrix_glob_align[i][j][1] = (i,j-1)
            if max(score_1 , score_2 , score_3) == score_3:
                matrix_glob_align[i][j][0] = score_3
                matrix_glob_align[i][j][1] = (i-1,j-1)
    """ Getting the cordinates of profils after finishing the calculs of the matrix """
    best_score_x = matrix_glob_align.shape[0]-1
    best_score_y = matrix_glob_align.shape[1]-1
    best_score_xy = (best_score_x, best_score_y)
    alignment_cordonnés = alignement_local.build_alignement_cordonnés(matrix_glob_align, best_score_xy, [])
    """ writing the amino acids according to the coordinates computed before """
    alignment_cordonnés = list(reversed(alignment_cordonnés))
    """ Building the alignement according to the coordinates """
    alignment_seqs = alignement_sequences_global(algnA, algnB, alignment_cordonnés)
    return alignment_seqs

def alignement_sequences_global(algnA, algnB, alignment_cordonnés):
    """ Uses the coordinates to generate the alignment between the sequences """
    if isinstance(algnA, tuple):
        pass
    else:
        algnA = tuple([algnA])
    if isinstance(algnB, tuple):
        pass
    else:
        algnB = tuple([algnB])
    longueur = len(alignment_cordonnés)
    nbr_seqs_A = len(algnA)
    nbr_seqs_B = len(algnB)
    algnA_aligned = []
    algnB_aligned = []
    for _ in range(nbr_seqs_A):
        algnA_aligned.append("")
    for _ in range(nbr_seqs_B):
        algnB_aligned.append("")
    for i in range(1,longueur):
        if alignment_cordonnés[i][0][0] == alignment_cordonnés[i-1][0][0]:
            for j in range(nbr_seqs_A):
                algnA_aligned[j] = algnA_aligned[j] + "-"
        if alignment_cordonnés[i][0][1] == alignment_cordonnés[i-1][0][1]:
            for j in range(nbr_seqs_B):
                algnB_aligned[j] = algnB_aligned[j] + "-"
        if alignment_cordonnés[i][0][0] != alignment_cordonnés[i-1][0][0]:
            for j in range(nbr_seqs_A):
                algnA_aligned[j] = algnA_aligned[j] + algnA[j][alignment_cordonnés[i-1][0][0]-1]
        if alignment_cordonnés[i][0][1] != alignment_cordonnés[i-1][0][1]:
            for j in range(nbr_seqs_B):
                algnB_aligned[j] = algnB_aligned[j] + algnB[j][alignment_cordonnés[i-1][0][1]-1]

    for i in range(nbr_seqs_B):
        algnA_aligned.append(algnB_aligned[i])
    return tuple(algnA_aligned)
def frequencies_calculator(algn):
    """ The function that computes the frequencies of the amino acids in the profils """
    nbr_seqs_algn = len(algn)
    freq_algn = np.zeros((21,len(algn[0])))
    for i in range(0, len(algn[0])):
        for aa in constants.AA_ARRAY:
            freq_aa = 0
            for j in range(0, nbr_seqs_algn):
                if aa == algn[j][i]:
                    freq_aa = freq_aa + 1/nbr_seqs_algn
            index_aa = np.where(constants.AA_ARRAY == aa)
            freq_algn[index_aa, i] = freq_aa
    return freq_algn
def matrix_conserved_positions(algn,sites):
    """ Computing the matrix that is going to be used by the algorithm neighbour joining """
    nbr_seq = len(algn)
    new_matrix = np.zeros((nbr_seq, nbr_seq))
    for i in range(nbr_seq):
        for j in range(i+1,nbr_seq):
            count = 0
            for k in sites:
                if algn[i][k] != algn[j][k]:
                    count = count + 1
            new_matrix[i,j] = count
    for i in range(nbr_seq):
        for j in range(i+1,nbr_seq):
            new_matrix[j,i] = new_matrix[i,j]
    return new_matrix

