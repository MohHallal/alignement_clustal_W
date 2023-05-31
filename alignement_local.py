import numpy as np
import constants as cst
def needlman_score(sequences):
    """ Computing the scores of Needelman Wunch in order to establish the guide tree that is going to be used to compute the UPGMA """
    seq_names = list(sequences.keys())
    scores = np.zeros((len(seq_names),len(seq_names)))
    for i in range(len(seq_names) - 1):
        seq_1 = sequences[seq_names[i]]
        for j in range(i + 1, len(seq_names)):
            seq_2 = sequences[seq_names[j]]
            scores[i,j] = calculs(seq_1,seq_2)
    for i in range(len(seq_names) - 1):
        for j in range(i + 1, len(seq_names)):
            scores[j,i] = scores[i,j]
    return scores
def calculs(seq_1, seq_2):
    """ The function that computes the values of the matrix thet is used by UPGMA """
    seq_1 = "-" + seq_1
    seq_2 = "-" + seq_2
    score_align = np.zeros((len(seq_1), len(seq_2)))
    for i in range(1,len(seq_1)):
        score_align[i,0] = score_align[i-1,0]+cst.SCORE_GAP
    for i in range(1,len(seq_2)):
        score_align[0,i] = score_align[0,i-1]+cst.SCORE_GAP
    for i in range(1,len(seq_1)):
        for j in range(1,len(seq_2)):
            aa1 = seq_1[i]
            aa2 = seq_2[j]
            aa1 = cst.DAA[aa1]
            aa2 = cst.DAA[aa2]
            score_1 = score_align[i-1 , j] + cst.SCORE_GAP
            score_2 = score_align[i , j-1] + cst.SCORE_GAP
            score_3 = score_align[i-1 , j-1] + cst.BLOSUM62[aa1,aa2]
            score_align[i,j] = max(score_1 , score_2 , score_3)
    return score_align[len(seq_1)-1,len(seq_2)-1]

def needlman_alignement(seq_1, seq_2):
    """ Computing the alignment of two sequences """
    seq_1 = "-" + seq_1
    seq_2 = "-" + seq_2
    """ Creating the matrix that is going to have the scores and the coordinates of the cell that has been used to compute the score """
    dt = np.dtype([('int', np.int32), ('list', 'O')])
    score_align = np.zeros((len(seq_1), len(seq_2)),dtype=dt)
    """ Initialisation of the matrix """
    for i in range(1,len(seq_1)):
        score_align[i][0][0] = score_align[i-1][0][0]+cst.SCORE_GAP
        score_align[i][0][1] = (i-1,0)
    for i in range(1,len(seq_2)):
        score_align[0][i][0] = score_align[0][i-1][0]+cst.SCORE_GAP
        score_align[0][i][1] = (0,i-1)
    """ Completing the rest of the matrix """
    for i in range(1,len(seq_1)):
        for j in range(1,len(seq_2)):
            aa1 = seq_1[i]
            aa2 = seq_2[j]
            aa1 = cst.DAA[aa1]
            aa2 = cst.DAA[aa2]
            """ Computing the three scores """
            score_1 = score_align[i-1][j][0] + cst.SCORE_GAP
            score_2 = score_align[i][j-1][0] + cst.SCORE_GAP
            score_3 = score_align[i-1][j-1][0] + cst.BLOSUM62[aa1,aa2]
            """ Identifying the best score and identifying the coordinates from where it comes """
            if max(score_1 , score_2 , score_3) == score_1:
                score_align[i][j][0] = score_1
                score_align[i][j][1] = (i-1,j)
            elif max(score_1 , score_2 , score_3) == score_2:
                score_align[i][j][0] = score_2
                score_align[i][j][1] = (i,j-1)
            else:
                score_align[i][j][0] = score_3
                score_align[i][j][1] = (i-1,j-1)
    """ Retrieving all the cordinates used to compute the matrix """
    best_score_x = score_align.shape[0]-1
    best_score_y = score_align.shape[1]-1
    best_score_xy = (best_score_x, best_score_y)
    alignment_cordonnés = build_alignement_cordonnés(score_align, best_score_xy,[])
    alignment_cordonnés = list(reversed(alignment_cordonnés))
    """ Using the cordinates to compute the represntation of the alignement between the two sequences """
    alignement_seqs = alignement_sequences(seq_1,seq_2,alignment_cordonnés)
    return alignement_seqs[0], alignement_seqs[1]
def alignement_sequences(seq_1,seq_2,alignment_cordonnés):
    """ Uses the coordinates to generate the alignment between the two sequences """
    seq_1_aligned = ""
    seq_2_aligned = ""
    longueur = len(alignment_cordonnés)
    for i in range(1,longueur):
        if alignment_cordonnés[i][0][0] == alignment_cordonnés[i-1][0][0]:
            seq_1_aligned = seq_1_aligned + "-"
        if alignment_cordonnés[i][0][1] == alignment_cordonnés[i-1][0][1]:
            seq_2_aligned = seq_2_aligned + "-"
        if alignment_cordonnés[i][0][0] != alignment_cordonnés[i-1][0][0]:
            seq_1_aligned = seq_1_aligned + seq_1[alignment_cordonnés[i-1][0][0]]
        if alignment_cordonnés[i][0][1] != alignment_cordonnés[i-1][0][1]:
            seq_2_aligned = seq_2_aligned + seq_2[alignment_cordonnés[i-1][0][1]]
    seq_1_aligned = seq_1_aligned + seq_1[-1]
    seq_2_aligned = seq_2_aligned + seq_2[-1]
    return seq_1_aligned,seq_2_aligned
def build_alignement_cordonnés(score_align, best_score_xy,coordonnés):
    """ The function the retrieves the coordinates """
    coordonnés.append([best_score_xy])
    if score_align[best_score_xy[0],best_score_xy[1]][1] != (0,0):
        return build_alignement_cordonnés(score_align,score_align[best_score_xy[0],best_score_xy[1]][1],coordonnés)
    else:
        return coordonnés
