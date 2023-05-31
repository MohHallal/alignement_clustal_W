import numpy as np
def arbre(matrice,sequences,size):
    if matrice.shape == (2,2):
        return sequences
    else:
        """ Detecting the closest sequences (which have the highest score from the matrix computed by Needleman Wunch """
        seq_1, seq_2 = closest_sequences(matrice, sequences)
        """ Retrieving their positions """
        pos_1 = sequences.index(seq_1)
        pos_2 = sequences.index(seq_2)
        """ Computing the matrix that is going to be used in the next round """
        updated_matrix = update_matrix(matrice, sequences.index(seq_1), sequences.index(seq_2),size, pos_1 , pos_2)
        """ updating the vector of the depths """
        size = update_size(size,sequences,seq_1,seq_2)
        """ Updating the tree by adding the closest sequences that we detected """
        updated_sequences = update_seq(list(sequences), seq_1, seq_2)
        return arbre(updated_matrix, updated_sequences,size)
def update_size(size, sequences, seq_1, seq_2):
    """ The function that updates the vector of the depths """
    position_seq_1 = sequences.index(seq_1)
    position_seq_2 = sequences.index(seq_2)
    new_element = size[position_seq_1] + size[position_seq_2]
    size.remove(size[position_seq_2])
    size.remove(size[position_seq_1])
    new_size = []
    for i in range(len(size)+1):
        if i == 0:
            new_size.append(new_element)
        else:
            new_size.append(size[i-1])
    return new_size
def closest_sequences(matrices, sequences):
    """ The function that detects the closest sequences """
    highest_score = matrices.max()
    vector = np.array(matrices)
    agglomeration = np.where(vector == highest_score)
    seq_1 = agglomeration[0][0]
    seq_2 = agglomeration[0][1]
    seq_1 = sequences[seq_1]
    seq_2 = sequences[seq_2]
    return seq_1, seq_2
def update_seq(sequences, seq_1, seq_2):
    """ Updating the tree """
    updated_sequences =[]
    sequences.remove(seq_1)
    sequences.remove(seq_2)
    seqs = [seq_1, seq_2]
    updated_sequences.append(tuple(seqs))
    for seq in sequences:
        updated_sequences.append(seq)
    updated_sequences = tuple(updated_sequences)
    return updated_sequences
def update_matrix(matrice, seq_1, seq_2, size, pos_1, pos_2):
    """ updating the matrix that is going to be used in the next round """
    dimension = matrice.shape
    temporary_matrix = np.delete(matrice, seq_1, 0)
    temporary_matrix = np.delete(temporary_matrix, seq_2 - 1, 0)
    temporary_matrix = np.delete(temporary_matrix, seq_1, 1)
    temporary_matrix = np.delete(temporary_matrix, seq_2 - 1, 1)
    new_matrix = np.zeros((dimension[0]-1, dimension[0]-1))
    dimension_2 = temporary_matrix.shape
    for i in range(0, dimension_2[0]):
        for j in range(0, dimension_2[0]):
            new_matrix[i+1,j+1] = temporary_matrix[i ,j]
    dimension_3 = new_matrix.shape
    profondeur = size[pos_1] + size[pos_2]
    j = 0
    for i in range(1, dimension_3[0]):
        while j in [seq_1, seq_2]:
            j = j+1
        new_matrix[0, i] = (matrice[j, seq_1] * size[pos_1] + matrice[j, seq_2] * size[pos_2])/ profondeur
        new_matrix[i, 0] = new_matrix[0 ,i]
        j = j+1
    return new_matrix
