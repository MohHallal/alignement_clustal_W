import numpy as np
import upgma
def nj_main(matrix_nj, species):
    if matrix_nj.shape == (2,2):
        species_str = str(species)
        species_str = species_str.replace("\\", "")
        return species_str
    else:
        nbr_seq = len(species)
        species = list(species)
        """ adding the additional column thah has the sum of ech row """
        matrix_nj_with_r = matrix_nj_r(matrix_nj, nbr_seq)
        """ Computing the matrix that is going to be used to identify the closest species """
        matrix_distances = nj_distances(matrix_nj_with_r,nbr_seq)
        """ Identifying the closest species """
        spec_1_index, spec_2_index = closest_sequences(matrix_distances)
        """ Computing the length of the branches """
        dist_spec_1_U = 0.5 * matrix_nj_with_r[spec_1_index, spec_2_index] + 0.5 * (matrix_nj_with_r[spec_1_index, -1] - matrix_nj_with_r[spec_2_index, -1])
        dist_spec_2_U = 0.5 * matrix_nj_with_r[spec_1_index, spec_2_index] + 0.5 * (matrix_nj_with_r[spec_2_index, -1] - matrix_nj_with_r[spec_1_index, -1])
        """ Updating the tree """
        name_with_dist_1 = str(species[spec_1_index]) + ": " + str(dist_spec_1_U)
        name_with_dist_2 = str(species[spec_2_index]) + ": " + str(dist_spec_2_U)
        species[spec_1_index] = name_with_dist_1
        species[spec_2_index] = name_with_dist_2
        """ Preparing for the newt round, a new matrix and a new tree """
        new_species = upgma.update_seq(species, name_with_dist_1, name_with_dist_2)
        new_matrix_nj = update_matrix_nj(matrix_nj, spec_1_index, spec_2_index)
        return nj_main(new_matrix_nj , new_species)

def update_matrix_nj(matrix_nj, spec_1_index, spec_2_index):
    """ Updating the matrix after finding the closest sequences """
    dimension = matrix_nj.shape
    temporary_matrix = np.delete(matrix_nj, spec_1_index, 0)
    temporary_matrix = np.delete(temporary_matrix, spec_2_index - 1, 0)
    temporary_matrix = np.delete(temporary_matrix, spec_1_index, 1)
    temporary_matrix = np.delete(temporary_matrix, spec_2_index - 1, 1)
    new_matrix = np.zeros((dimension[0]-1, dimension[0]-1))
    dimension_2 = temporary_matrix.shape
    for i in range(0, dimension_2[0]):
        for j in range(0, dimension_2[0]):
            new_matrix[i+1,j+1] = temporary_matrix[i ,j]
    dimension_3 = new_matrix.shape
    j = 0
    for i in range(1, dimension_3[0]):
        while j in [spec_1_index, spec_2_index]:
            j = j+1
        new_matrix[0, i] = (matrix_nj[spec_1_index,j] + matrix_nj[spec_2_index,j] - matrix_nj[spec_1_index, spec_2_index])/2
        new_matrix[i, 0] = new_matrix[0, i]
        j = j+1
    return new_matrix
def matrix_nj_r(matrix_nj, nbr_seq):
    """ Computing the additional column that has the sum of the each line that is used by Neighbour Joining """
    r = np.sum(matrix_nj, axis=1)
    r = r/(nbr_seq-2)
    r = r.reshape((-1, 1))
    matrix_nj_with_r = np.hstack((matrix_nj, r))
    return matrix_nj_with_r
def nj_distances(matrix_nj_with_r, nbr_seq):
    """ Computing the matrix that is used to identify the closest species """
    np.set_printoptions(suppress=True, formatter={'float_kind': '{:.1f}'.format})
    new_matrix = np.zeros((nbr_seq, nbr_seq))
    for i in range(nbr_seq):
        for j in range(i+1,nbr_seq):
            dist_i_j = matrix_nj_with_r[i,j] - matrix_nj_with_r[i,-1] - matrix_nj_with_r[j,-1]
            new_matrix[i,j] = dist_i_j
    for i in range(nbr_seq-1):
        for j in range(i+1,nbr_seq):
            new_matrix[j, i] = new_matrix[i,j]
    return new_matrix
def closest_sequences(matrix_dist):
    """ The function that identifies the minimum of  the matrix and retrieves the indexes of the species corresponding"""
    matrix_without_zero = np.where(matrix_dist == 0, np.inf, matrix_dist)
    agglomeration = np.unravel_index(np.argmin(matrix_without_zero), matrix_without_zero.shape)
    spec_1_index = agglomeration[0]
    spec_2_index = agglomeration[1]
    return spec_1_index, spec_2_index