import alignement_global
import alignement_local
import neibour_joining
import upgma
from collections import Counter
def seperating_sequences(file):
    """ Creates a dictionnary from the FASTA file, the species as keys and the protein sequence as its item """
    with open(file, "r", encoding="utf-8") as fasta_file:
        seq = ""
        seq_name = ""
        dico_seqs = {}
        for line in fasta_file:
            if line[0] == ">":
                if seq:
                    dico_seqs[seq_name] = seq
                    seq = ""
                seq_name = line[1:].strip()
            else:
                seq = seq + line.strip()
        if seq:
            dico_seqs[seq_name] = seq
        return dico_seqs
def get_names_from_tree(tree):
    """ It take the guide tree as argument and display all its components from right to left """
    if isinstance(tree, tuple):
        return sum((get_names_from_tree(species) for species in tree), ())
    else:
        return (str(tree),)
def consrv_pos(algn):
    """ It takes the multiple alignment as argument and guves back a tuple that have all the conserved positions """
    longueur = len(algn[0])
    nbr_seq = len(algn)
    positions = []
    for i in range(longueur):
        column = []
        for j in range(nbr_seq):
            column.append(algn[j][i])
        counter = Counter(column)
        if len(counter.keys()) == 2:
            positions.append(i)
    return tuple(positions)


if __name__=="__main__":
    """ Reading the FASTA file """
    input_file = "opsines"
    sequences = seperating_sequences(input_file)
    """ Computing the scores of needlemen wunch for the sequences """
    scores = alignement_local.needlman_score(sequences)
    print("Matrice des scores de NW, qui servira pour construire l'arbre guide:")
    print(scores)
    """ Computing the guide tree """
    sequences_names = list(sequences.keys())
    size = [1] * len(sequences_names)
    tree = upgma.arbre(scores,sequences_names,size)
    print("L'arbre guide:")
    print(tree)
    """ Computing the multiple alignment """
    global_alignment = alignement_global.align_glob(tree ,sequences)
    print("L'alignement des séquences est affiché selon cette ordre des espéces: ")
    ordered_species = get_names_from_tree(tree)
    print(ordered_species)
    for i in range(len(sequences_names)):
        print(global_alignment[i])
    """ Detecting the conserved positions in the sequences """
    conserved_positions = consrv_pos(global_alignment)
    print("Les postions des colonnes conservées:")
    print(conserved_positions)
    """ Computing the phylogenetic tree by the Neighbour Joining algorithm """
    matrix_NJ = alignement_global.matrix_conserved_positions(global_alignment, conserved_positions)
    print("La matrice selon les sites conservés avec le méme ordre des espéces:")
    print(matrix_NJ)
    nj_tree = neibour_joining.nj_main(matrix_NJ, ordered_species)
    print("L'arbre de l'algorithme neighbour joining: ")
    print(nj_tree)
