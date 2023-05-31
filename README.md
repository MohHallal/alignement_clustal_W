# Bonjour
Ce projet à pour but de réaliser alignement multiple a la facon de ClustalW pour les séquences qui sont dans le fichier "opsines", du coup on doit:  
- realisant d’abord tout les alignement par pair avec l’algorithme de Needleman-Wunsch (fait).  
- A partir de la matrice de score ainsi obtenu construire un arbre guide des sequences en utilisant l’algorithme UPGMA (fait).  
- Realiser les alignements multiples successifs pour cr´eer l’alignement multiple de toutes les sequences. En suivant l’arbre guide. Pour cela on utilisera Needleman-Wunsch (sur les profils) pour aligner entre eux des alignements multiples.  (fait).
- On recuperera les positions conserv´ees (deux AA diff´erents `a une position donee) dans l’alignement multiple et recreera un arbre en utilisant encore une fois l’algorithme de Neighbor-Joining  (en cours).
  
(Pour exécuter ce qui est déja fait, executer "main").
