import numpy as np
SCORE_GAP = -1

SCORE_GAP_vector = np.zeros(20)
SCORE_GAP_vector = np.append(SCORE_GAP_vector, 1)

AA_ARRAY = np.array(['C','S','T','A','G','P','D','E','Q','N','H','R','K','M','I','L','V','W','Y','F','-'])
NAA = len(AA_ARRAY)
DAA = {letter:i for i,letter in enumerate(AA_ARRAY)}

BLOSUM62 = np.array([[9,-1,-1,0,-3,-3,-3,-4,-3,-3,-3,-3,-3,-1,-1,-1,-1,-2,-2,-2,-1],
					  [0,4,1,1,0,-1,0,0,0,1,-1,-1,0,-1,-2,-2,-2,-3,-2,-2,-1],
					  [0,0,5,0,-2,-1,-1,-1,-1,0,-2,-1,-1,-1,-1,-1,0,-2,-2,-2,-1],
					  [0,0,0,4,0,-1,-2,-1,-1,-2,-2,-1,-1,-1,-1,-1,0,-3,-2,-2,-1],
					  [0,0,0,0,6,-2,-1,-2,-2,0,-2,-2,-2,-3,-4,-4,-3,-2,-3,-3,-1],
					  [0,0,0,0,0,7,-1,-1,-1,-2,-2,-2,-1,-2,-3,-3,-2,-4,-3,-4,-1],
					  [0,0,0,0,0,0,6,2,0,1,-1,-2,-1,-3,-3,-4,-3,-4,-3,-3,-1],
					  [0,0,0,0,0,0,0,5,2,0,0,0,1,-2,-3,-3,-2,-3,-2,-3,-1],
					  [0,0,0,0,0,0,0,0,5,0,0,1,1,0,-3,-2,-2,-2,-1,-3,-1],
					  [0,0,0,0,0,0,0,0,0,6,1,0,0,-2,-3,-3,-3,-4,-2,-3,-1],
					  [0,0,0,0,0,0,0,0,0,0,8,0,-1,-2,-3,-3,-3,-2,2,-1,-1],
					  [0,0,0,0,0,0,0,0,0,0,0,5,2,-1,-3,-2,-3,-3,-2,-3,-1],
					  [0,0,0,0,0,0,0,0,0,0,0,0,5,-1,-3,-2,-2,-3,-2,-3,-1],
					  [0,0,0,0,0,0,0,0,0,0,0,0,0,5,1,2,1,-1,-1,0,-1],
					  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,2,3,-3,-1,0,-1],
					  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,1,-2,-1,0,-1],
					  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,-3,-1,-1,-1],
					  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,2,1,-1],
					  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,3,-1],
					  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,-1],
					  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]])

for i in range(NAA-1):
	for j in range(i+1,NAA):
		BLOSUM62[j,i] = BLOSUM62[i,j]
