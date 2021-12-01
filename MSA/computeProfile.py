
UP = (-1, 0)
LEFT = (0, -1)
TOPLEFT = (-1, -1)

def compute_profile(alignment, alphabet):
    """
    Given an alphabet an a multiple sequence alignment in that alphabet,
    computes and returns its profile representation
    
    :param: alignment is a list of lists of characters in the alphabet
    :param: alphabet is a list of characters in the alphabet from which the strings are
            constructed
    :return: a dictionary where dict[x][i] is the frequency of the character
             x in the i-th position of the alignment.
    """
    
    if not alignment:
        return {}

    n = len(alignment)
    l = len(alignment[0])

    profile = {}
    
    for al in alphabet:
      profile[al] = [0] * l

    for i in range(n):
      for j in range(l):
        c = alignment[i][j]

        profile[c][j] += 1 


    for k,v in profile.items():
      profile[k] = [num/float(n) for num in v]
    
    return profile


def compute_tau(profile, alphabet, delta):
    """
    Given a profile, an alphabet and a scoring function for that alphabet,
    returns the scoring function for aligning a character in the alphabet
    to a column in the profile
    
    :param: profile is the profile representation of the multiple sequence
            we are aligning against
    :param: alphabet is the alphabet of characters that compose our sequences
    :param: delta is the scoring function between characters in our alphabet
    
    :return: The scoring function tau such that tau[x][i] is the score for aligning
             character x with column i of the profile.
    """
    
    tau = {}

    for al in alphabet:
      for k,v in profile.items():
        for j in range(len(v)):
          if not tau.get(al):
            tau[al] = [0] * len(v)
          tau[al][j] += v[j] * delta[al][k]


    return tau



def traceback(aln1, aln2, pointers):
    i = len(aln1[0])-1
    j = len(aln2[0])-1
    new_al1 = [list(v) for v in aln1]
    new_al2 = [list(w) for w in aln2]
    while True:
        di, dj = pointers[i][j]
        if (di, dj) == LEFT:
            for seq1 in new_al1:
                seq1.insert(i, '-')
        if (di, dj) == UP:
            for seq2 in new_al2:
                seq2.insert(j, '-')
        i, j = i + di, j + dj
        if (i <= 0 and j <= 0):
            break
    new_alignment = []
    for seq in new_al1:
        new_alignment.append(''.join(seq))
    for seq in new_al2:
        new_alignment.append(''.join(seq))
    return new_alignment


def align_sequence_profile(alignment, sequence, alphabet, delta):
    """
    This function aligns a sequence against a multiple sequence alignment
    
    :param: alignment is the multiple sequence alignment are aligning against.
            This is a list of list of characters
    :param: sequence is the new sequence we are aligning to the multiple alignment.
            This is a list of characters
    :param: alphabet is a list of characters that could compose the sequences in
            the alignments.
    :param: delta is the scoring function for aligning characters in our alphabet.
            delta[x][y] is the score for aligning the characters x and y.
    
    
    :return: a list of lists of characters in the alphabet, representing the 
             new multiple sequence alignment
    """
    # Base case when there is an empty multiple alignment
    if not alignment:
        return [sequence]
    M = [[0 for _ in range(len(alignment[0]))] for _ in range(len(sequence))] 
    pointers = [[(0,0) for _ in range(len(alignment[0]))] for _ in range(len(sequence))]
    score = None
    
    profile = compute_profile(alignment, alphabet)
    tau = compute_tau(profile, alphabet, delta)

    for i in range(len(sequence)):
        for j in range(len(alignment[0])):
            if i == 0 and j == 0:   
                M[i][j] = 0
            elif i == 0:
                M[i][j] = M[i][j-1] + tau['-'][j-1]
                pointers[i][j] = LEFT
            elif j == 0:
                sequence[i-1]
                M[i][j] = M[i-1][j] + delta[sequence[i-1]]['-']
                pointers[i][j] = UP
            else:
                best_sub = max([(LEFT, M[i][j-1] + tau['-'][j-1]), 
                               (UP, M[i-1][j] + delta[sequence[i-1]]['-']), 
                               (TOPLEFT, M[i-1][j-1] + tau[sequence[i-1]][j-1])], key = lambda x: x[1])
                pointers[i][j] = best_sub[0]
                M[i][j] = best_sub[1]

    score = M[-1][-1]
    return score, traceback([sequence], alignment, pointers)


def compute_sigma(p,q,alphabet, delta):
    """
    :param: p is the profile for the first multiple alignment
    :param: q is the profile for the second multiple alignment
    :param: alphabet is the list of all characters in our sequences
    :param: delta is the scoring function for aligning characters in our alphabet
    
    :returns: a list of lists sigma such that sigma[i][j] is the score for aligning column
              i of p with column j of q
    """
    
    size_p = len(p[alphabet[0]])
    size_q = len(q[alphabet[0]])
    sigma = [[0 for _ in range(size_q)] for _ in range(size_p)]
     
    for i in range(size_p):
      for j in range(size_q):
        for x in alphabet:
          for y in alphabet:
            sigma[i][j] += p[x][i]*q[y][j]*delta[x][y]
    return sigma

def compute_sigma(p,q,alphabet, delta):
    """
    :param: p is the profile for the first multiple alignment
    :param: q is the profile for the second multiple alignment
    :param: alphabet is the list of all characters in our sequences
    :param: delta is the scoring function for aligning characters in our alphabet
    
    :returns: a list of lists sigma such that sigma[i][j] is the score for aligning column
              i of p with column j of q
    """
    
    size_p = len(p[alphabet[0]])
    size_q = len(q[alphabet[0]])
    sigma = [[0 for _ in range(size_q)] for _ in range(size_p)]

     
    for i in range(size_p):
      for j in range(size_q):
        for x in alphabet:
          for y in alphabet:
            sigma[i][j] += p[x][i]*q[y][j]*delta[x][y]

    return sigma

def greedy_progressive_align(alignments, alphabet, delta):
    """
    :param: alignments is a list of list of strings representing the sequences to be aligned
            Note: This is because we need to represent our single sequences as multiple alignments
            ,and multiple alignments are lists of strings
    :param: alphabet is the alphabet from which the sequences are derived
    :param: delta is a scoring function. delta(x,y) gives us the score for aligning 
            character x with character y in our alphabet
            
    :returns: the greedy optimal multiple sequence alignment for a given set of sequences, and the score for that alignment
    """
    
    align_score, greedy_align = None, None
    while True: 

        if len(alignments) == 1:
          break
        

        best_score = -float("inf")
        best_alignment = None
        best_m = -1
        best_n = -1

        # Compute pairwise distances 
        for m in range(len(alignments)):
            for n in range(m):
                score, align = align_profile_profile(alignments[m], alignments[n], alphabet, delta)
                if not best_score or score > best_score:
                  best_score, best_alignment = score, align
                  best_m, best_n = m, n

        next_alignments = [best_alignment]
        for i in range(len(alignments)):
            if i!=best_m and i!=best_n:
                next_alignments.append(alignments[i])
        alignments = next_alignments
        align_score, greedy_align = best_score, best_alignment

    return align_score, greedy_align