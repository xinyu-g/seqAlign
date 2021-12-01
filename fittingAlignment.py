UP = (-1,0)
LEFT = (0, -1)
TOPLEFT = (-1, -1)
ORIGIN = (0, 0)

def traceback_fitting(v, w, M, init_j, pointers):
    i, j = len(v), init_j
    new_v = []
    new_w = []
    while True:
        di, dj = pointers[i][j]
        if (di,dj) == LEFT:
            new_v.append('-')
            new_w.append(w[j-1])
        elif (di,dj) == UP:
            new_v.append(v[i-1])
            new_w.append('-')
        elif (di,dj) == TOPLEFT:
            new_v.append(v[i-1])
            new_w.append(w[j-1])
        i, j = i + di, j + dj
        if (i <= 0):
            break
    return ''.join(new_v[::-1]) + '\n'+''.join(new_w[::-1])

def fitting_align(short, reference, delta):
    """
    Returns the score of the maximum scoring alignment of short and all 
    substrings of reference. 
    
    :param: short the shorter of the two strings we are trying to align    
    :param: reference the longer string among whose substrings we are doing global alignment
    :param: delta the scoring function for the alphabet of the two strings
    
    :returns: a tuple (score, alignment)
    """
    M = [[0 for j in range(len(reference)+1)] for i in range(len(short)+1)]     
    pointers = [[ORIGIN for j in range(len(reference)+1)] for i in range(len(short)+1)] 
    score = float('-inf')
    init_j = 0
    # YOUR CODE HERE
    # raise NotImplementedError()
    

    m = len(short)
    n = len(reference)
    
    for i in range(1, m+1):
      M[i][0] = M[i-1][0] + delta[short[i-1]]['-']
      pointers[i][0] = UP
   
    for i in range(1, m+1):
      for j in range(1, n+1):
        M[i][j] = M[i-1][j] + delta[short[i-1]]['-']
        pointers[i][j] = UP
        if M[i][j-1] + delta['-'][reference[j-1]] > M[i][j]:
          M[i][j] = M[i][j-1] + delta['-'][reference[j-1]]
          pointers[i][j] = LEFT
        if M[i-1][j-1] + delta[short[i-1]][reference[j-1]] > M[i][j]:
          M[i][j] = M[i-1][j-1] + delta[short[i-1]][reference[j-1]]
          pointers[i][j] = TOPLEFT
    
    
    for j in range(0, n+1):
      if M[m][j] >= score:
        score = M[m][j]
        init_j = j

    alignment = traceback_fitting(short,reference,M, init_j,pointers)
    return score, alignment