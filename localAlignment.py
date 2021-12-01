UP = (-1,0)
LEFT = (0, -1)
TOPLEFT = (-1, -1)
ORIGIN = (0, 0)

def traceback_local(v, w, M, init_i, init_j, pointers):
    i,j = init_i, init_j
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
        if (M[i][j] == 0):
            break
    return ''.join(new_v[::-1]) + '\n'+''.join(new_w[::-1])

def local_align(v, w, delta):
    """
    Returns the score of the maximum scoring alignment of all possible substrings of v and w. 
    
    :param: v
    :param: w
    :param: delta
    """
    M = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]        
    pointers = [[ORIGIN for j in range(len(w)+1)] for i in range(len(v)+1)] 
    score = 0
    init_i, init_j = 0,0
    # YOUR CODE HERE

    m = len(v)
    n = len(w)

    for i in range(1, m+1):
      for j in range(1, n+1):
        M[i][j] = M[i-1][j] + delta[v[i-1]]['-']
        pointers[i][j] = UP
        if M[i][j-1] + delta['-'][w[j-1]] > M[i][j]:
          M[i][j] = M[i][j-1] + delta['-'][w[j-1]]
          pointers[i][j] = LEFT
        if M[i-1][j-1] + delta[v[i-1]][w[j-1]] > M[i][j]:
          M[i][j] = M[i-1][j-1] + delta[v[i-1]][w[j-1]]
          pointers[i][j] = TOPLEFT
        if 0 > M[i][j]:
          M[i][j] = 0
          pointers[i][j] = ORIGIN

    for i in range(0, m+1):
      for j in range(0, n+1):
        if M[i][j] >= score:
          score = M[i][j]
          init_i, init_j = i, j

    # raise NotImplementedError()
    alignment = traceback_local(v, w, M, init_i, init_j , pointers)
    return score, alignment 