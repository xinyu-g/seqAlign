UP = (-1,0)
LEFT = (0, -1)
TOPLEFT = (-1, -1)
ORIGIN = (0, 0)

def traceback_global(v, w, pointers):
    i,j = len(v), len(w)
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
        if (i <= 0 and j <= 0):
            break
    return ''.join(new_v[::-1])+'\n'+''.join(new_w[::-1])
    


def global_align(v, w, delta):
    """
    Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment as 
    computed by traceback_global. 
    
    :param: v
    :param: w
    :param: delta
    """
    M = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    pointers = [[ORIGIN for j in range(len(w)+1)] for i in range(len(v)+1)]
    score, alignment = None, None
    # YOUR CODE HERE
    m = len(v)
    n = len(w)
    
    for i in range(1, m+1):
      M[i][0] = M[i-1][0] + delta[v[i-1]]['-']
      pointers[i][0] = UP
    for j in range(1, n+1):
      M[0][j] = M[0][j-1] + delta['-'][w[j-1]]
      pointers[0][j] = LEFT
   
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
        

    score = M[m][n]
    alignment = traceback_global(v,w, pointers)
    return score, alignment