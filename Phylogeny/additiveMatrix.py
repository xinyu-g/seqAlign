def is_additive(D):
    """
    Returns true if the square matrix D is additive
    
    :param: D is an nxn list of lists of ints 
    :return: true if D is an additive distance matrix
    """
    n = len(D)
    for i in range(n):
      for j in range(i, n):
        for k in range(j, n):
          for l in range(k, n):
            s1 = D[i][j] + D[k][l]
            s2 = D[i][k] + D[j][l]
            s3 = D[i][l] + D[j][k]
            lst = sorted([s1, s2, s3])
            if lst[1] != lst[2]:
              return False 
    return True