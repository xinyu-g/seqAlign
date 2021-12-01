def check_tree(D,T):
    """
    Implement floyd-warshall algorithm on the graph defined in T
    
    """
    
    V = len(T)
    dT = [[float("inf") for i in range(V)] for j in range(V)]
    n = (V + 2)//2
    
    # Length check
    if n!=len(D):
        return False
    
    # fill in existing edges
    for k in T:
        for l in T[k]:
            dT[k][l] = T[k][l]
            dT[l][k] = T[l][k]
    
    # fill in the diagonal elements 
    for i in range(len(dT)):
        dT[i][i] = 0.0
    
    # relax edges 
    for k in range(V):
        for i in range(V):
            for j in range(V):
                if dT[i][j] > dT[i][k] + dT[k][j]:
                    dT[i][j] = dT[i][k] + dT[k][j]
                    dT[j][i] = dT[i][j]
    
    for i in range(n) :
        print(dT[i][:n])
    # Check each value in dT
    for i in range(n):
        for j in range(n):
            if D[i][j] != dT[i][j]:
                return False
    return True

def min_S_value(D, u):
    """
    returns the value (i,j) for which 
    (m-2)*D[i][j] - u_i - u_j is minimum
    """
    m = len(D)
    min_S, min_i, min_j = float("inf"),-1,-1
    for k in D:
        for l in D[k]:
            if l!=k:
                crit = (m-2)*D[k][l] - u[k] - u[l]
                if crit < min_S:
                    min_S = crit
                    min_i = k
                    min_j = l
    return (min_i, min_j)


def neighbor_join(D):
    """
    Takes a distance matrix D, and returns the tree T 
    consistent with the closest additive matrix D' to D.
    
    :param: D is a dict of dicts representing pairwise distances between leaves
    :return: a dict of dicts that contains all the edges with their weights in the tree defined by D'.    
    """
    # clusters = {}
    T = {}
    r = len(D)
    while len(D)>2:
      u = {}
      for k,v in D.items():
        u[k] = sum(v.values())
      i,j = min_S_value(D,u)
      print(i, j)
      T[r] = {}
      if not T.get(i):
        T[i] = {}
      if not T.get(j):
        T[j] = {}

      T[r][i] = 0.5*(D[i][j] + (u[i] - u[j])/(len(D)-2))
      T[r][j] = 0.5*(D[i][j] + (u[j] - u[i])/(len(D)-2))
      T[i][r] = 0.5*(D[i][j] + (u[i] - u[j])/(len(D)-2))
      T[j][r] = 0.5*(D[i][j] + (u[j] - u[i])/(len(D)-2))
         
      D[r] = {}
      # print(D)
      for m, _ in D.items():
        if m != i and m != j and m != r: 
          # print(m)
          
          D[m][r] = 0.5*(D[i][m] + D[j][m] - D[i][j])
          D[r][m] = D[m][r]
          
      del D[i]
      del D[j]
      for k in D.keys():
        if D[k].get(i):
          del D[k][i]
        if D[k].get(j):
          del D[k][j]
      r = r + 1

    k1, k2 = list(D.keys())[0], list(D.keys())[1]
    if not T.get(k1):
      T[k1] = {}
    if not T.get(k2):
      T[k2] = {}

    T[k1][k2] = D[k1][k2]
    T[k2][k1] = D[k2][k1]

    return T