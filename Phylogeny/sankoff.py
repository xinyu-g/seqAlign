class Node:
    def __init__(self, name,  character, states_alphabet, root=False, leaf=False):
        self.name = name
        self.state = character
        self.leaf = leaf
        self.root = root
        self.parent = None
        self.children = []
        self.mu = {s:float("inf") for s in states_alphabet}
        self.min_val = None
        
    def add_child(self, child):
        self.children.append(child)
        child.assign_parent(self)
        
    def assign_parent(self, parent):
        self.parent = parent

        
def cost(s,t):
    """
    Computes the cost of a mutation between s and t
    """
    return 0 if s==t else 1


def sankoff_fill(root_node, states_alphabet):
    """
    Takes a reference to the root of a character based phylogeny tree,
    and fills in the mu values for all nodes
    
    """
    # states_alphabet = sorted(states_alphabet)
    if root_node.leaf:

        root_node.mu[root_node.state] = 0
    else:
        # recurse (post order traversal)
        for w in root_node.children:
            sankoff_fill(w, states_alphabet)


        
        # post order processing of current node
        # min_cost = float("inf")
        for s in states_alphabet:
          prev = []
          for w in root_node.children:
            costs = []
            for t in states_alphabet:
              c = cost(s, t) + w.mu[t]
              costs.append(c)
            prev.append(min(costs))
          root_node.mu[s] = sum(prev)


    root_node.min_val = min(root_node.mu.values())
    
def sankoff_backtrace(root_node, state_alphabet):
    """
    Takes a reference to the root node of a character based phylogeny tree and 
    fills in the state values with the appropriate character. 
    
    """
    # base case
    if root_node.leaf:
        return
    
    # process current node
    if root_node.root:
      min_cost = float('inf')
      for s in state_alphabet:
        if root_node.mu[s] < min_cost:
            min_cost = root_node.mu[s]
            root_node.state = s
      # for w in root_node.children:
      #   w.parent = root_node.state

    else:
      min_cost = float('inf')
      a = root_node.parent.state
      for b in state_alphabet:
        c = cost(a,b) + root_node.mu[b]
        if c < min_cost:
          min_cost = c
          root_node.state = b

    
    # recurse
    for w in root_node.children:
        sankoff_backtrace(w, state_alphabet)   