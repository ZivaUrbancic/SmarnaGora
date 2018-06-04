# gregorjerse/rt2. (2018). GitHub. Retrieved 1 June 2018, from https://github.com/gregorjerse/rt2/blob/master/2015_2016/lab13/Extending%20values%20on%20vertices.ipynb

from itertools import combinations, chain

def simplex_closure(a):    
    """Returns the generator that iterating over all subsimplices (of all dimensions) in the closure
    of the simplex a. The simplex a is also included.
    """
    return chain.from_iterable([combinations(a, l) for l in range(1, len(a) + 1)])
        
def closure(K):
    """Add all missing subsimplices to K in order to make it a simplicial complex."""
    return list({s for a in K for s in simplex_closure(a)})

def contained(a, b):
    """Returns True is a is a subsimplex of b, False otherwise."""
    return all((v in b for v in a))

def star(s, cx):
    """Return the set of all simplices in the cx that contais simplex s.
    """
    return {p for p in cx if contained(s, p)}

def intersection(s1, s2):
    """Return the intersection of s1 and s2."""
    return list(set(s1).intersection(s2))

def link(s, cx):
    """Returns link of the simplex s in the complex cx.
    """
    # Link consists of all simplices from the closed star that have 
    # empty intersection with s.
    return [c for c in closure(star(s, cx)) if not intersection(s, c)]

def simplex_value(s, f, aggregate):
    """Return the value of f on vertices of s
    aggregated by the aggregate function.
    """
    return aggregate([f[v] for v in s])

def lower_link(s, cx, f):
    """Return the lower link of the simplex s in the complex cx.
    The dictionary f is the mapping from vertices (integers)
    to the values on vertices.
    """
    sval = simplex_value(s, f, min)
    return [s for s in link(s, cx) 
            if simplex_value(s, f, max) < sval]

def join(a, b):
    """Return the join of 2 simplices a and b."""
    return tuple(sorted(set(a).union(b)))

def extend(K, f):
    """Extend the field to the complex K.
    Function on vertices is given in f.
    Returns the pair V, C, where V is the dictionary containing discrete gradient vector field
    and C is the list of all critical cells.
    """
    V = dict()
    C = []
    for v in (s for s in K if len(s)==1):
        ll = lower_link(v, K, f) 
        if len(ll) == 0:
            C.append(v)
        else:
            V1, C1 = extend(ll, f)
            mv, mc = min([(f[c[0]], c) for c in C1 if len(c)==1])
            V[v] = join(v, mc)
            for c in (c for c in C1 if c != mc):
                C.append(join(v, c))
            for a, b in V1.items():
                V[join(a, v)] = join(b, v)
    return V, C

if __name__ == "__main__":
    K = closure([(1, 2, 3), (2, 3, 4)])
    f = {1: 0, 2: 2, 3: 2, 4: 0}
    print(extend(K, f))