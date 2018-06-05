from random import choice
from collections import defaultdict
from numpy import infty
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


def extract_raw(K, fun, p):
    """Extend the field to the complex K.
    Function on vertices is given in fun.
    Returns the pair V, C, where V is the dictionary containing discrete gradient vector field
    and C is the list of all critical cells.
    """
    V = dict()
    C = []
    for v in (s for s in K if len(s) == 1):
        ll = lower_link(v, K, fun)
        if len(ll) == 0:
            C.append(v)
        else:
            V1, C1 = extract(ll, fun, p)
            mv, mc = min([(fun[c[0]], c) for c in C1 if len(c) == 1])
            V[v] = join(v, mc)
            for c in (c for c in C1 if c != mc):
                C.append(join(v, c))
            for a, b in V1.items():
                V[join(a, v)] = join(b, v)
    return V, C


def generate_all_sxs(K):
    """Returns a set of all simplices in simplicial complex K."""
    sxi = set()
    for sx in K:
        if type(sx[0]) != tuple:
            sx = [tuple(i) for i in sx]
        n = len(sx)
        for i in range(n):
            sxi = sxi.union(combinations(sx, i + 1))
    return sxi


def boundary(a):
    """Return the boundary (codimension one faces) of the simplex a.
    """
    return combinations(a, len(a)-1)


def expand_path(V, ends, path):
    """Expands a path with the given prefix which is stored in the
    list path. The last simplex in the list path in the end of the
    prefix. The path stops at the simplices in the set/list ends.
    Returns the list of all paths with the given prefix which end in
    critical cells.
    """
    if path[-1] in ends:      # path already ends at the critical simplex
        yield path
    elif path[-1] in V:       # can continue the path
        children = (s for s in boundary(V[path[-1]]) if s != path[-1])
        paths = (expand_path(V, ends, path + [V[path[-1]], c]) for c in children)
        for p in chain.from_iterable(paths):
            yield p


def find_paths_between_critical_simplices(V, C):
    """Find all paths connecting pairs of critical simplices a and b.
    (starting in the boundary of a and ending in b).
    Returns dictionary whose keys are pairs of connected critical simplices
    and value for the given key is a list of paths connecting those simplices.
    """
    paths = defaultdict(list)
    for candidate in (s for s in C if len(s) > 1):
        for start in boundary(candidate):
            for path in expand_path(V, C, [start]):
                paths[(candidate, path[-1])].append(path)
    return paths


def cancel_along_path(V, C, start, path):
    """Cancel two critical cells connected with the given path."""
    reversed_path = list(reversed(path)) + [start]
    for i in range(0, len(reversed_path), 2):
        V[reversed_path[i]] = reversed_path[i + 1]
    C.remove(start)
    C.remove(path[-1])


def cancel_all(V, C):
    """Cancel all possible critical pairs."""
    while True:
        paths = find_paths_between_critical_simplices(V, C)
        candidates = [(pair, path[0]) for pair, path in paths.items() if len(path) == 1]
        if len(candidates) == 0:
            break
        pair, path = choice(candidates)
        cancel_along_path(V, C, pair[0], path)


def max_f(sx, fun):
    """Returns the highest value of fun in 0-dimensional subsimplices of sx."""
    return max([fun[s] for s in list(sx)])


def extract_cancel(K, fun, p, j, V, C):
    for sigma in [c for c in C if len(c) == j]:
        paths = find_paths_between_critical_simplices(V, [c for c in C if len(c) in [j - 1, j]])
        candidates = {}
        ms = {}
        for pair in paths:
            starting, ending = pair
            if starting == sigma and len(paths[pair]) == 1 and max_f(ending, fun) > max_f(starting, fun) - p:
                candidates[pair] = paths[pair][0]
                ms[pair] = max_f(ending, fun)
        if len(ms) > 0:
            new_pair, _ = min(ms.items())
            cancel_along_path(V, C, sigma, candidates[new_pair])  # paths[(sigma, new_sigma)][0]))


def extract(K, fun, p):
    V, C = extract_raw(K, fun, infty)
    for j in range(2, 4):
        extract_cancel(K, fun, p, j, V, C)
    return V, C


def scalar_product(sx1, sx2):
    """Returns the coefficient of sx2 in the image of sx1 withe the usual boundary
    operator."""
    d = set(sx1).difference(set(sx2)).pop()
    if list(sx1).index(d) % 2 == 0:
        return 1
    else:
        return -1


def same_orientation(path):
    """If sx1 and sx2 are simplices and 'path' a path between them, this function
    returns weather the orientation preserves while moving along 'path'"""
    o = 1
    for i in range(0, len(path)-1, 2):
        o = o * scalar_product(path[i+1], path[i]) * scalar_product(path[i+1], path[i+2])
    return o


def morse_complex(V, C):
    """Given a pair (V, C), which defines a gradient vector field, the function
    returns a dictionary 'boundary' with critical simplexes as keys.
    boundary[c] = a list of pairs (critical simplex, coefficient) that describes
    the image of morse boundary operator on c."""
    boundary = defaultdict(list)
    paths = find_paths_between_critical_simplices(V, C)
    for pair in paths:
        coeff = 0
        for path in paths[pair]:
            coeff += same_orientation(path)
        if pair[0] in boundary and coeff != 0:
            boundary[pair[0]].append((pair[1], coeff))
        elif coeff != 0:
            boundary[pair[0]] = [(pair[1], coeff)]
    return boundary


torus = [((0, 0), (1, 0), (1, 1)), ((0, 0), (1, 1), (0, 1)), ((1, 0), (2, 0), (2, 1)),
         ((1, 0), (2, 1), (1, 1)), ((2, 0), (0, 0), (0, 1)), ((2, 0), (0, 1), (2, 1)),
         ((0, 1), (1, 1), (1, 2)), ((0, 1), (1, 2), (0, 2)), ((1, 1), (2, 1), (2, 2)),
         ((1, 1), (2, 2), (1, 2)), ((2, 1), (0, 1), (0, 2)), ((2, 1), (0, 2), (2, 2)),
         ((0, 2), (1, 2), (1, 0)), ((0, 2), (1, 0), (0, 0)), ((1, 2), (2, 0), (1, 0)),
         ((1, 2), (2, 2), (2, 0)), ((2, 2), (0, 2), (0, 0)), ((2, 2), (0, 0), (2, 0))]
T = generate_all_sxs(torus)
f = {(0, 0): 0, (1, 0): 80, (2, 0): 50, (0, 1): 10, (1, 1): 60, (2, 1): 30, (0, 2): 70,
     (1, 2): 20, (2, 2): 40}


# T = generate_all_sxs([((1, 1), (2, 2), (3, 0)), ((2, 2), (3, 0), (7, 0.5))])
# f = {(1, 1): 4, (2, 2): 5, (3, 0): 10, (7, 0.5): 3}
# print(extract(T, f, infty))
# extract_cancel(T, f, infty, 3, V, C)
#print()
