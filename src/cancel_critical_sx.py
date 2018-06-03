from readwrite import read_dgvf_from_file
from random import choice
from itertools import chain, combinations
from collections import defaultdict
from expand_function import extract_raw


def boundary(a):
    """Return the boundary (codimension one faces) of the simplex a.
    """
    return combinations(a, len(a)-1)


def reverse_path(path):
    return [path[i] for i in range(len(path)-1, -1, -1)]


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
        print(len(paths))
        candidates = [(pair, path[0]) for pair, path in paths.items() if len(path) == 1]
        if len(candidates) == 0:
            break
        pair, path = choice(candidates)
        cancel_along_path(V, C, pair[0], path)


def max_f(sx, fun):
    return max([fun[s] for s in list(sx)])


def extract_cancel(K, fun, p, j, V, C):
    paths = find_paths_between_critical_simplices(V, [c for c in C if len(c) in [j, j-1]])
    for sigma in [c for c in C if len(c) == j]:
        candidates = {}
        for pair in paths:
            starting, ending = pair
            if starting == sigma and len(paths[pair]) == 1 and max_f(ending, fun) > max_f(starting, fun) - p:
                candidates[pair] = paths[pair][0]
        ms = {}
        for pair1 in candidates:
            _, y = pair1
            special = True
            for pair2 in candidates:
                _, z = pair2
                if y == z:
                    special = False
            if special:
                ms[y] = max_f(y, fun)
        if len(ms) > 0:
            new_sigma, _ = min(ms.items())
            cancel_along_path(V, C, new_sigma, reverse_path(paths[(sigma, new_sigma)][0]))


def extract(K, fun, p, read_from_file=True):
    if read_from_file:
        V, C = read_dgvf_from_file("discrete_vector_field.txt")
    else:
        V, C = extract_raw(K, fun)
    for j in range(1, 3):
        extract_cancel(K, fun, p, j, V, C)
    return V, C

# extract(generate_all_sxs(T), f, 4)
