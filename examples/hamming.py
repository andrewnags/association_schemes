import functools as ft
import itertools as it
import operator as op
import sys

# Sage
import sage.all

from sage.groups.perm_gps.permgroup import direct_product_permgroups
from sage.groups.perm_gps.permgroup_named import CyclicPermutationGroup

# Local
import association_schemes.association as assoc
import association_schemes.p_polynomial as ppoly
import association_schemes.translation as trans
import association_schemes.bounds as bounds

def Zqn(q, n):
    Zq = CyclicPermutationGroup(q)
    return direct_product_permgroups([Zq]*n)

# FIXME q != 2 not working
def hamming_generators(n, q):
    group = Zqn(q, n)
    partition = [[] for _ in range(n + 1)]
    for elt in group.conjugacy_classes_representatives():
        partition[len(elt.cycles())].append(elt)
    return group, partition

def bound_hamming(n, q):
    try:
        scheme = trans.TranslationScheme.from_group_partition(
            *hamming_generators(n, q)
        )
        clique = set(range(2, scheme.diameter() + 1))
        return bounds.general_lp_bound(scheme, clique)
    except Exception:
        return "fail"

experiments = [
    (4, 2), (4, 3),
    (5, 2), (5, 3),
]

def run_experiments(params):
    for param in params:
        print("translation", param, file=sys.stderr)
        yield (
            "translation",
            "H" + str(param),
            bound_hamming(*param),
        )

if __name__ == "__main__":
    print("bound,graph,value")
    for bound, graph, value in run_experiments(experiments):
        print(f"{bound},\"{graph}\",{value}")
