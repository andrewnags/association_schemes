import functools as ft
import itertools as it
import logging
import operator as op
import sys

# Sage
import sage.all

from sage.groups.abelian_gps.abelian_group_gap import (  # type: ignore
        AbelianGroupGap
)
from sage.groups.perm_gps.permgroup import direct_product_permgroups
from sage.groups.perm_gps.permgroup_named import CyclicPermutationGroup

# Local
import association_schemes.association as assoc
import association_schemes.p_polynomial as ppoly
import association_schemes.translation as trans
import association_schemes.bounds as bounds

def Zqn(q, n):
    return AbelianGroupGap([q]*n)

# FIXME q != 2 not working
def hamming_generators(n, q):
    logging.info("Getting generators for H({}, {})".format(n, q))
    group = Zqn(q, n)
    partition = [[] for _ in range(n + 1)]
    for elt in group.conjugacy_classes_representatives():
        partition[sum(map(bool, elt.exponents()))].append(elt)
    return group, partition

def hamming_scheme(n, q):
    logging.info("Making scheme of H({}, {})".format(n, q))
    return trans.TranslationScheme.from_group_partition(
        *hamming_generators(n, q)
    )

def bound_scheme_cocliques(scheme):
    logging.info("Getting LP coclique bound.")
    clique = set(range(2, scheme.diameter() + 1))
    try:
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
        scheme = hamming_scheme(*param)
        yield (
            "translation",
            "H" + str(param),
            bound_scheme_cocliques(scheme),
        )
        eigvals = scheme.P().T[1]
        yield (
            "trans-ratio",
            "H" + str(param),
            scheme.len() / (1 - eigvals[0] / eigvals[-1])
        )

if __name__ == "__main__":
    print("bound,graph,value")
    for bound, graph, value in run_experiments(experiments):
        print(f"{bound},\"{graph}\",{value}")
