###############################################################################
#   IMPORTS
###############################################################################

# Python
import math
import sys

from typing import (
    Set,
)

import numpy as np  # type: ignore

from scipy import optimize  # type: ignore

# Sage
import sage.all  # type: ignore # noqa: F401

from sage.graphs.graph import Graph  # type: ignore

# Local
from association_schemes.association import AssociationScheme
from association_schemes.p_polynomial import PPolynomialScheme

###############################################################################
#   LINEAR PROGRAMMING (LP) BOUND
###############################################################################

def general_lp_bound(scheme: AssociationScheme, clique: Set[int]) -> int:
    Q = scheme.Q().numpy()
    nonzero_clique = set(clique) - {0}
    # Setting up the dual LP
    objective = Q[0, 1:]
    coefficients = Q[sorted(nonzero_clique), 1:]
    constraints = -np.ones(len(nonzero_clique))
    # Note that in the scipy implementation,
    # one *minimizes* the `objective @ x`
    # subject to `coefficients @ x <= constraints`.
    result = optimize.linprog(
        objective,
        coefficients,
        constraints,
    )
    if not result.success:
        # FIXME Change exception type
        raise Exception(result.message)
    # Reducing the precision to avoid fp-errors
    # where `floor` takes a value from `opt - eps` to `opt - 1`.
    return 1 + math.floor(np.float16(result.fun))

def lp_bound(graph: Graph) -> int:
    scheme = PPolynomialScheme.from_graph(graph)
    clique = set(range(2, scheme.diameter() + 1))
    return general_lp_bound(scheme, clique)

###############################################################################
#   HOFFMAN BOUND
###############################################################################

def hoffman_bound(graph: Graph) -> int:
    valency = graph.degree()[0]
    eigval = min(graph.adjacency_matrix().eigenvalues())
    return math.floor(len(graph.vertices()) * 1 / (1 - valency / eigval))
