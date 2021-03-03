###############################################################################
#   IMPORTS
###############################################################################

# Python
import abc
import functools as ft
import operator as op
from typing import (
    List,
    Tuple,
)

import numpy as np

# Sage
import sage.all  # noqa: F401

from sage.graphs.graph import Graph
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix
#from sage.rings.qqbar import QQbar  # noqa: E265
from sage.rings.real_double import RDF
from sage.structure.element import Matrix

mult_inv = np.vectorize(ft.partial(op.truediv, 1))

###############################################################################
#   SPECTRAL DECOMPOSITION
###############################################################################

def sort_decomp(decomp):
    return sorted(decomp, key=lambda x: x[0], reverse=True)

def group_eigspaces(vals, vecs):
    d = dict()
    for i, val in enumerate(vals):
        d[val] = d.get(val, []) + [i]
    return (
        (val, len(inds), vecs[:, inds])
        for val, inds in d.items()
    )

def spec_decomp_fast(A, float_type=None):
    A = np.array(A, dtype="uint8")
    assert (A == A.T).all()
    vals, vecs = np.linalg.eigh(A)
    if float_type is not None:
        vals = vals.astype(float_type)
        vecs = vecs.astype(float_type)
    return group_eigspaces(vals, vecs)

def spec_decomp_exact(A):
    A = matrix(A)
    return (
        (
            val,
            space.dimension(),
            space.basis_matrix()
                # FIXME Using QQbar here causes a `RecursionError`
                # to occur when computing square roots
                # for `gram_schmidt` below.
                #.change_ring(QQbar)  # noqa: E265
                .change_ring(RDF)
                .gram_schmidt(True)[0].T,
        )
        for val, space in A.right_eigenspaces()
    )

def spec_decomp(A, exact=True, float_type=None):
    decomp = None
    if exact:
        decomp = spec_decomp_exact(A)
    else:
        decomp = spec_decomp_fast(A, float_type)
    return sort_decomp(decomp)

###############################################################################
#   ASSOCIATION CLASSES
###############################################################################

class AssociationClass(abc.ABC):

    @abc.abstractmethod
    def adjacency(self) -> Matrix:
        pass

    @abc.abstractmethod
    def graph(self) -> Graph:
        pass

class AssociationGraph(AssociationClass):

    def __init__(self, graph):
        self.__graph = Graph(graph).to_undirected()
        self.__graph.remove_multiple_edges()

    @ft.lru_cache(None)
    def adjacency(self) -> Matrix:
        return self.__graph.adjacency_matrix()

    def graph(self) -> Graph:
        return self.__graph

class AssociationMatrix(AssociationClass):

    def __init__(self, adj_mat):
        self.__adjacency = matrix(adj_mat)

    def adjacency(self) -> Matrix:
        return self.__adjacency

    @ft.lru_cache(None)
    def graph(self) -> Graph:
        return Graph(self.__adjacency)


###############################################################################
#   ASSOCIATION SCHEMES
###############################################################################

def columnize_matrices(mats: List[Matrix]) -> Matrix:
    return matrix([
        mat.list()
        for mat in mats
    ]).T

def matrix_of_eigenvalues(
    adjacencies: List[Matrix],
    idempotents: List[Matrix],
) -> Matrix:
    adjs_mat = columnize_matrices(adjacencies)
    idem_mat = columnize_matrices(idempotents)
    return idem_mat.solve_right(adjs_mat)

class AssociationScheme(object):

    @staticmethod
    def from_adjacency(adj_mat) -> "AssociationScheme":
        assoc_scheme = AssociationScheme()
        assoc_scheme.root = AssociationMatrix(adj_mat)
        return assoc_scheme

    @staticmethod
    def from_graph(graph) -> "AssociationScheme":
        assoc_scheme = AssociationScheme()
        assoc_scheme.root = AssociationGraph(graph)
        return assoc_scheme

    # TODO
    def validate(self, raise_=True) -> bool:
        return True

    @ft.lru_cache(None)
    def len(self) -> int:
        return len(self.root.graph())

    @ft.lru_cache(None)
    def diameter(self) -> int:
        return self.root.graph().diameter()

    @ft.lru_cache(None)
    def classes(self) -> List[AssociationClass]:
        dist_mat = self.root.graph().distance_matrix().numpy()
        return [
            AssociationMatrix(np.array(dist_mat == dist, dtype="uint8"))
            for dist in range(self.diameter() + 1)
        ]

    def basis(self) -> List[Tuple[int, Matrix]]:
        return list(zip(self.degrees(), self.adjacencies()))

    def adjacencies(self) -> List[Matrix]:
        return [
            assoc_class.adjacency()
            for assoc_class in self.classes()
        ]

    @ft.lru_cache(None)
    def dual_basis(self) -> List[Tuple[int, Matrix]]:
        decomp = spec_decomp(self.root.adjacency())
        if len(decomp) != self.diameter() + 1:
            raise ValueError(
                "The spectral decomposition does not contain"
                + " the correct number of eigenspaces."
            )
        return [
            (dim, basis_mat * basis_mat.T)
            for (_, dim, basis_mat) in decomp
        ]

    def idempotents(self) -> List[Matrix]:
        return [
            proj_matrix
            for (_, proj_matrix) in self.dual_basis()
        ]

    @ft.lru_cache(None)
    def P(self) -> Matrix:
        return matrix_of_eigenvalues(
            self.adjacencies(),
            self.idempotents(),
        )

    def degrees(self) -> List:
        return self.P()[0, :].list()

    def multiplicities(self) -> List:
        return [
            dim
            for (dim, _) in self.dual_basis()
        ]

    @ft.lru_cache(None)
    def Q(self) -> Matrix:
        return (
            diagonal_matrix(mult_inv(self.degrees()))
            * self.P().T
            * diagonal_matrix(self.multiplicities())
        )
