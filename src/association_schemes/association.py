###############################################################################
#   IMPORTS
###############################################################################

# Python
import abc
import functools as ft
import operator as op
from typing import (
    Dict,
    Iterable,
    List,
    Sequence as Seq,
    Tuple,
    Type,
    TypeVar,
)

import numpy as np  # type: ignore

# Sage
import sage.all  # type: ignore # noqa: F401

from sage.graphs.graph import Graph  # type: ignore
from sage.matrix.constructor import matrix  # type: ignore
from sage.matrix.special import diagonal_matrix  # type: ignore
#from sage.rings.qqbar import QQbar  # noqa: E265
from sage.rings.real_double import RDF  # type: ignore
from sage.structure.element import (  # type: ignore
    CommutativeRingElement,
    Matrix,
)

mult_inv = np.vectorize(ft.partial(op.truediv, 1))

###############################################################################
#   SPECTRAL DECOMPOSITION
###############################################################################


def sort_decomp(
    decomp: Iterable[Tuple[CommutativeRingElement, ...]],
) -> Seq[Tuple[CommutativeRingElement, ...]]:
    return sorted(decomp, key=lambda x: x[0], reverse=True)

def group_eigspaces(
    vals: np.ndarray,
    vecs: np.ndarray,
) -> Iterable[Tuple[CommutativeRingElement, int, Matrix]]:
    d: Dict[CommutativeRingElement, List[int]] = dict()
    for i, val in enumerate(vals):
        d[val] = d.get(val, []) + [i]
    return (
        (val, len(inds), vecs[:, inds])
        for val, inds in d.items()
    )

def spec_decomp_fast(
    A,
    float_type=None,
) -> Iterable[Tuple[CommutativeRingElement, int, Matrix]]:
    A = np.array(A, dtype="uint8")
    assert (A == A.T).all()
    vals, vecs = np.linalg.eigh(A)
    if float_type is not None:
        vals = vals.astype(float_type)
        vecs = vecs.astype(float_type)
    return group_eigspaces(vals, vecs)

def spec_decomp_exact(A) -> Iterable[Tuple[CommutativeRingElement, int, Matrix]]:
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

# TODO Type this
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

def columnize_matrices(mats: Seq[Matrix]) -> Matrix:
    return matrix([
        mat.list()
        for mat in mats
    ]).T

def matrix_of_eigenvalues(
    adjacencies: Seq[Matrix],
    idempotents: Seq[Matrix],
) -> Matrix:
    adjs_mat = columnize_matrices(adjacencies)
    idem_mat = columnize_matrices(idempotents)
    return idem_mat.solve_right(adjs_mat)

class AssociationScheme(object):
    Self = TypeVar("Self", bound="AssociationScheme")

    __classes: Seq[AssociationClass]

    @classmethod
    def from_classes(cls: Type[Self], classes: Seq[AssociationClass]) -> Self:
        assoc_scheme = cls()
        assoc_scheme.__classes = list(classes)
        if len(classes) < 2:
            raise ValueError(
                "Need at least two classes to create an association scheme"
            )
        return assoc_scheme

    def classes(self) -> Seq[AssociationClass]:
        return self.__classes

    def diameter(self) -> int:
        return len(self.classes()) - 1

    def len(self) -> int:
        return len(self.classes()[1].graph())

    # TODO
    def validate(self, raise_: bool = True) -> bool:
        raise NotImplementedError()
        return True

    def basis(self) -> Seq[Tuple[int, Matrix]]:
        return list(zip(self.degrees(), self.adjacencies()))

    def adjacencies(self) -> Seq[Matrix]:
        return [
            assoc_class.adjacency()
            for assoc_class in self.classes()
        ]

    # TODO
    def dual_basis(self) -> Seq[Tuple[int, Matrix]]:
        # If the scheme is *not* P-polynomial
        # (i.e. derived from a DRG), then the
        # idempotents can't just be computed from
        # the "root" class.
        raise NotImplementedError()

    def idempotents(self) -> Seq[Matrix]:
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

    def degrees(self) -> Seq:
        return self.P()[0, :].list()

    def multiplicities(self) -> Seq:
        return [
            dim
            for (dim, _) in self.dual_basis()
        ]

    @ft.lru_cache(None)
    def Q(self) -> Matrix:
        deg = self.degrees()
        mult = self.multiplicities()
        return (
            diagonal_matrix(mult_inv(deg))
            * self.P().T
            * diagonal_matrix(mult)
        )
