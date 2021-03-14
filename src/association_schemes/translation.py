###############################################################################
#   IMPORTS
###############################################################################

# Python
import functools as ft

from typing import (
    Callable,
    cast,
    Dict,
    List,
    Optional,
    Sequence as Seq,
    Set,
    Tuple,
    Type,
    TypeVar,
)

# Sage
import sage.all  # type: ignore # noqa: F401

from sage.graphs.graph import Graph  # type: ignore
from sage.groups.group import Group  # type: ignore
from sage.matrix.constructor import matrix  # type: ignore
from sage.rings.real_double import RDF  # type: ignore
from sage.structure.element import (  # type: ignore
    CommutativeRingElement,
    Matrix,
    MultiplicativeGroupElement,
    Vector,
)

# Local
import association_schemes.association as asch


Part = Set[MultiplicativeGroupElement]

###############################################################################
#   ASSOCIATION CLASSES
###############################################################################

class CayleyClass(asch.AssociationClass):
    Self = TypeVar("Self", bound="CayleyClass")

    _group: Group
    __part: Optional[Part]

    def __init__(self, group: Group, part: Part):
        part = set(part)
        if part is None or part == {group.identity()}:
            self.__part = None
        else:
            self.__part = part
        self._group = group

    def group(self) -> Group:
        return self._group

    @ft.lru_cache(None)
    def part(self) -> Part:
        if self.__part is None:
            return {self._group.identity()}
        return set(self.__part) | {
            elt.inverse()
            for elt in self.__part
        }

    @ft.lru_cache(None)
    def charvec(self) -> List[int]:
        part = self.part()
        return [
            1 if elt in part
            else 0
            for elt in self.group().conjugacy_classes_representatives()
        ]

    @classmethod
    def from_class(cls: Type[Self], class_: asch.AssociationClass) -> Self:
        if isinstance(class_, cls):
            return class_
        raise NotImplementedError(
            "Conversion from `AssociationClass` to `CayleyClass`"
            + " not yet implemented."
        )

    def adjacency(self) -> Matrix:
        return self.graph().adjacency_matrix()

    @ft.lru_cache(None)
    def graph(self) -> Graph:
        if self.__part is None:
            g = Graph(loops=True)
            g.add_edges([
                (v, v)
                for v in self._group
            ])
            return g
        # This will automatically take
        # the closure under inversion.
        return self._group.cayley_graph(
            generators=self.__part
        ).to_undirected()


###############################################################################
#   ASSOCIATION SCHEMES
###############################################################################


# FIXME Inexactness sensitivity?
def group_rows_by_eq(mat: Matrix) -> Seq[Seq[int]]:
    classes: Dict[Vector, List[int]] = dict()
    for i, _ in enumerate(mat):
        row = mat[i]  # Do this so the rows will be hashable
        classes[row] = classes.get(row, []) + [i]
    return [
        indices
        for _, indices in asch.sort_decomp(
            (row[1], indices)
            for row, indices in classes.items()
        )
    ]

# TODO Use this in p_polynomial/association?
def basis_to_proj(basis: Matrix) -> Matrix:
    basis = basis.change_ring(RDF).gram_schmidt(True)[0]
    return basis.T * basis

T = TypeVar("T")

def precompute_character_partition(
    method: Callable[..., T]
) -> Callable[..., T]:
    @ft.wraps(method)
    def precompute(self: "TranslationScheme", *args, **kwargs) -> T:
        if not getattr(self, "_character_partition", []):
            self._characters = self.group().character_table()
            self._pre_P = self._characters * self.partition_matrix()
            self._character_partition = group_rows_by_eq(self._pre_P)
            if len(self._character_partition) != self.diameter() + 1:
                raise ValueError(
                    "The group partition does not induce the same number of "
                    + "cells in the dual partition, so this is not a valid "
                    + "association scheme."
                )
        return method(self, *args, **kwargs)
    return precompute

# TODO Derive from P-polynomial?
class TranslationScheme(asch.AssociationScheme):
    Self = TypeVar("Self", bound="TranslationScheme")

    _group: Group
    # The character table of the group
    _characters: Matrix
    # The matrix of eigenvalues, with some rows duplicated
    _pre_P: Matrix
    # Partition of the characters into eigenspaces
    _character_partition: Seq[Seq[int]]

    def validate(self, raise_=True) -> bool:
        raise NotImplementedError()

    @classmethod
    def from_scheme_with_group(
        cls: Type[Self],
        scheme: asch.AssociationScheme,
        group: Group,
    ) -> Self:
        # TODO Maybe cache some pre-computed values from `scheme`?
        ts = cls()
        ts.__classes = [
            CayleyClass.from_class(class_)
            for class_ in scheme.classes()
        ]
        ts._group = group
        ts._character_partition = []
        return ts

    @classmethod
    def from_group_partition(
        cls: Type[Self],
        group: Group,
        partition: Seq[Part],
    ) -> Self:
        return cls.from_scheme_with_group(
            asch.AssociationScheme.from_classes([
                CayleyClass(group, part)
                for part in partition
            ]),
            group,
        )

    def classes(self) -> Seq[CayleyClass]:
        return cast(Seq[CayleyClass], self.__classes)

    def group(self) -> Group:
        return self._group

    def partition(self) -> Seq[Part]:
        return [
            class_.part()
            for class_ in self.classes()
        ]

    def partition_matrix(self) -> Matrix:
        return matrix(
            class_.charvec()
            for class_ in self.classes()
        ).T

    @ft.lru_cache(None)
    @precompute_character_partition
    def dual_basis(self) -> Seq[Tuple[int, Matrix]]:
        # Note that the rows of the real parts of the character table
        # provide a real basis for the eigenspace.
        real_characters = self._characters.apply_map(
            lambda z: getattr(z, "real", lambda :z)()
        )
        return [
            (len(indices), basis_to_proj(real_characters[indices, :]))
            for indices in self._character_partition
        ]

    @ft.lru_cache(None)
    @precompute_character_partition
    def P(self) -> Matrix:
        prep = self._pre_P
        cp = self._character_partition
        # Deduplicate rows
        rows = [
            prep[indices[0]]
            for indices in cp
        ]
        return matrix(rows)
