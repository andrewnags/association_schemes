
import abc

from typing import (
    cast,
    Collection,
    Generic,
    Optional,
    Sequence as Seq,
    Set,
    Type,
    TypeVar,
    Union,
)

# Sage
import sage.all  # type: ignore # noqa: F401

from sage.combinat.set_partition import SetPartition  # type: ignore
from sage.matrix.constructor import matrix  # type: ignore
from sage.modules.free_module_element import vector  # type: ignore
from sage.structure.element import (  # type: ignore
    Matrix,
    Vector,
)

T = TypeVar("T", covariant=True)

def _subset_from_vector(
    cv: Vector,
    elements: Collection[T],
    check: bool = False,
) -> Set[T]:
    if check and len(elements) != len(cv):
        raise TypeError(
            "The number of elements must be the "
            + "length of the characteristic vector."
        )
    if check and not set(cv) <= {0, 1}:
        raise ValueError("The characteristic vector must be a 01 vector.")
    return set(cv.support())

def _vector_from_subset(
    subset: Collection[T],
    supset: Collection[T],
    check: bool = False,
) -> Vector:
    if check and not set(subset) <= set(supset):
        raise ValueError(
            "The given elements are not a superset of the subset."
        )
    return vector([
        1 if elt in subset
        else 0
        for elt in supset
    ])


def set_partition_from_matrix(
    matrix: Matrix,
) -> SetPartition:
    elements = list(range(len(matrix.rows())))
    return SetPartition(
        _subset_from_vector(row, elements)
        for row in matrix
    )

def characteristic_matrix(sp: SetPartition) -> Matrix:
    supset = sorted(sp.base_set())
    return matrix([
        _vector_from_subset(part, supset)
        for part in sp
    ]).T
