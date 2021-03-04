###############################################################################
#   IMPORTS
###############################################################################

# Python
import functools as ft

from typing import (
    List,
    Tuple,
    Type,
    TypeVar,
)

import numpy as np  # type: ignore

# Sage
import sage.all  # type: ignore # noqa: F401

from sage.structure.element import Matrix  # type: ignore

# Local
import association_schemes.association as asch

def make_classes(root: asch.AssociationClass) -> List[asch.AssociationClass]:
    dist_mat = root.graph().distance_matrix().numpy()
    diameter = dist_mat.max()
    return [
        asch.AssociationMatrix(np.array(dist_mat == dist, dtype="uint8"))
        for dist in range(diameter + 1)
    ]

class PPolynomialScheme(asch.AssociationScheme):
    Self = TypeVar("Self", bound="PPolynomialScheme")

    @classmethod
    def from_root(cls: Type[Self], root: asch.AssociationClass) -> Self:
        return cls.from_classes(
            make_classes(root)
        )

    @classmethod
    def from_adjacency(cls: Type[Self], adj_mat: Matrix) -> Self:
        return cls.from_root(
            asch.AssociationMatrix(adj_mat)
        )

    @classmethod
    def from_graph(cls: Type[Self], graph) -> Self:
        return cls.from_root(
            asch.AssociationGraph(graph)
        )

    def validate(self, raise_=True) -> bool:
        valid = self.classes()[1].graph().is_distance_regular()
        if raise_ and not valid:
            raise ValueError(
                "The first (note: not zeroeth) class in the scheme"
                + " must be a distance-regular graph"
            )
        return valid

    @ft.lru_cache(None)
    def dual_basis(self) -> List[Tuple[int, Matrix]]:
        # We can compute the idempotents more
        # efficiently in the DRG/P-polynomial case.
        decomp = asch.spec_decomp(self.classes()[1].adjacency())
        if len(decomp) != self.diameter() + 1:
            raise ValueError(
                "The spectral decomposition does not contain"
                + " the correct number of eigenspaces."
            )
        return [
            (dim, basis_mat * basis_mat.T)
            for (_, dim, basis_mat) in decomp
        ]
