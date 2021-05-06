import itertools as it

import pytest
import sage.all
from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap

from association_schemes import characters

EXPONENTS = [
    (2,),
    (3,),
    (4,),
    (2, 2),
    (3, 3),
    (4, 4),
    (2, 2, 2),
    (3, 3, 3),
    (4, 4, 4),
    (2, 2, 3),
    (2, 3, 2),
    (3, 2, 2),
    (2, 3, 3),
    (3, 2, 3),
    (3, 3, 2),
    (2, 3, 4),
    (2, 4, 3),
    (3, 2, 4),
    (3, 4, 2),
    (4, 2, 3),
    (4, 3, 2),
]

@pytest.mark.unit
@pytest.mark.parametrize("exponents", EXPONENTS)
def test_index_exponents(exponents):
    """ Test that `AbelianGroupElements` computes the correct element at each
    index.
    """
    actual_list = list(it.product(*map(range, exponents)))
    group = characters.AbelianGroupElements(exponents)
    for i, actual in enumerate(actual_list):
        assert group.element_at(i) == actual

def row_sums(M):
    return [sum(row) for row in M.rows()]

def coln_sums(M):
    return row_sums(M.T)

def elts_set(M):
    return {elt for vec in M for elt in vec}

def is_permutation_matrix(M):
    return (
        set(row_sums(M)) == {1}
        and set(coln_sums(M)) == {1}
        and elts_set(M) == {0, 1}
    )

class AbelianGroupElements(characters.AbelianGroupElements):

    def __init__(self, exponents):
        self.group = AbelianGroupGap(exponents)
        self.elements = self.group.conjugacy_classes_representatives()

    def element_at(self, index):
        return self.elements[index].exponents()

@pytest.mark.unit
@pytest.mark.parametrize("exponents", EXPONENTS)
def test_character_table(exponents):
    """ Test that the character table is computed correctly (up to a
    permutation of the characters/elements).
    """
    G = AbelianGroupElements(exponents)
    actual_CT = G.group.character_table()
    test_CT = characters.CharacterTable(G).matrix()
    assert is_permutation_matrix(test_CT * actual_CT**(-1))
