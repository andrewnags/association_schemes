import operator as op
import itertools as it

import numpy as np
cimport numpy as np

import sage.all
from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
from sage.matrix.constructor import matrix
from sage.rings.number_field.number_field import CyclotomicField


cdef class AbelianGroupElements:

    cdef readonly unsigned long[:] exponents
    cdef readonly unsigned long[:] moduli
    cdef readonly unsigned long n
    cdef readonly unsigned long size

    cdef readonly tuple elts

    def __cinit__(self, exponents):
        cdef list prods
        self.exponents = np.array(exponents, dtype=np.uint)
        self.n = len(exponents)
        prods = list(it.accumulate(self.exponents[::-1], op.mul))
        self.size = prods.pop()
        self.moduli = np.array([1,] + prods, dtype=np.uint)

    cpdef tuple element_at(self, unsigned long index):
        cdef unsigned long i, ni
        cdef unsigned long[:] element = np.zeros(self.n, dtype=np.uint)
        for i in range(self.n):
            ni = self.n - i - 1
            element[ni] = (index // self.moduli[i]) % self.exponents[ni]
        return tuple(element)

cdef class Partition:

    cpdef unsigned long n_elements(self):
        raise NotImplementedError()

    cpdef unsigned long n_parts(self):
        raise NotImplementedError()

    cpdef unsigned long part(self, unsigned long element):
        raise NotImplementedError()

    cpdef all_parts(self):
        return [
            Part(self, i)
            for i in range(self.n_parts)
        ]

    cpdef Product product(self, CharacterTable char_table):
        return Product(char_table, self)

ctypedef unsigned long (*partition_function)(unsigned long)

cdef class FunctionPartition(Partition):

    cdef readonly partition_function function
    cdef readonly _n_elements

    def __cinit__(self,
        partition_function function,
        unsigned long n_elements,
        unsigned long n_parts,
    ):
        self.function = function
        self._n_elements = n_elements
        self._n_parts = n_parts

    cpdef unsigned long n_elements(self):
        return self._n_elements

    cpdef unsigned long n_parts(self):
        return self._n_parts

    cpdef unsigned long part(self, unsigned long element):
        return self.function(element)

cdef class EnumeratedPartition(Partition):

    cdef readonly partition
    cdef dict lookup

    def __cinit__(self, partition)
        """
        Params:
            partition : a list of sets; each set is a part of the partition
        """
        self.partition = partition
        self.lookup = None

    cpdef unsigned long n_elements(self):
        return len(set.union(*self.partition))

    cpdef unsigned long n_parts(self):
        return len(self.partition)

    cdef _make_lookup(self):
        self.lookup = dict()
        for i, part in enumerate(self.partition):
            for elt in part:
                self.lookup[elt] = i

    cpdef unsigned long part(self, unsigned long element):
        if self.lookup is None:
            self._make_lookup()
        return self.lookup[element]

# TODO MatrixPartition type
# TODO Conversion between partition types

cdef class Part:

    cdef readonly Partition partition
    cdef readonly unsigned long part

    def __cinit__(self, Partition partition, unsigned long part):
        self.partition = partition
        self.part = part

cdef class CharacterTable:

    cdef readonly AbelianGroupElements group
    cdef readonly list fields
    cdef readonly list roots

    def __cinit__(self, AbelianGroupElements group):
        self.group = group
        self.fields = [
            CyclotomicField(i)
            for i in self.group.exponents
        ]
        self.roots = [
            field.gen()
            for field in self.fields
        ]

    cpdef index(self,
        unsigned long chrIx,
        unsigned long eltIx,
    ):
        cdef tuple character = self.group.element_at(chrIx)
        cdef tuple element = self.group.element_at(eltIx)
        cdef prod = self.roots[0].parent().one()
        for i, (chr, elt) in enumerate(zip(character, element)):
            prod *= self.roots[i]**(chr * elt)
        return prod

    cpdef matrix(self):
        return matrix([[
                self.index(i, j)
                for j in range(self.group.size)
            ] for i in range(self.group.size)
        ])

    cpdef Product product(self, Partition partition):
        return partition.product(self)

cdef class Product:

    cdef readonly CharacterTable char_table
    cdef unsigned long current_row

    def __cinit__(self, CharacterTable char_table, Partition partition):
        self.char_table = char_table
        self.partition = partition
        self.current_row = 0

    cpdef row(self, unsigned long i):
        raise NotImplementedError()

    def __iter__(self):
        return self

    def __next__(self):
        cdef row = self.row(self.current_row)
        self.current_row += 1
        return row

cdef class FunctionProduct(Product):

    cdef readonly FunctionPartition partition

    cdef row(self, unsigned long i):
        cdef unsigned long j, p
        cdef row = np.zeros(self.partition.n_parts(), dtype=object)
        for j in range(self.char_table.group.size()):
            p = self.partition.part(j)
            row[p] += self.char_table.index(i, j)
        return row

cdef class EnumeratedProduct(Product):

    cdef readonly EnumeratedPartition partition

    cpdef row(self, unsigned long i):
        return [sum(
                self.char_table.index(i, j)
                for j in part
            ) for part in self.partition.partition
        ]
