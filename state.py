"""
"""

import numpy as np

class DensityMatrix:
    """
    """

    cutoffs = [5]*4

    @classmethod
    def set_cutoff(cls, n1, n2, n3, n4):
        cls.cutoffs = [n1, n2, n3, n4]

    def __init__(self, array):
        self.validate_dims(array)
        self.array = array

    @classmethod
    def validate_dims(cls, array):
        if not array.shape == cls.shape():
            raise ValueError("Array dimensions do not match with %s" % cls)

    @classmethod
    def shape(cls):
        dims = 1
        for cut in cls.cutoffs:
            dims *= cut
        return (dims,)*2

    @classmethod
    def index2photons(cls, i):
        ns = []
        for dim in range(len(cls.cutoffs)):
            ns.append(i % cls.cutoffs[dim])
            i = i // cls.cutoffs[dim]
        return ns

    @classmethod
    def zero(cls):
        array = np.zeros(cls.shape())
        return cls(array)

    @classmethod
    def vacuum(cls):
        self = cls.zero()
        self.array[0,0] = 1
        return self

    def ketbras(self):
        return map(lambda x: (x[0], (self.index2photons(x[1][0]), self.index2photons(x[1][1]))), map(lambda x: (x[1], x[0]), np.ndenumerate(self.array)))

    @classmethod
    def ketbra(self, ketbra_tuple):
        pass
