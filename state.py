"""
"""

import numpy as np

class State:
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
        raise NotImplementedError("to be implemented in derived classes")

class Ket(State):

    @classmethod
    def shape(cls):
        D = 1
        for d in cls.cutoffs:
            D *= d
        return (D,)

    def dagger(self):
        return Bra(self.array.H)

class Bra(Ket):

    @classmethod
    def shape(cls):
        return (1,super(cls).shape())

    def dagger(self):
        return Ket(self.array.H)

class DensityMatrix(Ket)
    """
    """

    @classmethod
    def shape(cls):
        super(cls).shape()*2

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
        return map(lambda x: (x[1], x[0]), self.array.ndenumerate())
