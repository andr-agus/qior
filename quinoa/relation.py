"""
"""

import math

import numpy as np
import qutip

from . import state

class IORelation:
    """
    """

    def __init__(self, matrix, cutoffs):
        if not self.is_unitary(matrix):
            raise ValueError("input-output relations must be unitary")
        if not len(cutoffs) == 2:
            raise ValueError("cutoffs on all four modes must be provided")
        self.cutoffs = cutoffs
        self.matrix = matrix

    @classmethod
    def with_reflectivity(cls, R, cutoffs):
        array = np.matrix([[math.sqrt(1 - R), math.sqrt(R)],
                           [-math.sqrt(R), math.sqrt(1 - R)]])
        return cls(array, cutoffs)

    @staticmethod
    def is_unitary(matrix):
        return np.allclose(np.eye(matrix.shape[0]), matrix.H * matrix)

    def time_evolution_operator(self):
        U = 0 * qutip.tensor(qutip.qeye(self.cutoffs[0]), qutip.qeye(self.cutoffs[1]))
        for n1 in range(self.cutoffs[0]):
            for n2 in range(self.cutoffs[1]):
                bra = qutip.tensor(qutip.basis(self.cutoffs[0], n1), qutip.basis(self.cutoffs[1], n2)).dag()
                ket = self.evolve_photons(n1, n2)
                U += ket * bra
        return U

    def evolve(self, initial_state):
        return self.U * initial_state

    @property
    def U(self):
        u = getattr(self, "_U", None)
        if u is None:
            self._U = self.time_evolution_operator()
            return self._U
        return u

    def evolve_photons(self, n1, n2):
        if n1 == 0 and n2 == 0:
            return qutip.tensor(qutip.basis(self.cutoffs[0]), qutip.basis(self.cutoffs[1]))
        elif not n1 == 0:
            return self.output1() / math.sqrt(n1) * self.evolve_photons(n1 - 1, n2)
        elif not n2 == 0:
            return self.output2() / math.sqrt(n2) * self.evolve_photons(n1, n2 - 1)

    def input1(self):
        return qutip.tensor(qutip.create(self.cutoffs[0]), qutip.qeye(self.cutoffs[1]))

    def input2(self):
        return qutip.tensor(qutip.qeye(self.cutoffs[0]), qutip.create(self.cutoffs[1]))

    def output1(self):
        return self.matrix[0, 0] * self.input1() + self.matrix[0, 1] * self.input2()

    def output2(self):
        return self.matrix[1, 0] * self.input1() + self.matrix[1, 1] * self.input2()
