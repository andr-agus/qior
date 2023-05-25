"""
"""

import math

import numpy as np
import qutip

class IORelation:
    """
    """

    def __init__(self, matrix, dims, acting_on):
        """
        """

        if not self.is_unitary(matrix):
            raise ValueError("input-output relations must be unitary")

        if len(dims) < 2:
            raise ValueError("this is not a bigger system")

        if not len(acting_on) == 2:
            raise ValueError("IORelations act on two systems, not less, not more")

        if len(dims) == 2:
            if not acting_on == (0, 1):
                raise ValueError()

        if acting_on[0] < 0 or len(dims) <= acting_on[0]:
            raise ValueError()

        if acting_on[1] < 0 or len(dims) <= acting_on[1]:
            raise ValueError()

        self.matrix = matrix
        self.dims = dims
        self.acting_on = acting_on
        self._U = self.time_evolution_operator()

    def __call__(self, state):
        return self.evolve(state)

    @classmethod
    def with_reflectivity(cls, R, dims, acting_on):
        array = np.matrix([[math.sqrt(1 - R), math.sqrt(R)],
                           [-math.sqrt(R), math.sqrt(1 - R)]])
        return cls(array, dims, acting_on)

    @staticmethod
    def is_unitary(matrix):
        return np.allclose(np.eye(matrix.shape[0]), matrix.H * matrix)

    @property
    def U(self):
        u = getattr(self, "_U", None)
        if u is None:
            self._U = self.time_evolution_operator()
            return self._U
        return u

    def time_evolution_operator(self):
        U = 0 * qutip.tensor(qutip.qeye(self.dims[self.acting_on[0]]), qutip.qeye(self.dims[self.acting_on[1]]))
        for n1 in range(self.dims[self.acting_on[0]]):
            for n2 in range(self.dims[self.acting_on[1]]):
                bra = qutip.tensor(
                    qutip.basis(self.dims[self.acting_on[0]], n1),
                    qutip.basis(self.dims[self.acting_on[1]], n2)).dag()
                ket = self.evolve_photons(n1, n2)
                U += ket * bra
        return self.expand_to_bigger_system(U)

    def evolve_photons(self, n1, n2):
        if n1 == 0 and n2 == 0:
            return qutip.tensor(qutip.basis(self.dims[self.acting_on[0]]), qutip.basis(self.dims[self.acting_on[1]]))
        elif not n1 == 0:
            return self.output1() / math.sqrt(n1) * self.evolve_photons(n1 - 1, n2)
        elif not n2 == 0:
            return self.output2() / math.sqrt(n2) * self.evolve_photons(n1, n2 - 1)

    def output1(self):
        return self.matrix[0, 0] * self.input1() + self.matrix[0, 1] * self.input2()

    def output2(self):
        return self.matrix[1, 0] * self.input1() + self.matrix[1, 1] * self.input2()

    def input1(self):
        return qutip.tensor(qutip.create(self.dims[self.acting_on[0]]), qutip.qeye(self.dims[self.acting_on[1]]))

    def input2(self):
        return qutip.tensor(qutip.qeye(self.dims[self.acting_on[0]]), qutip.create(self.dims[self.acting_on[1]]))

    def expand_to_bigger_system(self, U):
        """
        """
        new_dims = self.permute_dims()

        for D in new_dims[2:]:
            U = qutip.tensor(U, qutip.qeye(D))

        return self.permute_back_unitary(U)

    def permute_dims(self):
        new_dims = list(self.dims)
        new_dims[0] = self.dims[self.acting_on[0]]
        new_dims[self.acting_on[0]] = self.dims[0]
        new_new_dims = list(new_dims)
        new_new_dims[1] = new_dims[self.acting_on[1]]
        new_new_dims[self.acting_on[1]] = new_dims[1]
        return new_new_dims

    def permute_back_unitary(self, U):
        identity = list(range(len(self.dims)))
        first_permutation = list(identity)
        first_permutation[1] = identity[self.acting_on[1]]
        first_permutation[self.acting_on[1]] = identity[1]
        second_permutation = list(first_permutation)
        second_permutation[0] = first_permutation[self.acting_on[0]]
        second_permutation[self.acting_on[0]] = first_permutation[0]
        return U.permute(second_permutation)

    def evolve(self, initial_state):
        if self.is_pure(initial_state):
            return self.U * initial_state
        else:
            return self.U * initial_state * self.U.dag()

    @classmethod
    def is_pure(cls, state):
        return not state.dims[0] == state.dims[1]
