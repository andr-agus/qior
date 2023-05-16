"""
"""

import math

import numpy as np

from . import state

class IORelation:
    """
    """

    def __init__(self, array):
        if not self.is_unitary(array):
            raise ValueError("input-output relations must be unitary")
        self.array = array

    @classmethod
    def with_reflectivity(cls, r):
        array = np.array([[math.sqrt(1 - r**2), r], [-r, math.sqrt(1 - r**2)]])
        return cls(array)

    @staticmethod
    def is_unitary(array):
        return np.allclose(np.eye(array.shape[0]), array.H * array)

    def __call__(self, state):
        return self.evolve(state)

    def evolve(self, initial_state):
        final_state = state.DensityMatrix.zero()
        for c, ketbra in initial_state.ketbras():
            if not c == 0:
                final_state += c*self.evolve_ketbra(ketbra)
        return final_state

    def evolve_ketbra(self, ketbra):
        ket, bra = ketbra
        return state.DensityMatrix.ketbra(
            self.evolve_number_ket(ket),
            self.evolve_number_bra(bra)
            )

    def evolve_number_ket(self, ket):
        final_ket = state.State.vacuum()
        self._evolve_number_ket(initial_ket, final_ket)
        return final_ket

    def _evolve_number_ket(self, initial_state, final_state):
        N1, N2, N3, N4 = initial_state
        if not N1 == 0:
            down = operator.Down(1, 4)
            up = operator.Up(3, 4) + operator.Up(4, 4)
        elif not N2 == 0:
            down = operator.Down(2, 4)
            up = operator.Up(3, 4) + operator.Up(4, 4)
        elif not N3 == 0:
            down = operator.Down(3, 4)
            up = operator.Up(1, 4) + operator.Up(2, 4)
        elif not N4 == 0:
            down = operator.Down(4, 4)
            up = operator.Up(1, 4) + operator.Up(2, 4)
        else:
            return final_state
        return self._evolve_Neigenstate(down*initial_state, up*final_state)

    def evolve_number_bra(self, initial_state):
        return self.evolve_number_ket(initial_state.dagger()).dagger()
