# Copyright 2023 and later, Andres Agusti Casado
# This file is part of the python package qior.
# qior is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any
# later version.
# qior is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
# You should have received a copy of the GNU General Public License along 
# with qior. If not, see <https://www.gnu.org/licenses/>. 
"""
"""
import math

import numpy
import qutip as qp
import qior
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

from qior.example.utils import *

def main():
    """
    This example introduces the most general way of creating input-output
    relations, that is, through a unitary numpy array that contains the
    coefficients of the input-output relation. In particular, we propose
    a relation that adds a pi phase to the transmitted photons from the
    first mode to the second one. This is a different definition of 
    reflectivity that is as natural as the one used in the examples before,
    so with this example we show the user how to change those conventions
    if they so desire.
    """
    d1, d2 = dims = (2, 2)
    input = qp.tensor(qp.basis(d1, 1), qp.basis(d2, 0))
    array = 1/math.sqrt(2)*numpy.array([[-1, 1], [1, 1]])
    relation = qior.InputOutputRelation(array, dims)
    output = relation(input)

    # prepare the figure
    fig = prepare_input_output_figure()
    plot_input_tomography(input, fig, colorbar = True)
    plot_output_tomography(output, fig, "¿?", colorbar = True)
    fig.savefig("qior-example5.png")

if __name__ == "__main__":
    main()
