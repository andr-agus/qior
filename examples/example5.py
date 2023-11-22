"""
"""
import math

import numpy
import qutip as qp
import quinoa
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

from utils import *

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
    relation = quinoa.InputOutputRelation(array, dims)
    output = relation(input)

    # prepare the figure
    fig = prepare_input_output_figure()
    plot_input_tomography(input, fig, colorbar = True)
    plot_output_tomography(output, fig, "¿?", colorbar = True)
    fig.savefig("quinoa-example5.png")

if __name__ == "__main__":
    main()
