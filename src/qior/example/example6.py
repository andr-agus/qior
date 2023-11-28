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
    This example shows that the calculations performed by this package are
    able to retain the coherences and correlations of the two input modes
    with other systems the input-output relations do not act upon. That is,
    if more systems are provided to the relations at creation with the
    "dims" argument being an iterable with more elements than two, then
    the time evolution operator is wrapped with identities acting on those
    extra systems. To indicate on which two systems the relations act upon,
    provide a tuple with the two system indices as the "acting_on" argument.
    """
    # set the initial parameters
    d1, d2, d3 = dims = (2, 2, 2)
    input = qp.tensor(qp.basis(d1, 1), qp.basis(d2, 0), qp.basis(d3, 0))
    Rs = linspace(start = 0, end = 1, steps = 20)

    # the input state is boring so far. We create some coherences between
    # the first and third mode with a half reflective beam splitter.
    relation1 = qior.with_reflectivity(.5, dims, acting_on = (0, 2))
    input = relation1(input) # the new input for future relations

    # prepare a figure in which to draw an animation
    fig = prepare_input_output_figure()

    # create a function that, in turn, creates each frame for the animation
    def frame_update(frame):
        r = Rs[frame]
        relation = qior.with_reflectivity(r, dims, acting_on = (0, 1))
        # this relation retains the coherences between the first and third
        # modes in the input state, but creates new ones with the second
        # mode
        output = relation(input)
        plot_input_tomography(input, fig, colorbar = frame_update.colorbar)
        plot_output_tomography(output, fig, r, colorbar = frame_update.colorbar)
        frame_update.colorbar = False
    frame_update.colorbar = True

    # run the animation and save it to disk
    options = dict(frames = len(Rs), interval = 3000/len(Rs))
    animation = FuncAnimation(fig, frame_update, **options)
    animation.save("qior-example6.mp4")

if __name__ == "__main__":
    main()
