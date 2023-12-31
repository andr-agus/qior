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

import qutip as qp
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

import qior

from qior.example.utils import *

def main():
    """
    Example 1 provides the user insight on the meaning of the qior 
    parameter labeled as "reflectivity". This function creates an two-figure
    plot.

    On the left subplot there is an histogram-like representation of the
    density matrix that represents the same initial state as in example 0.

    On the right subplot there is an animation of the output state for
    different values of the reflectivity. It is straightforward to see that
    with no reflectivity, the excitation in the first mode is transmitted to
    the second one, while with unit reflectivity the excitation stays in the
    first mode. For intermediate values, a coherent superposition of 
    transmitted and reflected excitations is produced.
    """
    # set the initial state, the same as in example 0
    d1, d2 = dims = (2, 2)
    psi = qp.tensor(qp.basis(d1, 1), qp.basis(d2, 0))
    input = psi*psi.dag()

    # create a tuple containing different values of the reflectivity
    Rs = linspace(start = 0, end = 1, steps = 20)

    # make an empty figure with two subplots
    fig = prepare_input_output_figure()

    # make a function that plots the tomography of the input and output
    # density matrix on each subplot of that figure
    def frame_update(frame_index):
        r = Rs[frame_index]
        dims = tuple(input.dims[0])
        relation = qior.with_reflectivity(r, dims)
        output = relation(input)
        plot_input_tomography(input, fig)
        plot_output_tomography(output, fig, r)

    # make an animation using that function and all the values of Rs
    # (first, some optional arguments for finetuning the animation)
    kw = dict(frames = len(Rs), interval = 2000/len(Rs))
    # Now run the animation
    animation = FuncAnimation(fig, frame_update, **kw)
    # and save it
    animation.save("qior-example1.mp4")

if __name__ == "__main__":
    main()
