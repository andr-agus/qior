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
import qior
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

from qior.example.utils import *

def main():
    """
    All the examples so far began with no photons on the second input mode.
    However, the second input mode has a relevant difference when compared to
    the first one. The photons that start in that second input mode and are
    reflected back to the second output mode gain a -1 phase. If they were to
    behave like photons on the first mode and gain no phase the global
    evolution would not be unitary.

    This example shows what happens when the input state contains a photon on
    each the first and the second modes. That -1 phase obtained in the
    reflection in the second mode is responsible for an interesting quantum
    interference phenomenon known as the Hong Ou Mandel effect (Phys. Rev.
    Lett. 59 2044 (1987) where the output state consists of a coherent 
    superposition of having both initial photons in the first output mode or
    both in the second output mode, but no state with only one photon per 
    mode.

    To showcase these intereference phenomena this function creates an
    animation for different values of the reflectivity as well as a static
    figure for reflectivity 0.5, where the HOM effect takes place.
    """
    # set the initial parameters
    d1, d2 = dims = (3, 3)
    psi = qp.tensor(qp.basis(d1, 1), qp.basis(d2, 1))
    input = psi*psi.dag()
    Rs = linspace(start = 0, end = 1, steps = 30)

    # prepare the figure
    fig = prepare_input_output_figure()

    # create a function that, in turn, creates each frame for the animation
    def frame_update(frame):
        r = Rs[frame]
        relation = qior.with_reflectivity(r, dims, acting_on = (0, 1))
        output = relation(input)
        plot_input_tomography(input, fig, colorbar = frame_update.colorbar)
        plot_output_tomography(output, fig, r, colorbar = frame_update.colorbar)
        frame_update.colorbar = False
    frame_update.colorbar = True

    # run the animation and save it to disk
    options = dict(frames = len(Rs), interval = 3000/len(Rs))
    animation = FuncAnimation(fig, frame_update, **options)
    animation.save("qior-example4.mp4")

    # now it is time for the static figure
    r = 0.5
    relation = qior.with_reflectivity(r, dims)
    output = relation(input)
    fig = prepare_input_output_figure()
    plot_input_tomography(input, fig, colorbar = True)
    plot_output_tomography(output, fig, r, colorbar = True)
    fig.savefig("qior-example4.png")


if __name__ == "__main__":
    main()
