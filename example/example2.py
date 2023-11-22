"""
"""

import qutip as qp

import quinoa
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

from utils import *

def main():
    """
    Example 2 shows that quinoa deals with multiphoton states consistently.

    Now, the input state contains three photons on the first mode. Then, the
    tranmission and reflection of those photons result in a final state that
    has non-zero projections onto |3, 0>, |2, 1>, |1, 2> and |0 ,3> with
    coefficients that follow the binomial distribution with the probability 
    of reflecting a single photon being the reflectivity.

    This function writes an animation equivalent to example 1 to disk, but
    since the photon number is larger than one, the animation contains too
    much information to the user. Thus, we save a static figure to disk for
    the specific value of reflectivity = 0.5 too. That figure shows that the
    probabilities of measuring each of the number eigenstates, with no 
    information regarding coherences or phases, see example 4 for more info.
    """
    # set the initial parameters
    d1, d2 = dims = (4, 4)
    psi = qp.tensor(qp.basis(d1, 3), qp.basis(d2, 0))
    input = psi*psi.dag()
    Rs = linspace(start = 0, end = 1, steps = 20)

    # prepare the figure
    fig = prepare_input_output_figure()

    # create a function that, in turn, creates each frame for the animation
    def frame_update(frame):
        r = Rs[frame]
        relation = quinoa.with_reflectivity(r, dims, acting_on = (0, 1))
        output = relation(input)
        plot_input_tomography(input, fig)
        plot_output_tomography(output, fig, r)

    # run the animation and save it to disk
    kw = dict(frames = len(Rs), interval = 2000/len(Rs))
    animation = FuncAnimation(fig, frame_update, **kw)
    animation.save("quinoa-example2.mp4")

    # now it is time for the static figure
    r = 0.5
    relation = quinoa.with_reflectivity(r, dims)
    output = relation(input)
    fig, _ = plt.subplots()
    plot_fock_distribution(output, fig)
    fig.savefig("quinoa-example2.png")

if __name__ == "__main__":
    main()
