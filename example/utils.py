"""
"""
import itertools

import qutip as qp
from matplotlib import pyplot as plt

def linspace(start, end, steps):
    delta = (end - start)/steps
    return tuple([start + i*delta for i in range(steps+1)])

def prepare_input_output_figure():
    fig, _ = plt.subplots(1, 2,
        figsize = (20, 10),
        subplot_kw = {"projection":"3d"},
    )
    fig.subplots_adjust(left=.01, right = .99, wspace=.1)
    return fig

def plot_input_tomography(state, fig):
    ax_in = fig.axes[0]
    ax_in.clear()
    xlabels = labels_for_state(state)
    hist_in = qp.visualization.matrix_histogram_complex(
        state, xlabels, xlabels, fig = fig, ax = ax_in,
        colorbar = False,
        )
    ax_in.set_title("input state")
    ax_in.set_xlabel("photon number")
    ax_in.set_ylabel("photon number")
    ax_in.set_zlim((0, 1))
    return fig

def plot_output_tomography(state, fig, reflectivity):
    xlabels = labels_for_state(state)
    ax_out = fig.axes[1]
    ax_out.clear()
    hist_out = qp.visualization.matrix_histogram_complex(
        state, xlabels, xlabels, fig = fig, ax = ax_out,
        colorbar = False,
        )
    ax_out.set_title("output when R = %f" % reflectivity)
    ax_out.set_xlabel("photon number")
    ax_out.set_ylabel("photon number")
    ax_out.set_zlim((0, 1))

def labels_for_state(state):
    labels = []
    for ints in itertools.product(*[range(D) for D in state.dims[0]]):
        labels.append("".join(map(str, ints)))
    return labels
    
def plot_compu_basis_projections(state, fig):
    xlabels = labels_for_state(state)
    ax = fig.axes[0]

def expect_computational_basis(state):
    """
    Return an iterable containing the expectation values of the projectors on
    each of the elements of the computational basis. This is the same as
    the square of the wavefunction if the state is pure.
    """
    return [qp.expect(vec*vec.dag(), state) for vec in computational_basis(state)]

def computational_basis(state):
    """
    Return an iterable containing all the states of the computational basis
    that have the same cutoffs as the argument state.
    """
    dims = state.dims[0]
    vecs = []
    for ints in itertools.product(*[range(d) for d in dims]):
        vecs.append([])
        for d, excitation in zip(dims, ints):
            vecs[-1].append(qp.basis(d, excitation))
        vecs[-1] = qp.tensor(*vecs[-1])
    return vecs
