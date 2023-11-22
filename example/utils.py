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

def plot_input_tomography(state, fig, colorbar = False):
    if is_pure(state):
        state = state*state.dag()
    ax_in = fig.axes[0]
    ax_in.clear()
    xlabels = labels_for_state(state)
    hist_in = qp.visualization.matrix_histogram_complex(
        state, xlabels, xlabels, fig = fig, ax = ax_in,
        colorbar = colorbar,
        )
    ax_in.set_title("input state")
    ax_in.set_xlabel("photon number")
    ax_in.set_ylabel("photon number")
    ax_in.set_zlim((0, 1))
    return fig

def is_pure(state):
    return state.type == "ket"

def plot_output_tomography(state, fig, reflectivity, colorbar = False):
    if is_pure(state):
        state = state*state.dag()
    xlabels = labels_for_state(state)
    ax_out = fig.axes[1]
    ax_out.clear()
    hist_out = qp.visualization.matrix_histogram_complex(
        state, xlabels, xlabels, fig = fig, ax = ax_out,
        colorbar = colorbar,
        )
    if isinstance(reflectivity, float):
        ax_out.set_title("output when R = %f" % reflectivity)
    else:
        ax_out.set_title("output when R = %s" % reflectivity)
    ax_out.set_xlabel("photon number")
    ax_out.set_ylabel("photon number")
    ax_out.set_zlim((0, 1))

def labels_for_state(state):
    labels = []
    for ints in itertools.product(*[range(D) for D in state.dims[0]]):
        labels.append("".join(map(str, ints)))
    return labels
    
def plot_fock_distribution(state, fig, ax = None):
    if ax is None:
        ax = fig.gca()
    xlabels = labels_for_state(state)
    qp.plot_fock_distribution(state, fig = fig, ax = ax)
    ax.set_xticks(range(state.shape[0]), xlabels)
    ax.set_xlabel("Fock numbers")
