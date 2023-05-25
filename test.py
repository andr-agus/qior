import math

import qutip
from matplotlib import pyplot as plt

from quinoa import *
import quinoa

def three_spdc(beta, cutoffs):
    alpha = math.sqrt(1 - abs(beta)**2)
    pure = alpha * qutip.tensor(
            qutip.basis(cutoffs[0]),
            qutip.basis(cutoffs[1]),
            qutip.basis(cutoffs[2]),
        ) + beta * qutip.tensor(
            qutip.basis(cutoffs[0], 1),
            qutip.basis(cutoffs[1], 1),
            qutip.basis(cutoffs[2], 1),
        )
    return pure*pure.dag()

def initial_radar_state(beta, nth, cutoffs):
    signal = three_spdc(beta, cutoffs[:3])
    background = qutip.tensor(
        qutip.thermal_dm(cutoffs[3], nth),
        qutip.thermal_dm(cutoffs[4], nth),
        )
    return qutip.tensor(signal, background)

def state_with_no_target(nth, cutoffs):
    r1 = relation.IORelation.with_reflectivity(0, cutoffs, acting_on = (1,3))
    r2 = relation.IORelation.with_reflectivity(0, cutoffs, acting_on = (2, 4))
    initial_state = initial_radar_state(0, nth, cutoffs)
    final_global_state = r2(r1(initial_state))
    final_local_state = final_global_state.ptrace((0,3,4))
    return final_local_state

def state_with_target(signal_strength, nth, beam_splitter1, beam_splitter2, cutoffs):
    initial_state = initial_radar_state(signal_strength, nth, cutoffs)
    return beam_splitter2(beam_splitter1(initial_state)).ptrace((0,3,4))

def bounds_for_many_parameters():
    signal_strengths = [0.1*i for i in range(11)]
    reflectivities = [0.1*i for i in range(11)]
    nths = [0.5*i for i in range(11)]
    cutoffs = (2, 4, 4, 4, 4)

    rho0s = [state_with_no_target(nth, cutoffs) for nth in nths]

    for r in reflectivities:
        bs1 = relation.IORelation(r, cutoffs, acting_on = (1, 3))
        bs2 = relation.IORelation(r, cutoffs, acting_on = (2, 4))
        for rho0, nth in zip(rho0s, nths):
            for s in signal_strengths:
                rho1 = state_with_target(s, nth, bs1, bs2, cutoffs)
                bound = qinfo.bhattacharyya_bound(rho0, rho1, 1)

def plot_tomography(rho, title):
    fig, ax = quinoa.graphic.tomography(rho)
    fig.suptitle(title)
    return fig, ax

def compute_bounds(beta, nth, Rs):
    bounds = []
    cutoffs = (2, 4, 4, 4, 4)
    plot_tomography(three_spdc(beta, cutoffs[:3]), "signal before target")
    rho0 = state_with_no_target(nth, cutoffs)
    plot_tomography(rho0.ptrace((1,)), "no target - mode 1")
    plot_tomography(rho0.ptrace((2,)), "no target - mode 2")
    plot_tomography(rho0, "no target - all modes")
    for R in Rs:
        if R == Rs[-1]:
            plot_tomography(rho1, "target - all modes")
        bs1 = relation.IORelation.with_reflectivity(R, cutoffs, acting_on = (1, 3))
        bs2 = relation.IORelation.with_reflectivity(R, cutoffs, acting_on = (2, 4))
        rho1 = state_with_target(beta, nth, bs1, bs2, cutoffs)
        bounds.append(qinfo.bhattacharyya_bound(rho0, rho1, 1))
    return bounds

def main():
    beta = 0.25
    nth = 100
    Rs = [0.01*i for i in range(11)]
    bounds = compute_bounds(beta, nth, Rs)
    fig, ax = plt.subplots()
    ax.plot(Rs, bounds)
    ax.set_xlabel("reflectivity")
    ax.set_ylabel("upper bound on error probability")
    plt.show()

if __name__ == "__main__":
    main()
