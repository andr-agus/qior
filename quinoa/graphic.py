from matplotlib import pyplot as plt
import qutip

def fock_distribution(initial_state, final_state):
    fig, ax = plt.subplots()
    qutip.plot_fock_distribution(final_state, fig = fig, ax = ax)
    ticks = range(final_state.shape[0])
    labels = []
    for i in range(final_state.dims[0][0]):
        for j in range(final_state.dims[0][1]):
            labels.append("|" + str(i) + " " + str(j) + ">")
    ax.set_xticks(ticks, labels)
    for label in ax.get_xticklabels():
        label.set(rotation=90)
    return fig, ax
