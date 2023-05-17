import math

import qutip
from matplotlib import pyplot as plt

from quinoa import *

def main():
    r = relation.IORelation.with_reflectivity(.5, (6, 6))
    initial_state = qutip.tensor(qutip.basis(6), qutip.basis(6,5))
    final_state = r.evolve(initial_state)

    graphic.fock_distribution(initial_state, final_state)

    plt.show()

if __name__ == "__main__":
    main()
