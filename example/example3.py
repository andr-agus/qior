"""
"""

import qutip as qp
import quinoa

def main():
    """
    Example 3 shows that the mode-by-mode cutoffs allow for some input modes
    that produce output modes that leak outside those cutoffs. That is, if
    we have two modes with dimension 2, that is, cutoff to have at most one
    photon each, then the initial state that has a photon on each mode leaks
    outside the cutoff to states having two photons on one mode and so on.

    quinoa checks whether this happens and throws an exception if that is 
    the case.
    """
    d1, d2 = dims = (2, 2)
    input = qp.tensor(qp.basis(d1, 1), qp.basis(d2, 1))
    relation = quinoa.with_reflectivity(.5, dims)
    output = relation(input)
    print(input)
    print(output)


if __name__ == "__main__":
    main()
