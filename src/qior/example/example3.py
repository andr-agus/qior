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

def main():
    """
    Example 3 shows that the mode-by-mode cutoffs allow for some input modes
    that produce output modes that leak outside those cutoffs. That is, if
    we have two modes with dimension 2, that is, cutoff to have at most one
    photon each, then the initial state that has a photon on each mode leaks
    outside the cutoff to states having two photons on one mode and so on.

    qior checks whether this happens and throws an exception if that is 
    the case.
    """
    d1, d2 = dims = (2, 2)
    input = qp.tensor(qp.basis(d1, 1), qp.basis(d2, 1))
    relation = qior.with_reflectivity(.5, dims)
    output = relation(input)
    print(input)
    print(output)


if __name__ == "__main__":
    main()
