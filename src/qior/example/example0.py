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
    Example 0 provides the user with the minimal working example of qior.
    The idea is that all the code the user needs to know is contained here,
    and they do not need to read the implementation of qior.

    This function creates an initial state in a truncated bipartite system
    that can be interpreted as a two-mode state. The state contains an
    excitation on the first mode and vacuum in the other one.

    Then, create an input output relation with zero reflectivity and compute
    the output state. Since there is no reflectivity, the initial excitation
    in the first mode is completely transmited to the second one.
    """
    # set the dimensions to work with
    d1, d2 = dims = (2, 2)
    # set the initial state with a single excitation on the first mode
    input = qp.tensor(qp.basis(d1, 1), qp.basis(d2, 0))
    # create the input-output relation with zero reflectivity and taking
    # into account the cutoff in the Hilbert space
    relation = qior.with_reflectivity(0, dims)
    # compute the final state
    output = relation(input)
    print("initial state: ", input)
    print("final state: ", output)
    # for more interesting examples see example1

if __name__ == "__main__":
    main()
