"""
"""

import qutip as qp
import quinoa

def main():
    d1, d2 = dims = (2, 2)
    input = qp.tensor(qp.basis(d1, 1), qp.basis(d2, 1))
    relation = quinoa.with_reflectivity(.5, dims)
    output = relation(input)
    print(input)
    print(output)


if __name__ == "__main__":
    main()
