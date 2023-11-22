"""
This module instantiates the class InputOutputRelation, as well as a 
function "with_reflectivity" that conveniently creates a relation from a 
single float representing the reflectivity
"""

import math

import numpy as np
import qutip as qp

def with_reflectivity(*a, **kw):
    """
    See InputOutputRelation.with_reflectivity docstring
    """
    return InputOutputRelation.with_reflectivity(*a, **kw)

class InputOutputRelation:
    """
    Input-output relations are often written as formulas equating a creation
    operator to a weighted sum of other creation operators. It is often
    consfusing interpreting those equalities as such. What they mean is
    that the time-evolution of a physical system has been discretized in time
    so that an initial state becomes, after a particular time step, a final
    state according to the following rule:

        - Write the initial state as the action of a bunch of input creation 
          operators acting on vacuum.
        - Substitute each of those creation operators by the weighted sum
          of the output creation operators these input-output relations set
          them "equal" to.
        - The actions of the new output creation operators on vacuum returns
          the final state

    The rule looks simple, but understanding the properties of the final 
    states analitically is, often, non-trivial. Therefore, this class allows 
    its user to create arbitrary linear and reversible input-output relations
    from numpy arrays, as well as to compute the final state when calling
    this object with the initial state as argument.
    """

    some_instances = dict()

    def __init__(self, matrix, dims, acting_on = (0,1)):
        """
        Initialize an input-output relation.

        Arguments:
            - matrix: a numpy array with shape (2, 2) containing the
              coefficients of the analytical input-output relation desired.
              To fix notation, the matrix argument is considered to follow

                    a_i1 -> matrix[0,0] a_o1 + matrix[0, 1] a_o2
                    a_i2 -> matrix[1,0] a_o1 + matrix[1, 1] a_o2

              where a_i1 and a_i2 are the input modes anihilation operators
              while a_o1 and a_o2 the output modes. Computing the final
              state resulting from the application of this relation is done
              as indicated in InputOutputRelation.__doc__. Note that for
              that rule to prescribe a unitary evolution in the states, the
              argument "matrix" has to be a unitary matrix too!

            - acting_on: an iterable of two integers indicating the indices
              of the two systems this relation acts on in the order prescribed
              by the "dims" argument.

            - dims: input-output relations describe the dynamics of quantum
              states that belong to infinite dimensional Hilbert spaces.
              Moreover, these spaces are, at least, bipartite. The argument
              "dims" is an iterable of integers indicating the maximum number
              of excitations to be considered in each single-body space, so
              that their dimension is not infinite. The input-output relation
              acts on only two of those single-body spaces, while acting as
              the identity on the rest. Those two single-body spaces start
              as the input modes and end up representing the output modes.
              Note that input-output relations represent passive processes
              where energy is conserved, so considering a bigger output state
              space is not necessary since energy is not going to be created.
        """

        if not self.is_unitary(matrix):
            raise ValueError("input-output relations must be unitary")

        if len(dims) < 2:
            raise ValueError("input-output relations act on multipartite systems!")

        if not len(acting_on) == 2:
            raise ValueError("input-output relations act on two systems, not less, not more")

        if acting_on[0] < 0 or len(dims) <= acting_on[0]:
            raise ValueError("input-output relations must act on systems whose state dimension is provided with the 'dims' argument, not outside its indices")

        if acting_on[1] < 0 or len(dims) <= acting_on[1]:
            raise ValueError("input-output relations must act on systems whose state dimension is provided with the 'dims' argument, not outside its indices")

        self.matrix = matrix
        self.dims = dims
        self.acting_on = acting_on
        self.U = self.time_evolution()

    @staticmethod
    def is_unitary(array):
        """ Return True iff matrix is close to a numpy unitary"""
        identity = np.eye(array.shape[0])
        UdU = np.matmul(array.conj().T, array)
        return np.allclose(identity, UdU)

    def time_evolution(self):
        """
        Return a qp.Qobj representing the unitary evolution that turns
        any initial state into a final state.
        """
        return self.expand_to_bigger_system(self.local_time_evolution())

    def local_time_evolution(self):
        """
        Return a qp.Qobj representing the unitary evolution that turns a 
        reduced initial state of the input modes alone into a local final
        state.
        """
        U = self.local_zero_operator()
        d1, d2 = self.local_dims()
        for n1 in range(d1):
            for n2 in range(d2):
                bra = qp.tensor(qp.basis(d1, n1), qp.basis(d2, n2)).dag()
                ket = self.evolve_photon_numbers(n1, n2)
                U += ket * bra # map to the ket |n1n2> the ket "ket"
        return U

    def local_zero_operator(self):
        """
        Return the zero operator that has dimensions matching to the input modes
        """
        identities = tuple(map(qp.qeye, self.local_dims()))
        return 0 * qp.tensor(*identities)

    def local_dims(self):
        """
        Return a tuple of two integers indicating the cutoffs regarding the input-output modes
        """
        return self.dims[self.acting_on[0]], self.dims[self.acting_on[1]]

    def evolve_photon_numbers(self, n1, n2):
        """
        Return the reduced state containing only the output modes' state 
        resulting from the evolution of an initial state with n1 photons
        on the first input mode and n2 photons on the second input mode.
        """
        if n1 == 0 and n2 == 0:
            return self.local_vacuum() # vacuum maps to vacuum
        elif not n1 == 0:
            output_operator = self.output1() / math.sqrt(n1)
            return output_operator * self.evolve_photon_numbers(n1 - 1, n2)
        elif not n2 == 0:
            output_operator = self.output2() / math.sqrt(n2)
            return output_operator * self.evolve_photon_numbers(n1, n2 - 1)

    def local_vacuum(self):
        d1, d2 = self.local_dims()
        return qp.tensor(qp.basis(d1), qp.basis(d2))

    def output1(self):
        """
        Return the creation operator of the first output mode as a weighted
        sum of the creation operators of the input modes.

        Remember that this operator is not supposed to create a photon from
        vacuum in the output modes, but to create a photon from a
        superposition of photons in the input modes that end up as a single
        photon on this output mode.
        """
        return self.matrix[0, 0] * self.input1() + \
               self.matrix[0, 1] * self.input2()

    def output2(self):
        """
        Return the creation operator of the second output mode as a wighted
        sum of the creation operators fothe input modes.

        See output1.__doc__ for more information.
        """
        return self.matrix[1, 0] * self.input1() + \
               self.matrix[1, 1] * self.input2()

    def input1(self):
        """
        Return the creation operator on the first input mode.
        """
        d1, d2 = self.local_dims()
        return qp.tensor(qp.create(d1), qp.qeye(d2))

    def input2(self):
        """
        Return the creation operator on the second input mode.
        """
        d1, d2 = self.local_dims()
        return qp.tensor(qp.qeye(d1), qp.create(d2))

    def expand_to_bigger_system(self, U):
        """
        Take a qutip.Qobj operator defined on the input-output modes alone
        and tensor-wrap it with identities so that it becomes a global
        operator.
        """
        new_dims = self.permute_dims()

        for D in new_dims[2:]:
            U = qp.tensor(U, qp.qeye(D))

        return self.permute_back_unitary(U)

    def permute_dims(self):
        """
        Return a rearrangement of self.dims so that the local dimensions
        are in the two first entries.
        """
        new_dims = list(self.dims)
        new_dims[0] = self.dims[self.acting_on[0]]
        new_dims[self.acting_on[0]] = self.dims[0]
        new_new_dims = list(new_dims)
        new_new_dims[1] = new_dims[self.acting_on[1]]
        new_new_dims[self.acting_on[1]] = new_dims[1]
        return new_new_dims

    def permute_back_unitary(self, U):
        """
        Return U with the tensor order permuted so that the two first
        subsystems end up at indices self.acting_on[0] and self.acting_on[1]
        """
        identity = list(range(len(self.dims)))
        first_permutation = list(identity)
        first_permutation[1] = identity[self.acting_on[1]]
        first_permutation[self.acting_on[1]] = identity[1]
        second_permutation = list(first_permutation)
        second_permutation[0] = first_permutation[self.acting_on[0]]
        second_permutation[self.acting_on[0]] = first_permutation[0]
        return U.permute(second_permutation)


    @classmethod
    def with_reflectivity(cls, R, dims, acting_on = (0, 1)):
        """
        A typical parameter used to enunciate input-output relations is 
        reflectivity. This method allows the user to create relations with
        that parameter instead of having to write down the complete matrix.

        However, there are a couple of conventions to fix in order to
        properly define reflectivity. Everyone agrees that the process of 
        reflection is dual to that of transmission. If we have two systems
        whose global evolution is to be defined by an input-output relation,
        we define reflectivity as follows:

        1. Consider an initial state in which the first system contains an
        excitation while the second one contains none.

        2. Because of the basic properties of input-output relations, the
        final state will have the same global number of excitations as the
        initial one, just not in the same subsystems.

        3. Transmission takes place if the excitation in the first system
        leaks into the second. A complete Rabi oscillation is an input
        output relation with transmitivity of 1

        4. Reflection takes place if the excitation in the first system
        stays in that system. Two decoupled systems evolve under an input-
        output relation with reflectivity 1

        5. Reflectivity is considered a probability of not-transmission, 
        therefore it is a real number between 0 and 1, both included. This
        definition takes from the user the ability to change some relative
        phases, but this is the approach that makes the most sense to 
        the inmediate applications of this program.

        We understand that the user may want to define reflectivity the other
        way around: the probability of getting an initial excitation on one
        system to the other system. The user can write a wrapper around this
        method to change its behaviour that way.

        With these statements, it is safe to define the argument "R" as the
        reflectivity. The arguments "dims" and "acting_on" have the same 
        meaning as in cls.__init__
        """
        key = (R, dims, acting_on)
        self = cls.some_instances.get(key, None)
        if self is None:
            array = np.matrix([[math.sqrt(R), math.sqrt(1-R)],
                               [math.sqrt(1-R), -math.sqrt(R)]])
            return cls(array, dims, acting_on)
        else:
            return self

    def __call__(self, state):
        """
        Return the final state resulting of applying "self" to "state"
        """
        return self.evolve(state)

    def evolve(self, initial_state):
        """
        Apply the unitary matrix computed in self.time_evolution_operator()
        to an initial_state and return the final state. The initial state must
        have a null projection in some of the high energy eigenstates below
        the cutoff so that the output state is sure to be contained below the
        cutoff. Specifically,

        1. The cutoff of the first input mode is 
                
                self.dims[self.acting_on[0]]

           or N1_max for short. Define N2_max in an analogous manner.

        2. If the initial state has a non-zero projection onto both the 
           eigenstates with N1_max and N2_max eigenvalues then, depending on
           reflectivity, the final state may leak outside of the cutoff if
           those excitations add up into the first or second output modes up
           to photon numbers of N_1max + N2_max

        3. To avoid this, this method checks which of the projections of the
           initial state onto the number states are null. If all the number
           state with non-zero projection and maximum number of photons is n1
           for the first mode and n2 for the second mode, these inequalities
           must hold
                
                n1 + n2 <= N1_max        and       n1 + n2 <= N2_max

           so that the output state does not evolve outside of the cutoff. If
           it does leak outside the cutoff, and exception is thrown
        """
        if self.output_leaks_outside_dims(initial_state):
            dims = initial_state.dims[0]
            raise ValueError(("given the input output relation %s and its" + \
            " cutoffs %s, the output state is not contained within those " + \
            "cutoffs") % (self, dims))
        if self.is_pure(initial_state):
            return self.U * initial_state
        else:
            return self.U * initial_state * self.U.dag()

    def output_leaks_outside_dims(self, initial_state):
        """
        Return True iff the initial_state would result in a final state
        that leaks ouside of the dimensions specified in self.dims
        """
        N0 = self.max_number_of_photons(initial_state, self.acting_on[0])
        N1 = self.max_number_of_photons(initial_state, self.acting_on[1])
        if N0 + N1 >= self.dims[self.acting_on[0]]:
            return True
        if N0 + N1 >= self.dims[self.acting_on[1]]:
            return True
        return False

    def max_number_of_photons(self, state, subsystem):
        """
        Return the photon number of the photon state that has a non-zero
        projection onto "state" with highest photon number.
        """
        d = self.dims[subsystem]
        nmax = d - 1
        while not nmax == 0:
            proj = self.projector_on_photon_number_eigenspace(state, subsystem, nmax)
            exp_val = qp.expect(proj, state)
            if exp_val:
                return nmax
            nmax -= 1
        return nmax

    @staticmethod
    def projector_on_photon_number_eigenspace(state, subsystem, eigenvalue):
        dims = state.dims[0]
        opers = []
        for i, d in enumerate(dims):
            if i == subsystem:
                vec = qp.basis(d, eigenvalue)
                proj = vec * vec.dag()
                opers.append(proj)
            else:
                opers.append(qp.qeye(d))
        return qp.tensor(*opers)

    @classmethod
    def is_pure(cls, state):
        """
        Return True iff the state is pure
        """
        return not state.dims[0] == state.dims[1]
