"""
This module contains an obviously incomplete set of functions for the
computation of quantities relevant in quantum information theory.
"""

import qutip

def bhattacharyya_bound(state_0, state_1, copies):
    return ((state_0.sqrtm()*state_1.sqrtm()).tr()**copies)/2

def chernoff_bound(state_0, state_1, copies):
    raise NotImplementedError("someone knows how to do this?")
