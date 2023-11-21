"""quinoa: QUantum INput Output relAtions

Quinoa is a Python package for the evalutation of the evolution of arbitrary
initial states under quantum input-output relations. The user interface of
this module is:
    - InputOutputRelation: Class representing an input-output relation
    - InputOutputRelation.__call__: The objects of this class can be used
      as functions that take as input the input state and return the output
      state
    - with_reflectivity: Convenience builder for a particular kind of input-
      output relation where phases are not important but only the magnitude
      of the reflectivity.

For more information see the doc strings of those objects as well as the
examples provided by running this module with "python -m quinoa"
"""

from .relation import InputOutputRelation, with_reflectivity
