"""
Test use of FunctionExpression bindings
"""
import pathlib
import sys
import unittest
import numpy as np
from .tools import debug_on
import matplotlib.pyplot as plt
import time
import functools
import quadpy

import sympy as sp

np.random.seed()

try:
    from gsplines.collocation import legendre_gauss_lobatto_points_weights, \
        q_and_evaluation
except ImportError:
    MOD_PATH = pathlib.Path(__file__).parent.absolute()
    MOD_PATH_PYGSPLINES = pathlib.Path(MOD_PATH, '..', 'build')
    sys.path.append(str(MOD_PATH_PYGSPLINES))
    from gsplines.collocation import legendre_gauss_lobatto_points_weights, \
        q_and_evaluation


class MyTest(unittest.TestCase):
    """ Test elemental function and operations"""
    @debug_on()
    def __init__(self, *args, **kwargs):
        unittest.TestCase.__init__(self, *args, **kwargs)
        np.random.seed()

    @ debug_on()
    def pol_test(self):

        for pol_degre in range(2, 30):
            scheme = quadpy.c1.gauss_legendre(pol_degre)
            points_nom = scheme.points
            for root in points_nom:
                ln_test = q_and_evaluation(root, pol_degre)[2]
                self.assertLessEqual(np.abs(ln_test), 1.0e-9)

        for pol_degre in range(4, 30):
            scheme = quadpy.c1.gauss_lobatto(pol_degre)
            points_nom = scheme.points
            for root in points_nom:
                ln_test = q_and_evaluation(root, pol_degre-1)[0]
                self.assertLessEqual(np.abs(ln_test), 1.0e-9)

    @ debug_on()
    def points_test(self):

        for number_of_points in range(3, 30):
            scheme = quadpy.c1.gauss_lobatto(number_of_points)

            points_nom = scheme.points

            points_test, _ = \
                legendre_gauss_lobatto_points_weights(number_of_points)

            self.assertLessEqual(np.linalg.norm(
                points_nom - points_test.T), 1.0e-9)

    @ debug_on()
    def weights_test(self):

        for number_of_points in range(2, 30):
            scheme = quadpy.c1.gauss_lobatto(number_of_points)

            weights_nom = scheme.weights

            _, weights_test = \
                legendre_gauss_lobatto_points_weights(number_of_points)

            self.assertLessEqual(np.linalg.norm(
                weights_nom - weights_test.T), 1.0e-9)

    @ debug_on()
    def test(self):

        self.pol_test()

        self.points_test()
        self.weights_test()


if __name__ == '__main__':
    unittest.main()
