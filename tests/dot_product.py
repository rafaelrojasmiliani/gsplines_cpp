import pathlib
import os
import sys
import unittest
import numpy as np
from .tools import debug_on
from numpy.polynomial.legendre import Legendre
import matplotlib.pyplot as plt
import time
try:
    from gsplines.optimization import minimum_jerk_path
except ImportError:
    MOD_PATH = pathlib.Path(__file__).parent.absolute()
    MOD_PATH_PYGSPLINES = pathlib.Path(MOD_PATH, '..', 'build')
    sys.path.append(str(MOD_PATH_PYGSPLINES))
    from gsplines.basis import BasisLegendre
    from gsplines.optimization import minimum_jerk_path


class MyTest(unittest.TestCase):
    # @debug_on()
    def __init__(self, *args, **kwargs):
        super(MyTest, self).__init__(*args, **kwargs)
        np.random.seed()

    # @debug_on()
    def test(self):
        intervals = 3
        codom_dim = 3
        delta_time = 0.001
        time_spam = np.arange(0, 1+delta_time, delta_time)
        for _ in range(5):
            waypoints = np.random.rand(intervals+1, codom_dim)*6.14
            function_1 = minimum_jerk_path(waypoints)
            waypoints = np.random.rand(intervals+1, codom_dim)*6.14
            function_2 = minimum_jerk_path(waypoints)

            test_function = function_1.dot(function_2)

            test_result = test_function(time_spam)

            nom_result = np.einsum(
                'ij,ij->i', function_1(time_spam),
                function_2(time_spam)).reshape(-1, 1)

            test_value = np.max(np.abs(test_result-nom_result))

            self.assertLessEqual(test_value, 1.0e-9)

            test_function = test_function.deriv()

            test_result = test_function(time_spam)

            nom_result = np.einsum(
                'ij,ij->i', function_1(time_spam),
                function_2.deriv()(time_spam)).reshape(-1, 1) + \
                np.einsum(
                'ij,ij->i', function_1.deriv()(time_spam),
                function_2(time_spam)).reshape(-1, 1)

            test_value = np.max(np.abs(test_result-nom_result))

            self.assertLessEqual(test_value, 1.0e-9)
        # for q in q_array:
        #    show_piecewisefunction(q, 5, 0.001)
        # plt.show()


if __name__ == '__main__':
    unittest.main()
