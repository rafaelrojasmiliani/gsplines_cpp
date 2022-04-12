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
    from gsplines.basis import BasisLegendre
    from gsplines.optimization import optimal_sobolev_norm
    from gsplines.optimization import broken_lines_path
    from gsplines.optimization import minimum_acceleration_path
    from gsplines.optimization import minimum_jerk_path
    from gsplines.optimization import minimum_snap_path
    from gsplines.optimization import minimum_crackle_path
except ImportError:
    MOD_PATH = pathlib.Path(__file__).parent.absolute()
    MOD_PATH_PYGSPLINES = pathlib.Path(MOD_PATH, '..', 'build')
    sys.path.append(str(MOD_PATH_PYGSPLINES))
    from gsplines.basis import BasisLegendre
    from gsplines.optimization import optimal_sobolev_norm
    from gsplines.optimization import broken_lines_path
    from gsplines.optimization import minimum_acceleration_path
    from gsplines.optimization import minimum_jerk_path
    from gsplines.optimization import minimum_snap_path
    from gsplines.optimization import minimum_crackle_path


class MyTest(unittest.TestCase):
    # @debug_on()
    def __init__(self, *args, **kwargs):
        super(MyTest, self).__init__(*args, **kwargs)
        np.random.seed()

    # @debug_on()
    def test(self):
        basis = BasisLegendre(4)
        dim = 7  # np.random.randint(1, 10)
        intervals = 4
        waypoints = np.random.rand(intervals+1, dim)*6.14
        exec_time = intervals
        res = optimal_sobolev_norm(waypoints, basis, [(1, 1)], exec_time)
        q_array = [res, broken_lines_path(waypoints),
                   minimum_acceleration_path(waypoints),
                   minimum_jerk_path(waypoints)]


if __name__ == '__main__':
    unittest.main()
