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


def show_piecewisefunction(_q, _up_to_deriv=3, _dt=0.1, _title=''):
    dim = _q.get_codom_dim()
    fig, ax = plt.subplots(_up_to_deriv + 1, dim)
    if dim == 1:
        ax = np.array([[ax[i]] for i in range(_up_to_deriv + 1)])
    if _title:
        fig.suptitle(_title)
    t = np.arange(0.0, _q.get_exec_time(), _dt)

    for i in range(0, _up_to_deriv + 1):
        q = _q.deriv(i)
        qt = q(t)
        for j in range(0, dim):
            ax[i, j].plot(t, qt[:, j])
            ax[i, j].grid()
            if i == 0:
                ax[i, j].set_title('coordinate {:d}'.format(j + 1), fontsize=8)

            if hasattr(_q, 'get_domain_breakpoints'):
                for ti in _q.get_domain_breakpoints():
                    ax[i, j].axvline(ti, alpha=0.1, color='red')

    plt.subplots_adjust(
        left=0.025,
        bottom=0.05,
        right=0.975,
        top=0.95,
        wspace=0.25,
        hspace=0.15)


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

        for q in q_array:
            show_piecewisefunction(q, 5, 0.001)
        plt.show()


if __name__ == '__main__':
    unittest.main()
