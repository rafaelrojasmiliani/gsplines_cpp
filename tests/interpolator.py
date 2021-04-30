import pathlib
import os
import sys
import unittest
import numpy as np
from .tools import debug_on
from numpy.polynomial.legendre import Legendre
import matplotlib.pyplot as plt
try:
    from pygsplines import BasisLegendre
    from pygsplines import PyInterpolator as Interpolator
except ImportError:
    MOD_PATH = pathlib.Path(__file__).parent.absolute()
    MOD_PATH_PYGSPLINES = pathlib.Path(MOD_PATH, '..', 'build')
    sys.path.append(str(MOD_PATH_PYGSPLINES))
    from pygsplines import BasisLegendre
    from pygsplines import PyInterpolator as Interpolator


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
    plt.show()


class MyTest(unittest.TestCase):
    @debug_on()
    def __init__(self, *args, **kwargs):
        super(MyTest, self).__init__(*args, **kwargs)

    @debug_on()
    def test(self):
        basis = BasisLegendre(6)
        dim = 2
        intervals = 2
        wp = np.random.rand(intervals+1, dim)
        tau = np.array(intervals*[1])
        inter = Interpolator(dim, intervals, basis)
        res = inter.interpolate(tau, wp)
        inter.print_interpolating_matrix()
        inter.print_interpolating_vector()

        show_piecewisefunction(res, 5, 0.001)


if __name__ == '__main__':
    unittest.main()
