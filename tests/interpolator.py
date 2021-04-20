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


class cMyTest(unittest.TestCase):
    @debug_on()
    def __init__(self, *args, **kwargs):
        super(cMyTest, self).__init__(*args, **kwargs)

    @debug_on()
    def test(self):
        basis = BasisLegendre(6)
        dim = 3
        intervals = 1
        wp = np.random.rand(intervals+1, dim)
        tau = np.array([1])
        inter = Interpolator(dim, intervals, basis)
        res = inter.interpolate(tau, wp)

        time_span = np.arange(0, 1, 0.001)
        res = res(time_span)

        plt.plot(res[:, 0])
        plt.show()


if __name__ == '__main__':
    unittest.main()
