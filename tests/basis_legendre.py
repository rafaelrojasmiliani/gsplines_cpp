
import pathlib
import os
import sys
import unittest
import numpy as np
from .tools import debug_on
from numpy.polynomial.legendre import Legendre
try:
    from pygsplines import BasisLegendre
except ImportError:
    MOD_PATH = pathlib.Path(__file__).parent.absolute()
    MOD_PATH_PYGSPLINES = pathlib.Path(MOD_PATH, '..', 'build')
    sys.path.append(str(MOD_PATH_PYGSPLINES))
    from pygsplines import BasisLegendre


class cMyTest(unittest.TestCase):
    @debug_on()
    def __init__(self, *args, **kwargs):
        super(cMyTest, self).__init__(*args, **kwargs)

    @debug_on()
    def test(self):
        n = np.random.randint(2, 5)
        basis = BasisLegendre(n)

        vec = np.zeros((n,))

        basis_nominal = [
            Legendre([1 if i == j else 0 for j in range(n)]) for i in range(n)]

        for _ in range(100):
            tval = np.random.rand()*2 - 1

            basis.eval_on_window(tval, 2.0, vec)
            for j in range(n):
                assert abs(basis_nominal[j](tval)-vec[j]) < 1.0e-7

        for d in range(1, n):
            print('n={:d} d={:d}'.format(n, d))
            basis.eval_derivative_on_window(tval, 2.0, d, vec)
            for j in range(n):
                assert abs(basis_nominal[j].deriv(d)(tval)-vec[j]) < 1.0e-7


if __name__ == '__main__':
    unittest.main()
