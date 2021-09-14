
import pathlib
import sys
import unittest
import numpy as np
from numpy.polynomial.legendre import Legendre
#from .tools import debug_on
try:
    from gsplines.basis import BasisLegendre
except ImportError:
    MOD_PATH = pathlib.Path(__file__).parent.absolute()
    MOD_PATH_PYGSPLINES = pathlib.Path(MOD_PATH, '..', 'build')
    sys.path.append(str(MOD_PATH_PYGSPLINES))
    from gsplines.basis import BasisLegendre


class MyTest(unittest.TestCase):
    # @debug_on()
    def __init__(self, *args, **kwargs):
        super(MyTest, self).__init__(*args, **kwargs)

    # @debug_on()
    def test(self):
        n = 2*np.random.randint(1, 8)
        basis = BasisLegendre(n)

        vec = np.zeros((n,))

        basis_nominal = [
            Legendre([1 if i == j else 0 for j in range(n)]) for i in range(n)]

        for _ in range(100):
            tval = np.random.rand()*2 - 1

            basis.eval_on_window(tval, 2.0, vec)
            for j in range(n):
                self.assertLess(abs(basis_nominal[j](
                    tval)-vec[j]), 1.0e-7, "Error on value")

        for d in range(1, n):
            basis.eval_derivative_on_window(tval, 2.0, d, vec)
            for j in range(n):
                nom_val = basis_nominal[j].deriv(d)(
                    tval)
                test_val = vec[j]

                if abs(nom_val) < 1.0e-5:
                    self.assertTrue(abs(nom_val - test_val) <
                                    1.0e-7, 'error on derivative')
                else:
                    self.assertTrue(abs(nom_val - test_val) /
                                    abs(nom_val) < 1.0e-7,
                                    'error on derivative')


if __name__ == '__main__':
    unittest.main()
