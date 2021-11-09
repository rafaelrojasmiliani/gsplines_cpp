
import matplotlib.pyplot as plt
import time
import pathlib
import sys
import unittest
import numpy as np
from numpy.polynomial import Polynomial
import scipy
from scipy.special import factorial
from scipy.interpolate import lagrange
from scipy.linalg._expm_frechet import _diff_pade3, _diff_pade5, \
    _diff_pade7, _diff_pade9
from scipy.linalg import expm, expm_frechet


try:
    from gsplines.basis import BasisLagrange
    from gsplines.basis import BasisLegendre
    from gsplines.optimization import optimal_sobolev_norm
    from gsplines.manifolds import MatrixExponential

except ImportError:
    MOD_PATH = pathlib.Path(__file__).parent.absolute()
    MOD_PATH_PYGSPLINES = pathlib.Path(MOD_PATH, '..', 'build')
    sys.path.append(str(MOD_PATH_PYGSPLINES))
    from gsplines.basis import BasisLagrange
    from gsplines.basis import BasisLegendre
    from gsplines.optimization import optimal_sobolev_norm
    from gsplines.manifolds import MatrixExponential


class MyTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(MyTest, self).__init__(*args, **kwargs)
        np.set_printoptions(linewidth=500, precision=4)

    def functions(self):
        n = 12
        A = np.random.rand(n, n)
        E = np.random.rand(n, n)
        M1 = np.array(np.random.rand(n, n))
        M2 = np.array(np.random.rand(n, n))
        M3 = np.array(np.random.rand(n, n))
        M4 = np.array(np.random.rand(n, n))

        me = MatrixExponential(n)

        accum = 0

        accum_2 = 0

        for nom_fun, test_fun in zip([_diff_pade3, _diff_pade5, _diff_pade7,
                                      _diff_pade9],
                                     [me.diff_pade3, me.diff_pade5,
                                         me.diff_pade7, me.diff_pade9]):

            t0 = time.time()
            U, V, Lu, Lv = nom_fun(A, E, np.identity(n))
            t1 = time.time()

            t_nom = t1 - t0
            t0 = time.time()
            test_fun(A, E, M1, M2, M3, M4)
            t1 = time.time()
            t_test = t1 - t0
            #print('test', t_test)
            #print('nom ', t_nom)

            accum = (t_test - t_nom)/t_nom
            accum_2 += 1
            for test_val, nom_val, label in zip([U, V, Lu, Lv],
                                                [M1, M2, M3, M4],
                                                ['U', 'V', 'Lv', 'Lu']):
                res = np.max(np.abs(test_val - nom_val))
                if res < 1.0e-6:
                    err = res
                else:
                    err = max(res / np.max(np.abs(test_val)),
                              res / np.max(np.abs(nom_val)))
                self.assertLessEqual(err, 5.0e-8, label +
                                     " " + str(test_fun.__name__))

        #print('time improvement =  {:.3f} \%'.format(accum/accum_2))

    def value(self):
        n = 3
        A = np.random.rand(n, n)
        E = np.random.rand(n, n)
        E /= np.linalg.norm(E)

        A_norm_1 = scipy.linalg.norm(A, 1)
        me = MatrixExponential(n)

        expm_nom, expm_frechet_nom = expm_frechet(A, E)
        print(expm_nom)

        expm_test = me(A)
        expm_frechet_test = me.gataux_derivative(A, E)

        for test_val, nom_val, label in zip([expm_test, expm_frechet_test],
                                            [expm_nom, expm_frechet_nom],
                                            ['expm_nom', 'frechet']):
            res = np.max(np.abs(test_val - nom_val))
            if res < 1.0e-6:
                err = res
            else:
                err = max(res / np.max(np.abs(test_val)),
                          res / np.max(np.abs(nom_val)))
            print(err)
            tol = 5.0e-8
            if(A_norm_1 >= 1.78):
                tol = 5.0e-3
            self.assertLessEqual(err, tol, label)

    def test(self):
        for _ in range(10):
            self.functions()
            self.value()


def main():
    unittest.main()


if __name__ == '__main__':
    main()
