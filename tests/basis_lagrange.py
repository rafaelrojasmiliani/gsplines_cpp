
import pathlib
import sys
import unittest
import numpy as np
from numpy.polynomial import Polynomial
from scipy.special import factorial
from scipy.interpolate import lagrange

try:
    from gsplines.basis import BasisLagrange
except ImportError:
    MOD_PATH = pathlib.Path(__file__).parent.absolute()
    MOD_PATH_PYGSPLINES = pathlib.Path(MOD_PATH, '..', 'build')
    sys.path.append(str(MOD_PATH_PYGSPLINES))
    from gsplines.basis import BasisLagrange


def get_vector_with_one_one(_dim: int, _one_index: int) -> np.array:
    return np.array([1.0 if i == _one_index else 0.0 for i in range(_dim)])


class MyTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(MyTest, self).__init__(*args, **kwargs)
        np.set_printoptions(linewidth=500, precision=4)

    def value(self):
        """ Test the approximation quality of the approximation of the cosine
        function on [-1, 1] by the lagrange intepolators using the Chebychev
        points"""

        n_nodes = np.random.randint(3, 10)
        nodes = np.array([np.cos((2*k-1)/2/n_nodes*np.pi)
                          for k in range(1, n_nodes+1)])

        basis = BasisLagrange(nodes)

        nom_lagrange_pol = [Polynomial(lagrange(
            nodes, get_vector_with_one_one(n_nodes, i)).coef[::-1])
            for i in range(n_nodes)]

        def nom_basis(_val_in_window: float):
            return np.array([pol(_val_in_window) for pol in nom_lagrange_pol])

        def nom_basis_deriv(_val_in_window: float, _deg: int):
            return np.array([pol.deriv(_deg)(_val_in_window)
                             for pol in nom_lagrange_pol])

        val_in_window = 2*np.random.rand()-1
        buff = np.zeros((n_nodes, ))
        basis.eval_on_window(val_in_window, 2, buff)

        err = np.abs(buff-nom_basis(val_in_window))

        self.assertLessEqual(np.max(err), 1.0e-9)

        for i in range(1, 9):
            basis.eval_derivative_on_window(val_in_window, 2, i, buff)

            test_inf = np.max(np.abs(buff))
            nom_inf = np.max(np.abs(nom_basis_deriv(val_in_window, i)))
            delta_inf = np.max(np.abs(buff-nom_basis_deriv(val_in_window, i)))
            epsilon = 1.0e-9
            if test_inf < epsilon or nom_inf < epsilon:
                err = delta_inf
            else:
                err = max(delta_inf/test_inf, delta_inf/nom_inf)

            self.assertLessEqual(err, 5.0e-6, "deg = {:d}".format(i))

    def test(self):
        for _ in range(10):
            self.value()


def main():
    unittest.main()


if __name__ == '__main__':
    main()
