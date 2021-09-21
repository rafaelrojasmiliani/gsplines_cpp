
import pathlib
import sys
import unittest
import numpy as np
from numpy.polynomial import Polynomial
from scipy.special import factorial
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt

try:
    from gsplines.basis import BasisLagrange
    from gsplines.basis import BasisLegendre
    from gsplines.optimization import optimal_sobolev_norm

except ImportError:
    MOD_PATH = pathlib.Path(__file__).parent.absolute()
    MOD_PATH_PYGSPLINES = pathlib.Path(MOD_PATH, '..', 'build')
    sys.path.append(str(MOD_PATH_PYGSPLINES))
    from gsplines.basis import BasisLagrange
    from gsplines.basis import BasisLegendre
    from gsplines.optimization import optimal_sobolev_norm


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

    def plot_optimization(self):

        n_nodes = 6
        nodes = np.array([np.cos((2*k-1)/2/n_nodes*np.pi)
                          for k in range(1, n_nodes+1)])

        dim = 7  # np.random.randint(1, 10)
        intervals = 4

        exec_time = intervals
        waypoints = np.random.rand(intervals+1, dim)*6.14
        basis_lagrange = BasisLagrange(nodes)
        basis_legendre = BasisLegendre(n_nodes)
        curve_1 = optimal_sobolev_norm(
            waypoints, basis_lagrange, [(3, 1)], exec_time)
        curve_2 = optimal_sobolev_norm(
            waypoints, basis_legendre, [(3, 1)], exec_time)
        import pdb
        pdb.set_trace()
        # curve_1 = op
        show_piecewisefunction(curve_1, 5, 0.001)
        show_piecewisefunction(curve_2, 5, 0.001)
        plt.show()
        pass

    def test(self):
        for _ in range(10):
            self.value()

        self.plot_optimization()


def main():
    unittest.main()


if __name__ == '__main__':
    main()
