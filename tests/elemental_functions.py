"""
Test use of FunctionExpression bindings
"""
import pathlib
import sys
import unittest
import numpy as np
import matplotlib.pyplot as plt
import time
import functools

from numpy.polynomial.polynomial import Polynomial
import sympy as sp
from .tools import debug_on

np.random.seed()

try:
    from pygsplines import Sin, Cos, Exponential, Identity, ConstFunction, CanonicPolynomial
except ImportError:
    MOD_PATH = pathlib.Path(__file__).parent.absolute()
    MOD_PATH_PYGSPLINES = pathlib.Path(MOD_PATH, '..', 'build')
    sys.path.append(str(MOD_PATH_PYGSPLINES))
    from pygsplines import Sin, Cos, Exponential, Identity, ConstFunction, CanonicPolynomial


def sp_identity(_var):
    """ identity wrapper """
    return _var


class SpConst:
    """ Wrapper to a constant function """

    def __init__(self, _value):
        self.value_ = _value

    def __call__(self, _t):
        """ evaluates """
        return self.value_

    def deriv(self, _deg=1):
        """ returns derivative"""
        return 0.0


def get_random_elemental_fuction():
    """ returns a random funcion """
    option = np.random.randint(0, 5)
    if option == 0:
        result = Sin, sp.sin
    elif option == 1:
        result = Cos, sp.cos
    elif option == 2:
        result = Exponential, sp.exp
    elif option == 3:
        result = Identity, sp_identity
    elif option == 4:
        value = np.random.rand()

        def aux(_dom):
            return ConstFunction(_dom, 1, value)

        result = aux, SpConst(value)

    return result


def get_random_function(_dom, sp_var):
    """ REturns a random function expression """
    num = np.random.randint(5, 10)

    result_gs = Identity((-1, 1))

    result_sp = sp_var

    for _ in range(num):
        operation = np.random.randint(0, 4)
        fun_gs, fun_sp = get_random_elemental_fuction()

        if operation == 0:
            result_sp = result_sp + fun_sp(sp_var)
            result_gs = result_gs + fun_gs(_dom)
        elif operation == 1:
            result_sp = result_sp * fun_sp(sp_var)
            result_gs = result_gs * fun_gs(_dom)
        elif operation == 2:
            result_sp = result_sp - fun_sp(sp_var)
            result_gs = result_gs - fun_gs(_dom)
        elif operation == 3:
            result_sp = result_sp.subs({sp_var: fun_sp(sp_var)})
            result_gs = result_gs.compose(fun_gs(_dom))

    return result_sp, result_gs


def show_piecewisefunction(_q, _up_to_deriv=3, _dt=0.1, _title=''):
    dim = _q.get_codom_dim()
    fig, axes = plt.subplots(_up_to_deriv + 1, dim)
    if dim == 1:
        axes = np.array([[axes[i]] for i in range(_up_to_deriv + 1)])
    if _title:
        fig.suptitle(_title)
    t = np.arange(0.0, _q.get_exec_time(), _dt)

    for i in range(0, _up_to_deriv + 1):
        q = _q.deriv(i)
        qt = q(t)
        for j in range(0, dim):
            axes[i, j].plot(t, qt[:, j])
            axes[i, j].grid()
            if i == 0:
                axes[i, j].set_title(
                    'coordinate {:d}'.format(j + 1), fontsize=8)

            if hasattr(_q, 'get_domain_breakpoints'):
                for ti in _q.get_domain_breakpoints():
                    axes[i, j].axvline(ti, alpha=0.1, color='red')

    plt.subplots_adjust(
        left=0.025,
        bottom=0.05,
        right=0.975,
        top=0.95,
        wspace=0.25,
        hspace=0.15)
    plt.show()


class MyTest(unittest.TestCase):
    """ Test elemental function and operations"""
    @debug_on()
    def __init__(self, *args, **kwargs):
        unittest.TestCase.__init__(self, *args, **kwargs)
        np.random.seed()
        self.exp = Exponential((-1, 1))
        self.sin = Sin((-1, 1))
        self.cos = Cos((-1, 1))
        self.identity = Identity((-1, 1))
        self.const_function = functools.partial(ConstFunction, (-1, 1))
        self.time_spam = np.reshape(np.arange(-1, 1, 0.2), (-1, 1))

    @debug_on()
    def error_test(self, _f_nom, _f_test):
        """ calls np.linalg.norm """
        nom_val = _f_nom(self.time_spam)
        nom_val_max = np.abs(np.max(nom_val))
        if nom_val_max > 1.0e-9:
            assert(np.linalg.norm(
                _f_test(self.time_spam) -
                nom_val)/nom_val_max < 1.0e-9)
        else:
            assert(np.linalg.norm(
                _f_test(self.time_spam) -
                nom_val) < 1.0e-9)

    @debug_on()
    def add_test(self):
        """ test sum """

        f_test = self.exp + self.sin

        def f_nom(_t):
            return np.exp(_t) + np.sin(_t)

        self.error_test(f_nom, f_test)

    @debug_on()
    def mul_test(self):
        """ TEst multiplication """

        f_test = self.exp * self.sin

        def f_nom(_t):
            return np.multiply(np.exp(_t), np.sin(_t))

        self.error_test(f_nom, f_test)

    @debug_on()
    def subs_test(self):
        """ TEst multiplication """

        f_test = self.exp - self.sin

        def f_nom(_t):
            return np.exp(_t) - np.sin(_t)

        self.error_test(f_nom, f_test)

    @debug_on()
    def comp_test(self):
        """ TEst multiplication """

        f_test = self.exp.compose(self.sin)

        def f_nom(_t):
            return np.exp(np.sin(_t))

        self.error_test(f_nom, f_test)

    @debug_on()
    def deriv_test(self):
        """ TEst multiplication """

        sp_var = sp.symbols('t', real=True)

        for _ in range(100):
            f_nom, f_test = get_random_function((-1, 1), sp_var)

            f_nom_eval = sp.lambdify(sp_var, f_nom)

            self.error_test(f_nom_eval, f_test)

            f_d_nom = f_nom.diff()
            f_d_test = f_test.deriv()
            f_d_nom_eval = sp.lambdify(sp_var, f_d_nom)

            self.error_test(f_d_nom_eval, f_d_test)

    @debug_on()
    def polynomial_test(self):
        for _ in range(100):
            coeff = np.random.rand(2+np.random.randint(1, 5))
            pol_nom = Polynomial(coeff)
            pol_test = CanonicPolynomial((-1, 1), coeff)

            self.error_test(pol_nom, pol_test)

            for _ in range(6):
                deg = np.random.randint(0, 10)
                pol_d_nom = pol_nom.deriv(deg)
                pol_d_test = pol_test.deriv(deg)
                self.error_test(pol_d_nom, pol_d_test)

    @ debug_on()
    def test(self):
        """ runs all tests"""
        self.add_test()
        self.mul_test()
        self.subs_test()
        self.comp_test()
        self.deriv_test()
        self.polynomial_test()


if __name__ == '__main__':
    unittest.main()
