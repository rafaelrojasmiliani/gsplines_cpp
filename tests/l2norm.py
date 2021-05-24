"""
    Test the cost function from the problem 1010
"""
import numpy as np
import sympy as sp
from scipy.special import roots_legendre
import unittest
import sys
import pathlib
from gsplines.functionals.l2norm import cJerkL2Norm
from gsplines.basis.basis0010 import cBasis0010

from itertools import tee


def pairwise(iterable):
    '''s -> (s0,s1), (s1,s2), (s2, s3), ...'''
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


class cMyTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(cMyTest, self).__init__(*args, **kwargs)
        self.N_ = N = 20  # np.random.randint(2, 10)
        self.dim_ = dim = 8  # np.random.randint(1, 8)
        self.wp_ = wp = np.random.rand(N + 1, dim)
        self.basis_ = basis = cBasis0010()
        self.cost_ = cJerkL2Norm(wp, basis)

    def test_run(self):
        cost = self.cost_
        N = self.N_

        res = np.zeros((N, ))
        for i in range(5):
            tauv = np.random.rand(N)
            a = cost(tauv)
            g = cost.gradient(tauv, res)


    def test_value(self):
        wp = self.wp_
        N = self.N_
        tauv = np.random.rand(N)
        T = np.sum(tauv)

        cost = self.cost_
        q = cost.splcalc_.getSpline(tauv, wp)

        qddd = q.deriv(3)

        def qd3norm(t):
            return np.einsum('ij,ij->i', qddd(t), qddd(t))

        def runningcost(t):
            return qd3norm(t)

        Inom = cost(tauv)
        err = 1.e100
        badtrentCounter = 0
        for Ngl in range(100, 500, 100):
            lr, lw = roots_legendre(Ngl)
            time_partition = np.linspace(0, T, 50)
            Itest_1 = 0.0
            for t0, tf in pairwise(time_partition):
                Itest_1 += sum([w*runningcost(s)[0]*(tf-t0)/2.0 for w, s in zip(lw, (lr+1.0)/2.0*(tf-t0)+t0)])
            Itest_2 = q.l2_norm(3)
            if abs(Itest_1 - Inom) > err:
                badtrentCounter += 1
            else:
                badtrentCounter = 0
            assert badtrentCounter < 3
            err = abs(Itest_1 - Inom)
            err2 = abs(Itest_2 - Inom)
            print('Error w.r.t quadrature = {:.3f}'.format(err))
            print('Error w.r.t gspline impl, = {:.3f}'.format(err2))
            print('Nomial cost  = {:14.7e}, Quadrature = {:14.7e}, spl. imp = {:14.7e}'.format(Inom, Itest_1, Itest_2))
            if err < 1.0e-4:
                break


    def test_gradient(self):
        wp = self.wp_
        N = self.N_

        cost = self.cost_
        nom_grad = np.zeros((N, ))
        test_grad = np.zeros((N, ))
        dt = 1.0e-10
        for _ in range(3):
            tauv = np.random.rand(N)
            cost.gradient(tauv, nom_grad)
            for i in range(N):
                tauv_aux = tauv.copy()
                tauv_aux[i] -= dt
                cost0 = cost(tauv_aux)
                tauv_aux[i] += 2*dt
                cost1 = cost(tauv_aux)
                test_grad[i] = (cost1 - cost0)/(2*dt)
            print(test_grad)
            print(nom_grad)

            print(np.linalg.norm(test_grad - nom_grad))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
