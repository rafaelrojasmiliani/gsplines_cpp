import time
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import Legendre
from .tools import debug_on
import numpy as np
import unittest
import sys
import os
import pathlib
from gsplines.basis import BasisLegendre
from gsplines.basis import BasisLagrange
from gsplines import Interpolator


##


class MyTest(unittest.TestCase):
    # @debug_on()
    def __init__(self, *args, **kwargs):
        super(MyTest, self).__init__(*args, **kwargs)
        np.random.seed()

    # @debug_on()
    def compute_gspline_test(self):
        """ Test the computation of an interpolating gpline """
        basis_dim = 8
        basis = BasisLegendre(basis_dim)
        dim = 6  # np.random.randint(1, 10)
        intervals = np.random.randint(3, 6)
        waypoints = np.random.rand(intervals+1, dim)
        interval_lengths = 1.0 + np.random.rand(intervals) * 2.0
        inter = Interpolator(dim, intervals, basis)
        res = inter.interpolate(interval_lengths, waypoints)

        # plot(res, 5, 0.10, waypoints)

        nodes = np.array([np.cos((2*k-1)/2/basis_dim*np.pi)
                          for k in range(1, basis_dim+1)])

        basis = BasisLagrange(nodes)
        inter = Interpolator(dim, intervals, basis)
        res = inter.interpolate(interval_lengths, waypoints)
        # plot(res, 5, 0.10, waypoints)

#    #@debug_on()
#    def interpolation_test_computation_time(self):
#        """ copmutes the time it takes to interpolate """
#        basis = BasisLegendre(6)
#        dim = 7
#        intervals = 10
#        waypoints = np.random.rand(intervals+1, dim)
#        intervals_length = np.array(intervals*[1])
#
#        iters = 100
#        tmean = 0
#        for _ in range(iters):
#            time_0 = time.time()
#            inter = Interpolator(dim, intervals, basis)
#            _ = inter.interpolate(intervals_length, waypoints)
#            time_1 = time.time()
#            tmean += (time_1 - time_0)
#
#        print('mean time = ', tmean/iters)

    # @debug_on()
    def derivative_y(self):
        ''' Compare the numerical derivate of y w.r.t tau with the nominal one
        '''
        for _ in range(40):
            print('--- begin --')
            b_dim = 2*np.random.randint(1, 10)
            basis = BasisLegendre(b_dim)
            np.random.seed()
            dim = np.random.randint(1, 6)
            intervals = np.random.randint(2, 10)
            waypoints = (np.random.rand(intervals + 1, dim) - 0.5) * 2 * np.pi
            interval_lengths = 1.0 + np.random.rand(intervals) * 2.0
            print('interplator constrcutor codo_dim, N, bdim',
                  dim, intervals, b_dim)
            inter = Interpolator(dim, intervals, basis)
            y_nom = inter.solve_interpolation(
                interval_lengths, waypoints)

            dtau = 1.0e-5
            err = 0.0
            errp = 0.0
            for iinter in range(0, intervals):
                dydtauNom = inter.get_coeff_derivative_wrt_tau(
                    y_nom, interval_lengths, iinter)

                tauv_aux = interval_lengths.copy()
                tauv_aux[iinter] += -2 * dtau
                y0 = inter.solve_interpolation(
                    tauv_aux, waypoints).copy() * (1.0 / 12.0)
                tauv_aux[iinter] += dtau
                y1 = inter.solve_interpolation(
                    tauv_aux, waypoints).copy() * (-2.0 / 3.0)
                tauv_aux[iinter] += 2 * dtau
                y2 = inter.solve_interpolation(
                    tauv_aux, waypoints).copy() * (2.0 / 3.0)
                tauv_aux[iinter] += dtau
                y3 = inter.solve_interpolation(
                    tauv_aux, waypoints).copy() * (-1.0 / 12.0)
                dydtauiTest = (y0 + y1 + y2 + y3) / dtau
                ev = np.abs(dydtauiTest - dydtauNom)
                e = np.max(ev)
                eidx = np.argmax(ev)

                if(abs(dydtauiTest[eidx]) > 1.0e-5):
                    ep = e / abs(dydtauiTest[eidx])
                else:
                    ep = e

                if e > err:
                    err = e
                if ep > errp:
                    errp = ep

                self.assertTrue(ep < 5.0e-2, '''
                error on dydtau = {:10.7e}
                relative error  = {:10.7e}
                value of dydtau test = {:10.7e}
                value of dydtau nom = {:10.7e}
                '''.format(e, ep, dydtauiTest[eidx], dydtauNom[eidx]))

    def test(self):

        self.compute_gspline_test()


if __name__ == '__main__':
    unittest.main()
