"""
    Test functions in the space of solutions of the
    Euler Lagrange equations of

    \int_{-1}^{1} (2/tau) \alpha dq/ds + (2/tau)^5 (1-\alpha) d^3 q / ds^3 dt
"""
import unittest
from scipy.sparse import csc_matrix
import numpy as np
from gsplines.interpolator.gspline import cSplineCalc
from gsplines.basis.basis1010 import cBasis1010
from gsplines.basis.basis0010 import cBasis0010


class cMyTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(cMyTest, self).__init__(*args, **kwargs)
        import sys
        np.set_printoptions(
            linewidth=5000000,
            formatter={'float': '{:+10.3e}'.format},
            threshold=sys.maxsize)
        pass

    def testInversion(self):
        import time
        print('Testinf inversion of matrix')
        dim = 6  # np.random.randint(2, 6)
        N = 50  # np.random.randint(3, 120)
        a = np.random.rand()
        splcalc = cSplineCalc(dim, N, cBasis1010(a))
        for i in range(50):
            tauv = np.random.rand(N)
            A1 = splcalc.eval_A(tauv)
#            A0 = self.eval_A(tauv, dim, N, cBasis1010(a))
#            e = np.max(np.abs(A1 - A0.todense()))
#            #            print(A0)
#            #            print('----------------------------------------')
#            #            print(A1)
#            #            print(dim, N)
#            assert e < 1.0e-8

        splcalc.printPerformace()
        pass

    def testcontinuity(self):
        print('Test continuity constraints with plot')
        for i in range(3):
            dim = np.random.randint(2, 3)
            N = np.random.randint(3, 10)
            a = np.random.rand()
            wp = (np.random.rand(N + 1, dim) - 0.5) * 2 * np.pi
            tauv = 0.5 + np.random.rand(N) * 3.0
            tis = [np.sum(tauv[0:i]) for i in range(0, N + 1)]
            T = np.sum(tauv)
            splcalc = cSplineCalc(dim, N, cBasis1010(a))
            spln = splcalc.getSpline(tauv, wp)
            from matplotlib import pyplot as plt

            t = np.arange(0, T, 0.005)
            q_list = [spln.deriv(i)(t) for i in range(0, 6)]

            fig, axs = plt.subplots(6, dim)

            for i in range(0, 6):
                for j in range(0, dim):
                    axs[i, j].plot(t, q_list[i][:, j])
                    axs[i, j].grid()
                    for ti in tis:
                        axs[i, j].axvline(x=ti, color='b', linestyle='--')

            plt.show()

    def test_eval_b(self):
        import time
        print('Test evaluation of b vector')
        for i in range(20):
            dim = np.random.randint(1, 8)
            N = np.random.randint(3, 200)
            a = np.random.rand()
            wp = (np.random.rand(N + 1, dim) - 0.5) * 2 * np.pi
            dwp0 = np.zeros((dim, ))
            ddwp0 = np.zeros((dim, ))
            dwpT = np.zeros((dim, ))
            ddwpT = np.zeros((dim, ))
            splcalc = cSplineCalc(dim, N, cBasis1010(a))
            b1 = splcalc.eval_b(wp)
            b2 = self.eval_b(wp, N, dim, dwp0, ddwp0, dwpT, ddwpT)

            e = np.max(np.abs(b1 - b2))

            assert e < 1.0e-8

    def non_zero_diagonals_A(self):
        for i in range(0, 1):
            dim = np.random.randint(1, 8)
            N = np.random.randint(3, 200)
            a = np.random.rand()
            splcalc = cSplineCalc(dim, N, cBasis1010, a)

            tauv = 10.0 * np.random.rand(N) + 0.5
            splcalc.eval_A(tauv)
            A = splcalc.Aeval

            upper_diags = 0
            flag = 0
            for i in range(A.shape[0]):
                if np.max(np.abs(np.diagonal(A, i))) > 1.0e-10:
                    assert flag != 1, 'Matrix is Not Banded!!!'
                    upper_diags += 1
                else:
                    flag = 1

            lower_diags = 0
            flag = 0
            for i in range(A.shape[0]):
                if np.max(np.abs(np.diagonal(A, -i))) > 1.0e-10:
                    assert flag != 1, 'Matrix is Not Banded!!!'
                    lower_diags += 1
                else:
                    flag = 1

#            print('upper diagonas = {:d}  lower diagonals = {:d}'.format(
#                upper_diags, lower_diags))

            assert 4 * dim + 4 == max(upper_diags, lower_diags)

#        wp = (np.random.rand(N + 1, dim) - 0.5) * 2.0 * np.pi
#        b = splcalc.eval_b(wp)
#
#        wp = np.random.rand(N + 1, dim)
#
#        spline = splcalc.solve(wp, tauv)
#
#        tis = [np.sum(tauv[0:i]) for i in range(0, N + 1)]
#
#        t = np.arange(0, tis[-1], 0.1)
#
#        plt.plot(t, spline(t)[:, 0])
#
#        plt.show()

    def test_derivative_b(self):
        '''
            Here we rest the correctness of the numerical output of the basis
            class comparing it with its analitical form optained using sympy
        '''
        print('Test derivative of b w.r.t. waypoint components')
        np.random.seed()
        dim = np.random.randint(1, 8)
        N = np.random.randint(2, 60)
        a = np.random.rand()
        wp = (np.random.rand(N + 1, dim) - 0.5) * 2 * np.pi
        splcalc = cSplineCalc(dim, N, cBasis1010(a))

        dwp = 0.0005
        for i in range(N + 1):
            for j in range(dim):
                wpidx = i
                i = j

                wp_aux = wp.copy()
                wp_aux[wpidx, i] += -dwp
                b1 = splcalc.eval_b(wp_aux).copy()
                wp_aux[wpidx, i] += 2 * dwp
                b2 = splcalc.eval_b(wp_aux).copy()

                dbdwpij_num = 0.5 * (b2 - b1) / dwp

                dbdwpij = splcalc.eval_dbdwpij(wpidx, i)

                e = np.max(np.abs(dbdwpij_num - dbdwpij))

                if e > 1.0e-8:
                    print('Erroe in db_dwpij:')
                    print('implementation:')
                    print(dbdwpij)
                    print('(b1-b2)/dwp:')
                    print(dbdwpij_num)
                    print('component', i)
                    print('waypoint ', wpidx)
                    print('dimension ', dim)
                    print('number of intervals ', N)
                    raise AssertionError('Error of {:14.7e}'.format(e))

    def test_derivative_y(self):
        ''' Compare the numerical derivate of y w.r.t tau with the nominal one
        '''
        for i in range(40):
            np.random.seed()
            dim = np.random.randint(1, 3)
            N = np.random.randint(2, 6)
            a = np.random.rand()
            wp = (np.random.rand(N + 1, dim) - 0.5) * 2 * np.pi
            tauv = 1.0 + np.random.rand(N) * 2.0
            splcalc = cSplineCalc(dim, N, cBasis1010(a))
            dydtauNom, y = splcalc.eval_dydtau(tauv, wp)

            y = y.copy()

            dtau = 1.0e-8
            err = 0.0
            errp = 0.0
            for iinter in range(0, N):
                tauv_aux = tauv.copy()
                tauv_aux[iinter] += -2 * dtau
                y0 = splcalc.eval_y(tauv_aux, wp).copy() * (1.0 / 12.0)
                tauv_aux[iinter] += dtau
                y1 = splcalc.eval_y(tauv_aux, wp).copy() * (-2.0 / 3.0)
                tauv_aux[iinter] += 2 * dtau
                y2 = splcalc.eval_y(tauv_aux, wp).copy() * (2.0 / 3.0)
                tauv_aux[iinter] += dtau
                y3 = splcalc.eval_y(tauv_aux, wp).copy() * (-1.0 / 12.0)
                dydtauiTest = (y0 + y1 + y2 + y3) / dtau
                ev = np.abs(dydtauiTest - dydtauNom[:, iinter])
                e = np.max(ev)
                eidx = np.argmax(ev)

                ep = e / abs(dydtauiTest[eidx])

                if e > err:
                    err = e
                if ep > errp:
                    errp = ep

                assert ep < 5.0e-2, '''
                error on dydtau = {:10.7e}
                value of dydtau = {:10.7e}
                relative error  = {:10.7e}
                '''.format(e, dydtauiTest[eidx], ep)

    def test_derivative_wp(self):
        ''' Compare the numerical derivate of y w.r.t waypoints with the nominal one
        '''
        for _ in range(4):
            np.random.seed()
            dim = np.random.randint(1, 8)
            N = np.random.randint(2, 20)
            a = np.random.rand()
            wp = (np.random.rand(N + 1, dim) - 0.5) * 2 * np.pi
            tauv = 0.5 + np.random.rand(N) * 2.0
            splcalc = cSplineCalc(dim, N, cBasis1010(a))
            _, y = splcalc.eval_dydtau(tauv, wp)
            y = y.copy()

            err = 0.0
            errp = 0.0

            err = 0.0
            errp = 0.0
            dwp = 1.0e-5

            wpidx = [(i, j) for i in range(N + 1) for j in range(dim)]
            dydwpNom = np.zeros((y.shape[0], len(wpidx)))
            dydwpNom, _ = splcalc.eval_dydu(tauv, wp, wpidx, dydwpNom)

            for k, (i, j) in enumerate(wpidx):
                wp_aux = wp.copy()
                wpidx = i
                wpcom = j

                wp_aux[wpidx, wpcom] += -3 * dwp
                y0 = splcalc.eval_y(tauv, wp_aux).copy() * (-1.0 / 60.0)
                wp_aux[wpidx, wpcom] += dwp
                y1 = splcalc.eval_y(tauv, wp_aux).copy() * (3.0 / 20.0)
                wp_aux[wpidx, wpcom] += dwp
                y2 = splcalc.eval_y(tauv, wp_aux).copy() * (-3.0 / 4.0)
                wp_aux[wpidx, wpcom] += 2 * dwp
                y3 = splcalc.eval_y(tauv, wp_aux).copy() * (3.0 / 4.0)
                wp_aux[wpidx, wpcom] += dwp
                y4 = splcalc.eval_y(tauv, wp_aux).copy() * (-3.0 / 20.0)
                wp_aux[wpidx, wpcom] += dwp
                y5 = splcalc.eval_y(tauv, wp_aux).copy() * (1.0 / 60.0)

                dydwpTest = (y0 + y1 + y2 + y3 + y4 + y5) / dwp

                ev = np.abs(dydwpNom[:, k] - dydwpTest)
                e = np.max(ev)
                eidx = np.argmax(ev)
#                print('{:14.7e} {:14.7e} {:14.7e}'.format(
#                    e, dydwpNom[eidx, k], dydwpTest[eidx]))

                ep = e / dydwpTest[eidx]

                if e > err:
                    err = e
                if ep > errp:
                    errp = ep

                if e > 1.0e-4:
                    assert ep < 1.0e-8, '''
                    Relative Error   = {:10.3e}
                    Absolute Error   = {:10.3e}
                    '''.format(ep, e)


#            print('Maximum Error for dy dwp           = {:14.7e}'.format(err))
#            print('Maximum Relative Error for dy dwp  = {:14.7e}'.format(errp))

#                assert e < 5.0e-2, 'error = {:14.7e}'.format(e)

    def test_derivative_y_2(self):
        ''' Second test for the derivative of y wr.t. tau.
        First test the identity A*dydtau + dAdtau y = 0A
        Second test the identity of above using basis
        '''
        for i in range(40):
            np.random.seed()
            dim = np.random.randint(1, 3)
            N = np.random.randint(2, 6)
            a = np.random.rand()
            wp = (np.random.rand(N + 1, dim) - 0.5) * 2 * np.pi
            tauv = 0.5 + np.random.rand(N) * 2.0
            splcalc = cSplineCalc(dim, N, cBasis0010())
            basis = splcalc.basis_
            dydtauNom, y = splcalc.eval_dydtau(tauv, wp)
            A = splcalc.eval_A(tauv).todense()
            y = y.copy()

            # A*dydtau + dAdtau y = 0
            for i, _ in enumerate(tauv):
                v0 = A.dot(dydtauNom[:, i])
                v1 = splcalc.eval_dAdtiy(tauv, i, y).todense()
                res = v0 + v1.T
                e = np.abs(res)

                assert np.max(e) < 1.0e-10, '''
                e = {:14.7e}
                '''.format(e)

            for i, taui in enumerate(tauv):

                B = basis.evalDerivOnWindow(-1, taui, 0)
                dB_dtau = basis.evalDerivWrtTauOnWindow(-1, taui, 0)

                for idim in range(dim):
                    i0 = i * 6 * dim + 6 * idim
                    i1 = i0 + 6
                    e = dydtauNom[i0:i1, i].dot(B) + y[i0:i1].dot(dB_dtau)

                    assert np.abs(e) < 1.0e-10, '''
                    error computing dydtau^\\top B + y^\\top dB_dtau
                    e            = {:14.7e}
                    index of tau = {:d}
                    index of q   = {:d}
                    i0           = {:d}
                    i1           = {:d}
                    '''.format(e, i, idim, i0, i1)

    def eval_A(self, tau, _dim, _N, _basis):
        """
        Alternative way to fill the Matrix A
        WARNING:  This work ok for _N<120
        """

        A = np.zeros(2 * (6 * _dim * _N, ))

        dim = _dim
        nzv = np.zeros((((_N - 1) * dim + (3) * dim) * 4 * 3 + (3 - 1) *
                        (_N - 1) * dim * 8 * 3, ))
        idxv = np.zeros(
            (((_N - 1) * dim + (3) * dim) * 4 * 3 + (3 - 1) *
             (_N - 1) * dim * 8 * 3, ),
            dtype=np.int16)
        ptrv = np.zeros((2 * 3 * _N * dim + 1, ), dtype=np.int16)

        basis = _basis

        Pl = [basis.evalDerivOnWindow(-1.0, tau[0], i) for i in range(0, 5)]
        Pr = [basis.evalDerivOnWindow(1.0, tau[0], i) for i in range(0, 5)]
        # Fill the content for the derivatives at boundaries
        Cpl = -np.vstack(Pl[1:5])
        Cpr = np.vstack(Pr[1:5])

        App1 = np.vstack([Pl[:3], Pr[0]])

        nnz = 0
        nptr = 0
        # --------------------------------------------
        for j in range(0, 6 * dim):
            ptrv[nptr] = nnz
            i0 = (j // 6) * 4
            for i in range(i0, i0 + 4):
                A[i, j] = App1[(i - i0), j % 6]
                nzv[nnz] = App1[(i - i0), j % 6]
                idxv[nnz] = i
                nnz += 1
            i0 += 4 * dim
            for i in range(i0, i0 + 4):
                A[i, j] = Cpr[(i - i0), j % 6]
                nzv[nnz] = Cpr[(i - i0), j % 6]
                idxv[nnz] = i
                nnz += 1
            nptr += 1
        # --------------------------------------------
        for iinter in range(1, _N - 1):
            j0 = iinter * 6 * dim
            i0 = 4 * dim + 6 * dim * (iinter - 1)

            Pl = [
                basis.evalDerivOnWindow(-1.0, tau[iinter], i)
                for i in range(0, 5)
            ]
            Pr = [
                basis.evalDerivOnWindow(1.0, tau[iinter], i)
                for i in range(0, 5)
            ]
            Cpl = -np.vstack(Pl[1:5])
            Cpr = np.vstack(Pr[1:5])
            Appi = np.vstack([Pl[0], Pr[0]])
            # --------------------------------------------
            for j in range(j0, j0 + 6 * dim):
                ptrv[nptr] = nnz
                i1 = i0 + ((j - j0) // 6) * 4
                for i in range(i1, i1 + 4):
                    A[i, j] = Cpl[(i - i1), j % 6]
                    nzv[nnz] = Cpl[(i - i1), j % 6]
                    idxv[nnz] = i
                    nnz += 1
                i1 += (dim - (j - j0) // 6) * 4 + ((j - j0) // 6) * 2
                for i in range(i1, i1 + 2):
                    A[i, j] = Appi[i - i1, j % 6]
                    nzv[nnz] = Appi[i - i1, j % 6]
                    idxv[nnz] = i
                    nnz += 1
                i1 += (dim - (j - j0) // 6) * 2 + ((j - j0) // 6) * 4
                for i in range(i1, i1 + 4):
                    A[i, j] = Cpr[(i - i1), j % 6]
                    nzv[nnz] = Cpr[(i - i1), j % 6]
                    idxv[nnz] = i
                    nnz += 1
                nptr += 1
        # -------------------------------------
        #     Last column of the matrix A
        # ------------------------------------
        i0 = 4 * dim + 6 * dim * (_N - 2)
        j0 = (_N - 1) * 6 * dim
        Pl = [basis.evalDerivOnWindow(-1.0, tau[-1], i) for i in range(0, 5)]
        Pr = [basis.evalDerivOnWindow(1.0, tau[-1], i) for i in range(0, 5)]
        Cpl = -np.vstack(Pl[1:5])
        AppN = np.vstack([Pl[0], Pr[1:3], Pr[0]])
        # --------------------------------------------
        for j in range(j0, j0 + 6 * dim):
            ptrv[nptr] = nnz
            i1 = i0 + ((j - j0) // 6) * 4
            for i in range(i1, i1 + 4):
                A[i, j] = Cpl[(i - i1), j % 6]
                nzv[nnz] = Cpl[(i - i1), j % 6]
                idxv[nnz] = i
                nnz += 1
            i1 += (dim - (j - j0) // 6) * 4 + ((j - j0) // 6) * 4
            for i in range(i1, i1 + 4):
                A[i, j] = AppN[i - i1, j % 6]
                nzv[nnz] = AppN[(i - i1), j % 6]
                idxv[nnz] = i
                nnz += 1
            nptr += 1

        ptrv[nptr] = nnz

        return csc_matrix((nzv, idxv, ptrv), shape=2 * (2 * 3 * _N * dim, ))

    def eval_b(self, _wp, _N, _dim, _dwp0, _dwpT, _ddw0, _ddwT):
        assert _wp.shape[0] == _N + 1 and _wp.shape[1] == _dim
        dim = _dim
        b = np.zeros((2 * _N * 3 * _dim, ))

        for i in range(dim):
            b[4 * i:4 * i + 4] = (_wp[0][i], _dwp0[i], _ddw0[i], _wp[1][i])

        i0 = 8 * dim

        for idxwp, _ in enumerate(_wp[1:-2]):
            idxwp += 1
            for i in range(dim):
                b[i0 + 2 * i:i0 + 2 * i + 2] = (_wp[idxwp][i],
                                                _wp[idxwp + 1][i])
            i0 += 6 * dim

        for i in range(dim):
            b[i0 + 4 * i:i0 + 4 * i + 4] = (_wp[-2][i], _dwpT[i], _ddwT[i],
                                            _wp[-1][i])
        return b

    def test_derivative_A(self):
        for i in range(5):
            dim = 6  # np.random.randint(2, 6)
            N = 20  # np.random.randint(3, 120)
            a = np.random.rand()
            splcalc = cSplineCalc(dim, N, cBasis1010(a))

            tauv0 = 0.5 + np.random.rand(N)
            dtau = 0.0001
            z = np.ones((splcalc.linsys_shape_, ))
            for iinter, taui in enumerate(tauv0):
                tauv1 = tauv0.copy()
                tauv1[iinter] += -dtau
                A0 = splcalc.eval_A(tauv1).copy()
                tauv1[iinter] += 2 * dtau
                A1 = splcalc.eval_A(tauv1)

                dAdtaui = 0.5 * (A1 - A0) / dtau

                a = np.max(np.abs(dAdtaui))

                testVal = dAdtaui.dot(z)
                nomVal = np.ravel(
                    splcalc.eval_dAdtiy(tauv0, iinter, z).todense())

                ev = np.abs(testVal - nomVal)
                testValMax = np.max(testVal)
                nomValMax = np.max(nomVal)
                e = np.max(ev)
                ep = e / (nomValMax)
                assert ep < 5.0e-6, '''
                error              = {:10.3e}
                max nominal value  = {:10.3e}
                max test value     = {:10.3e}
                tau compoenet      = {:d}
                interation         = {:d}
                '''.format(e, nomValMax, testValMax, iinter, i)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
