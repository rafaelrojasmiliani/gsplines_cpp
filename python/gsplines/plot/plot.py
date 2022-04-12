import gsplines

import matplotlib.pyplot as plt
import numpy as np


def plot(_q, _up_to_deriv=3, _dt=0.1, _wp=None, _title='', _show=True):

    if hasattr(_q, 'get_waypoints'):
        _wp = _q.get_waypoints()

    dim = _q.get_codom_dim()
    fig, axis = plt.subplots(_up_to_deriv + 1, dim)
    if dim == 1:
        axis = np.array([[axis[i]] for i in range(_up_to_deriv + 1)])
    if _title:
        fig.suptitle(_title)
    time_spam = np.arange(0.0, _q.get_domain_length(), _dt)

    for i in range(0, _up_to_deriv + 1):
        current_curve = _q.deriv(i)
        current_curve_evaluated = current_curve(time_spam)
        for j in range(0, dim):
            axis[i, j].plot(time_spam, current_curve_evaluated[:, j])
            axis[i, j].grid()
            if i == 0:
                axis[i, j].set_title(
                    'coordinate {:d}'.format(j + 1), fontsize=8)

            if hasattr(_q, 'get_domain_breakpoints'):
                t_bp = _q.get_domain_breakpoints()
                for time_i in t_bp:
                    axis[i, j].axvline(time_i, alpha=0.1, color='red')
                if _wp is not None and i == 0:
                    axis[i, j].plot(t_bp, _wp[:, j], 'b*')

    plt.subplots_adjust(
        left=0.025,
        bottom=0.05,
        right=0.975,
        top=0.95,
        wspace=0.25,
        hspace=0.15)

    if _show:
        plt.show()
