import gsplines
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from .LinesContainer import VerticalLines, CurveVsTime
import matplotlib as mpl


def fix_to_lines(_ax):
    lines = _ax.get_lines()

    xmax = max([np.max(line.get_xdata()) for line in lines])
    xmin = min([np.min(line.get_xdata()) for line in lines])

    ymax = max([np.max(line.get_ydata()) for line in lines])
    ymin = min([np.min(line.get_ydata()) for line in lines])

    _ax.set_xlim(1.01*xmin, 1.01*xmax)
    _ax.set_ylim(1.01*ymin, 1.01*ymax)


def plot_derivatives_in_axes(_q: gsplines.functions.FunctionBase,
                             _axis, _dt: float = 0.01,
                             _up_to_deriv: int = 3, color='blue'):

    for i in range(_up_to_deriv+1):
        curve = _q.deriv(i)
        cvt = CurveVsTime(curve.get_codom_dim(), color=color)
        cvt.associate_axis(_axis[i, :])
        cvt.update(curve, _dt)
        for j in range(_q.get_codom_dim()):
            fix_to_lines(_axis[i, j])

        if hasattr(curve, 'get_domain_breakpoints'):
            dbp = VerticalLines(
                _q.get_number_of_intervals()+1, curve.get_codom_dim(),
                linestyle='--', alpha=0.4, color=color)
            dbp.associate_axis(_axis[i, :])
            dbp.update(curve.get_domain_breakpoints())

    if hasattr(_q, 'get_waypoints'):
        wp = _q.get_waypoints()
        cvt = CurveVsTime(wp.shape[1], marker='.',
                          linestyle='', markersize=10, color=color)
        cvt.associate_axis(_axis[0, :])
        cvt.update_data(_q.get_domain_breakpoints(), wp)


def plot(_q: gsplines.functions.FunctionBase, _up_to_deriv: int = 3,
         _dt: float = 0.01, _wp=None, _title='', _show=True):

    fig, axis = plt.subplots(_up_to_deriv+1, _q.get_codom_dim())

    plot_derivatives_in_axes(_q, axis)

    plt.show()


def plot_compare(_q: List, _colors: List, _up_to_deriv=3, _dt=0.1, _wp=None, _title='', _show=True):

    fig, axis = plt.subplots(_up_to_deriv + 1, _q[0].get_codom_dim())

    for curve in _q:
        plot_derivatives_in_axes(curve, axis)

    plt.subplots_adjust(
        left=0.025,
        bottom=0.05,
        right=0.975,
        top=0.95,
        wspace=0.25,
        hspace=0.15)

    if _show:
        plt.show()
