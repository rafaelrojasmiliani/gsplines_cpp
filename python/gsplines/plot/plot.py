from typing import List

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from .LinesContainer import VerticalLines, CurveVsTime
import matplotlib as mpl


def fix_to_lines(_ax):
    lines = _ax.get_lines()

    xmax = max([np.max(line.get_xdata()) for line in lines])
    xmin = min([np.min(line.get_xdata()) for line in lines])

    ymax = max([np.max(line.get_ydata()) for line in lines])
    ymin = min([np.min(line.get_ydata()) for line in lines])

    xcenter = 0.5*(xmin + xmax)
    dx = xmax - xmin
    ycenter = 0.5*(ymin + ymax)
    dy = ymax - ymin

    _ax.set_xlim(xcenter - 1.02*dx/2, xcenter + 1.02*dx/2)
    _ax.set_ylim(ycenter - 1.02*dy/2, ycenter + 1.02*dy/2)


def plot_derivatives_in_axes(_q,
                             _axis: np.array, _dt: float = 0.01,
                             _up_to_deriv: int = 3, color='blue'):

    _axis = np.atleast_2d(_axis)
    for i in range(_up_to_deriv+1):
        curve = _q.deriv(i)
        cvt = CurveVsTime(curve.get_codom_dim(), color=color)
        cvt.associate_axis(_axis[i, :])
        cvt.update(curve, _dt)

        if hasattr(curve, 'get_domain_breakpoints'):
            dbp = VerticalLines(
                _q.get_number_of_intervals()+1, curve.get_codom_dim(),
                linestyle='--', alpha=0.4, color=color)
            dbp.associate_axis(_axis[i, :])
            dbp.update(curve.get_domain_breakpoints())

        for j in range(_q.get_codom_dim()):
            fix_to_lines(_axis[i, j])

    if hasattr(_q, 'get_waypoints'):
        wp = _q.get_waypoints()
        cvt = CurveVsTime(wp.shape[1], marker='.',
                          linestyle='', markersize=10, color=color)
        cvt.associate_axis(_axis[0, :])
        cvt.update_data(_q.get_domain_breakpoints(), wp)

    for j in range(0, _q.get_codom_dim()):
        _axis[0, j].set_title('coordinate {:d}'.format(j))

    for i in range(1, _up_to_deriv+1):
        _axis[i, 0].set_ylabel('derivative {:d}'.format(i))


def plot(_q, _up_to_deriv: int = 3,
         _dt: float = 0.01, _wp=None, _title='', _show=True):

    fig, axis = plt.subplots(_up_to_deriv+1, _q.get_codom_dim())

    plot_derivatives_in_axes(_q, axis, _dt=_dt, _up_to_deriv=_up_to_deriv)

    plt.subplots_adjust(
        left=0.05,
        bottom=0.05,
        right=0.975,
        top=0.95,
        wspace=0.25,
        hspace=0.15)
    if (_show):
        plt.show()


def plot_compare(_q: List, _colors: List = [], _legends: List = [],
                 _up_to_deriv=3,
                 _dt=0.1, _wp=None, _title='', _show=True):

    if not _colors:
        cmap = plt.get_cmap('viridis')
        _colors = [cmap(val) for val in np.linspace(0, 1, len(_q))]

    fig, axis = plt.subplots(_up_to_deriv + 1, _q[0].get_codom_dim())
    if _q[0].get_codom_dim() == 1:
        axis = np.atleast_2d(axis).transpose()

    for curve, color in zip(_q, _colors):
        plot_derivatives_in_axes(
            curve, axis, color=color, _up_to_deriv=_up_to_deriv, _dt=_dt)

    for j in range(0, _q[0].get_codom_dim()):
        axis[0, j].set_title('coordinate {:d}'.format(j))

    for i in range(1, _up_to_deriv+1):
        axis[i, 0].set_ylabel('derivative {:d}'.format(i))

    if len(_legends) != 0:
        for i in range(0, _up_to_deriv+1):
            for j in range(0, _q[0].get_codom_dim()):
                axis[i, j].legend(axis[i, j].get_lines(),
                                  _legends, loc='lower right', fontsize=8)
    plt.subplots_adjust(
        left=0.05,
        bottom=0.05,
        right=0.975,
        top=0.95,
        wspace=0.25,
        hspace=0.15)

    if _show:
        plt.show()


def plot2d(_q, _up_to_deriv: int = 3,
           _dt: float = 0.01, _wp=None, _title='', _show=True):

    if _q.get_codom_dim() != 2:
        return None

    fig = plt.figure(figsize=(10, 5))

    # Define the grid spec for the layout
    gs = gridspec.GridSpec(1 + _up_to_deriv, 3, width_ratios=[2, 1, 1])
    ax_big = fig.add_subplot(gs[:, 0])

    time_spam = np.arange(*_q.get_domain(), _dt)
    points_to_plot = _q(time_spam)
    ax_big.plot(points_to_plot[:, 0], points_to_plot[:, 1])
    if hasattr(_q, 'get_waypoints'):
        waypoints = _q.get_waypoints()
        ax_big.scatter(waypoints[:, 0],
                       waypoints[:, 1], color='red', s=100, marker='o',
                       label='Points')  # 'o' marker for circles

    axes = np.array([[fig.add_subplot(gs[row, col + 1])
                      for col in range(2)] for row in range(1+_up_to_deriv)])
    plot_derivatives_in_axes(_q, axes)

    plt.subplots_adjust(
        left=0.05,
        bottom=0.05,
        right=0.975,
        top=0.95,
        wspace=0.25,
        hspace=0.15)
    if (_show):
        plt.show()


def plot2d_compare(_q: List, _colors: List = [], _legends: List = [],
                   _up_to_deriv: int = 3,
                   _dt: float = 0.01, _wp=None, _title='', _show=True):

    if not _colors:
        cmap = plt.get_cmap('viridis')
        _colors = [cmap(val) for val in np.linspace(0, 1, len(_q))]

    if not all(q.get_codom_dim() == 2 for q in _q):
        return None

    fig = plt.figure(figsize=(10, 5))

    # Define the grid spec for the layout
    gs = gridspec.GridSpec(1 + _up_to_deriv, 3, width_ratios=[2, 1, 1])
    ax_big = fig.add_subplot(gs[:, 0])

    for q, color in zip(_q, _colors):
        time_spam = np.arange(*q.get_domain(), _dt)
        points_to_plot = q(time_spam)
        ax_big.plot(points_to_plot[:, 0], points_to_plot[:, 1], color=color)
        if hasattr(q, 'get_waypoints'):
            waypoints = q.get_waypoints()
            ax_big.scatter(waypoints[:, 0],
                           waypoints[:, 1], color='red', s=100, marker='o',
                           label='Points')  # 'o' marker for circles

    if len(_legends) != 0:
        ax_big.legend(ax_big.get_lines(),
                      _legends, loc='lower right', fontsize=8)
    axis = np.array([[fig.add_subplot(gs[row, col + 1])
                      for col in range(2)] for row in range(1+_up_to_deriv)])

    for curve, color in zip(_q, _colors):
        plot_derivatives_in_axes(
            curve, axis, color=color, _up_to_deriv=_up_to_deriv)

    if len(_legends) != 0:
        for i in range(0, _up_to_deriv+1):
            for j in range(0, _q[0].get_codom_dim()):
                axis[i, j].legend(axis[i, j].get_lines(),
                                  _legends, loc='lower right', fontsize=8)
    for j in range(0, _q[0].get_codom_dim()):
        axis[0, j].set_title('coordinate {:d}'.format(j))

    for i in range(1, _up_to_deriv+1):
        axis[i, 0].set_ylabel('derivative {:d}'.format(i))
    plt.subplots_adjust(
        left=0.05,
        bottom=0.05,
        right=0.975,
        top=0.95,
        wspace=0.25,
        hspace=0.15)
    if (_show):
        plt.show()
