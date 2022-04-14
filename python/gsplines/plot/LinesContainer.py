from matplotlib.lines import Line2D
import gsplines
import numpy as np
import copy


class VerticalLines:
    def __init__(self, _number, _codom_dim, **args):
        self.lines_ = np.array(
            [[Line2D([], [], **args) for i in range(_number)]
                for j in range(_codom_dim)])
        self.axes_ = None

    def associate_axis(self, _axes):
        if self.lines_.shape[:-1] != _axes.shape:
            raise ValueError(" Axis and lines must have same dimensions")
        self.axes_ = _axes
        for i in range(self.lines_.shape[0]):
            for j, line in enumerate(self.lines_[i]):
                _axes[i].add_line(line)

    def update(self, _bplist: np.ndarray):
        for i in range(self.lines_.shape[0]):
            for k, bp in enumerate(_bplist):
                lower, upper = self.axes_[i].get_ylim()
                self.lines_[i, k].set_data([bp, bp], [lower, upper])


class CurveVsTime:
    """ Container of lines """

    def __init__(self, _codom_dim, **args):

        self.lines_ = np.array([Line2D([], [], **args)
                                for i in range(_codom_dim)])

    def associate_axis(self, _axes):
        """ call axis.add_line"""
        if self.lines_.shape != _axes.shape:
            raise ValueError(
                " Axis and lines must have same dimensions",
                self.lines_.shape, _axes.shape)
        self.axes_ = _axes
        for i in range(self.lines_.shape[0]):
            _axes[i].add_line(self.lines_[i])

    def update(self, _curve, _dt):
        """ Update data """
        time_spam = np.arange(*_curve.get_domain(), _dt)
        points_to_plot = _curve(time_spam)
        for j, line in enumerate(self.lines_):
            line.set_data(time_spam, points_to_plot[:, j])

    def update_data(self, _time_spam, _points):
        """ Update data """
        for j, line in enumerate(self.lines_):
            line.set_data(_time_spam, _points[:, j])
