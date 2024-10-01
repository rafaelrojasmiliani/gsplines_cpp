#!/usr/bin/python3
"""
Example of how to get a curve that stops in minimum time with acceleration
bounds.
"""
import unittest
import numpy as np
import pathlib
import yaml

try:
    import gsplines
    import gsplines.serialization
    import gsplines.plot as gplot
except ImportError:

    cwd = pathlib.Path(__file__).parent.absolute()
    setup_script = pathlib.Path(
        cwd, '../../', 'setup_python_env.py')
    exec(setup_script.read_text())
    try:
        import gsplines
        import gsplines.serialization
        import gsplines.plot as gplot
    except ImportError as e:
        import sys
        print('fail', e)
        sys.exit(1)


def main():
    """
    Main example
    """
    yaml_file_path = pathlib.Path(__file__).absolute(
    ).parents[0]/'message.yaml'
    with open(yaml_file_path, 'r') as yaml_file:
        data = yaml.safe_load(yaml_file)

    original_gspline_gspline = data.get(
        'original_gspline', {}).get('gspline', {})
    gspline_gspline = data.get(
        'gspline', {}).get('gspline', {})

    gs1 = gsplines.serialization.dict_to_gspline(gspline_gspline)

    gs = gsplines.serialization.dict_to_gspline(original_gspline_gspline)
    gplot.plot(gs, _up_to_deriv=2)
    gplot.plot(gs1, _up_to_deriv=2)


#     stop_ti = 1.29

#     alpha_acc_bound = 1.98

#     optimal_parametrization = opstop.minimum_time_bounded_acceleration(
#         trj, stop_ti, alpha_acc_bound, str(model_file), 8)

#     stop_trj = trj.compose(optimal_parametrization)

#     # print(alpha)
    # gplot.plot_compare([gs, gs1], ['red', 'blue'], [
    #                    'Emergency Stop Trajectory',
    #                    'Original Trajectory'], _show=True, _up_to_deriv=2)

if __name__ == '__main__':
    main()
