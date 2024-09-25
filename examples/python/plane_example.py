"""
This examples shows how to use the interpolator to joint waypoints with a
gspline.
- Input to the interpolator
    - interval lengths an array with doubles representing the lenth of the time
      interval taken to go from one waypoint to the next
    - waypoints: actual waypoints to interpolate a a matrix which rows are the
      waypoints
    - basis: The basis to compute the interpolation. e.g. Lagrange polynomials

The resulting gsplines is a curve [0, sum(interval lengths)] -> R^n

Note tha the the execution time of the trajectory is sum(interval lengths)
"""
import numpy as np
try:
    import gsplines
    import gsplines.plot
    import gsplines.optimization
except ImportError:
    import pathlib
    import sys
    cwd = pathlib.Path(__file__).parent.absolute()
    setup_script = pathlib.Path(
        cwd, '../../', 'setup_python_env.py')

    if not setup_script.exists():
        print('Import Error: Seems that you haven\'t compile the project')
        sys.exit(1)

    exec(setup_script.read_text())
    try:
        import gsplines
        import gsplines.plot
        import gsplines.optimization

    except ImportError:
        print('Import Error: Seems that you failed to compile properly')
        sys.exit(1)


def main():
    numberOfWaypoints = 4
    dimensionOfAmbientSpace = 2

    waypoints = np.random.rand(numberOfWaypoints, dimensionOfAmbientSpace)

    c1 = gsplines.optimization.broken_lines_path(waypoints)
    c2 = gsplines.optimization.minimum_acceleration_path(waypoints)
    c3 = gsplines.optimization.minimum_jerk_path(waypoints)
    c4 = gsplines.optimization.minimum_snap_path(waypoints)
    c5 = gsplines.optimization.minimum_crackle_path(waypoints)

    gsplines.plot.plot2d_compare([c1, c2, c3, c4, c5], [
                                 'green', 'blue', 'magenta', 'red', 'black'],
                                 ['min vel', 'min acceleration', 'min jerk',
                                  'min snap', 'min crackle'])


if __name__ == "__main__":
    main()
