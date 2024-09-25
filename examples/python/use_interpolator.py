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

    except ImportError:
        print('Import Error: Seems that you failed to compile properly')
        sys.exit(1)


def main():
    numberOfWaypoints = 4
    dimensionOfAmbientSpace = 4

    intervalLengths = np.random.rand(numberOfWaypoints-1)+0.5

    basis = gsplines.basis.Basis0101(0.999)
    waypoints = np.random.rand(numberOfWaypoints, dimensionOfAmbientSpace)

    curve = gsplines.interpolate(intervalLengths, waypoints, basis)

    gsplines.plot.plot(curve)

    basis = gsplines.basis.BasisLegendre(6)

    curve2 = gsplines.interpolate(intervalLengths, waypoints, basis)

    gsplines.plot.plot_compare([curve, curve2], ['blue', 'magenta'])


if __name__ == "__main__":
    main()
