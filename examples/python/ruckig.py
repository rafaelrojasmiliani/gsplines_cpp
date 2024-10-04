

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
    import gsplines.ruckig
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
        import gsplines.ruckig

    except ImportError:
        print('Import Error: Seems that you failed to compile properly')
        sys.exit(1)


def main():
    dimensionOfAmbientSpace = 1

    waypoints = [[0], [2]]

    # https://forum.universal-robots.com/t/maximum-axis-speed-acceleration/13338
    max_acc = 800*np.pi/180
    max_vel = 360*np.pi/180

    max_vel = dimensionOfAmbientSpace*[max_vel]
    max_acc = dimensionOfAmbientSpace*[max_acc]
    max_jerk = dimensionOfAmbientSpace*[400.5]

    r = gsplines.ruckig.interpolator(waypoints, max_vel, max_acc, max_jerk)
    z = gsplines.optimization.minimum_jerk_path(
        waypoints).linear_scaling_new_execution_time(r.get_domain()[1])

    print('r codom dim', r.get_codom_dim())
    if r is None:
        print('error')
        return
    print(r.get_domain())
    # gsplines.plot.plot_compare([r, z], ['red', 'magenta'])

    z = gsplines.optimization.rojas_path(
        waypoints, 0.3).linear_scaling_new_execution_time(r.get_domain()[1])
    gsplines.plot.plot_compare([r, z], ['red', 'magenta'], _dt=0.01)


if __name__ == "__main__":
    main()
