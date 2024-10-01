"""
Test use of FunctionExpression bindings
"""
import unittest
import numpy as np

try:
    import gsplines
except ImportError:
    import pathlib
    import sys
    cwd = pathlib.Path(__file__).parent.absolute()
    setup_script = pathlib.Path(
        cwd, '../', 'setup_python_env.py')

    if not setup_script.exists():
        print('Import Error: Seems that you haven\'t compiled the project')
        sys.exit(1)

    exec(setup_script.read_text())
    try:
        import gsplines

    except ImportError:
        print('Import Error: Seems that you failed to compile properly')
        sys.exit(1)


import gsplines.serialization


class MyTest(unittest.TestCase):
    """ Test elemental function and operations"""

    def __init__(self, *args, **kwargs):
        unittest.TestCase.__init__(self, *args, **kwargs)
        pass

    def test_basis(self):
        basis_ground = gsplines.basis.get_basis('legendre', 4, [])
        my_dict = gsplines.serialization.basis_to_dict(basis_ground)
        basis_test = gsplines.serialization.dict_to_basis(my_dict)
        self.assertEqual(basis_ground, basis_test)

    def test_gsplines(self):

        wp = np.random.rand(4, 6)
        gspline_ground = gsplines.optimization.broken_lines_path(wp)
        my_dict = gsplines.serialization.gspline_to_dict(gspline_ground)
        gspline_test = gsplines.serialization.dict_to_gspline(my_dict)
        self.assertEqual(gspline_ground, gspline_test)


if __name__ == '__main__':
    unittest.main()
