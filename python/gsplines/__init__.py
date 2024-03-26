import sys
import inspect
try:
    from .pygsplines import *  # noqa
except ImportError:
    import pygsplines
    from pygsplines import *  # noqa


submodules = inspect.getmembers(pygsplines, inspect.ismodule)
for module_info in submodules:
    sys.modules['gsplines.' + module_info[0]] = module_info[1]
