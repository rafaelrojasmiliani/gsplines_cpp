try:
    import inspect
    import sys
    from .pygsplines import *
    # ---
    from .plot import *
except ImportError:
    import pygsplines
    from pygsplines import *
    from .plot import *


submodules = inspect.getmembers(pygsplines, inspect.ismodule)
for module_info in submodules:
    print(module_info)
    sys.modules['gsplines.' + module_info[0]] = module_info[1]
