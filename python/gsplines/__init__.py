from .pygsplines import *
import sys
import inspect

submodules = inspect.getmembers(pygsplines, inspect.ismodule)
for module_info in submodules:
    sys.modules['gsplines.' + module_info[0]] = module_info[1]
