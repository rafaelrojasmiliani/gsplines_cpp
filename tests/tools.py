
import os
import unittest
import functools
import traceback
import sys
import pdb

import numpy as np
from numpy.random import rand
from numpy.linalg import norm
from numpy import pi, array
from numpy.polynomial.chebyshev import Chebyshev


def debug_on(*exceptions):
    ''' Decorator for entering in debug mode after exceptions '''
    if not exceptions:
        exceptions = (Exception, )

    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            try:
                return f(*args, **kwargs)
            except exceptions:
                info = sys.exc_info()
                traceback.print_exception(*info)
                pdb.post_mortem(info[2])
                sys.exit(1)

        return wrapper

    return decorator
