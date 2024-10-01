
import gsplines
from typing import Dict
from typing import List


def basis_to_dict(_basis: gsplines.basis.Basis) -> Dict:
    result = {}
    result['name'] = _basis.get_name()
    result['dim'] = _basis.get_dim()
    result['parameters'] = list(_basis.get_parameters())
    return result


def dict_to_basis(_dict: Dict) -> gsplines.basis.Basis:
    required_keys = {'name': str, 'dim': int, 'parameters': List}
    for key, expected_type in required_keys.items():
        if not isinstance(_dict.get(key), expected_type):
            raise TypeError(
                f"Incorrect type for key '{key}'. Expected {expected_type.__name__}, got {type(_dict.get(key)).__name__}")
    return gsplines.basis.get_basis(
        _dict['name'], _dict['dim'], _dict['parameters'])


def gspline_to_dict(_gspline: gsplines.GSpline) -> Dict:
    result = {}
    dom = _gspline.get_domain()
    result['domain_left_boundary'] = dom[0]
    result['domain_right_boundary'] = dom[1]

    result['codom_dim'] = _gspline.get_codom_dim()

    result['number_of_intervals'] = _gspline.get_number_of_intervals()

    result['basis'] = basis_to_dict(_gspline.get_basis())

    result['coefficients'] = list(_gspline.get_coefficients().copy())
    result['interval_lengths'] = list(_gspline.get_interval_lengths().copy())
    return result


def dict_to_gspline(_dict: Dict) -> gsplines.GSpline:
    required_keys = {'domain_left_boundary': float,
                     'domain_right_boundary': float,
                     'codom_dim': int,
                     'number_of_intervals': int,
                     'basis': Dict,
                     'coefficients': List,
                     'interval_lengths': List}
    for key, expected_type in required_keys.items():
        if not isinstance(_dict.get(key), expected_type):
            raise TypeError(
                f"Incorrect type for key '{key}'. Expected {expected_type.__name__}, got {type(_dict.get(key)).__name__}")

    domain = (_dict['domain_left_boundary'], _dict['domain_right_boundary'])
    codom_dim = _dict['codom_dim']
    number_of_intervals = _dict['number_of_intervals']
    basis = dict_to_basis(_dict['basis'])
    coefficients = _dict['coefficients']
    interval_lengths = _dict['interval_lengths']

    return gsplines.GSpline(domain, codom_dim, number_of_intervals,
                            basis, coefficients, interval_lengths, "gspline")
