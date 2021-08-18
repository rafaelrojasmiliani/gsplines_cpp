#include <gsplines/Basis.hpp>
#include <gsplines/BasisLagrange.hpp>
#include <gsplines/BasisLegendre.hpp>
#include <iostream>
#include <math.h>
namespace gsplines {

namespace basis {

std::unique_ptr<Basis> string_to_basis(const std::string &_basis_name) {
  std::string::size_type pos = _basis_name.find("_");
  if (_basis_name.find("legendre") != std::string::npos)
    return std::make_unique<BasisLegendre>(
        std::stoul(_basis_name.substr(pos + 1)));

  throw std::invalid_argument("The basis" + _basis_name + " is unknwon");
  return nullptr;
}
} // namespace basis
} // namespace gsplines
