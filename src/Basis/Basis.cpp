#include <chrono>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/Basis/BasisLagrange.hpp>
#include <gsplines/Basis/BasisLegendre.hpp>
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

const Eigen::SparseMatrix<double, Eigen::RowMajor> &Basis::continuity_matrix(
    std::size_t _number_of_intervals, std::size_t _codom_dim,
    std::size_t _deriv_order,
    Eigen::Ref<const Eigen::VectorXd> _interval_lengths) const {

  if (not(continuity_matrix_buff_.count(_number_of_intervals) and
          continuity_matrix_buff_[_number_of_intervals].count(_codom_dim) and
          continuity_matrix_buff_[_number_of_intervals][_codom_dim].count(
              _deriv_order))) {

    // for each derivative degree, this fills
    // (_number_of_intervals-1)*_codom_dim rows

    Eigen::SparseMatrix<double, Eigen::RowMajor> result(
        (_number_of_intervals - 1) * _codom_dim * (_deriv_order + 1),
        _number_of_intervals * _codom_dim * get_dim());

    Eigen::MatrixXd left_buffer(_deriv_order + 1, get_dim());
    Eigen::MatrixXd right_buffer(_deriv_order + 1, get_dim());

    for (std::size_t der = 0; der <= _deriv_order; der++) {
      eval_derivative_on_window(-1.0, 2.0, der, left_buffer.row(der));
      eval_derivative_on_window(1.0, 2.0, der, right_buffer.row(der));
    }

    std::size_t i0, j0;

    for (std::size_t der_coor = 0; der_coor <= _deriv_order; der_coor++) {

      for (std::size_t interval_coor = 0;
           interval_coor < _number_of_intervals - 1; interval_coor++) {

        i0 = _codom_dim * (_number_of_intervals - 1) * der_coor +
             interval_coor * _codom_dim;
        // fill components relative to the rhs value od the interval
        j0 = interval_coor * get_dim() * _codom_dim;
        for (std::size_t codom_coor = 0; codom_coor < _codom_dim;
             codom_coor++) {

          for (std::size_t basis_coor = 0; basis_coor < get_dim();
               basis_coor++) {

            result.insert(i0 + codom_coor,
                          j0 + basis_coor * 1.0 + get_dim() * codom_coor) =
                right_buffer(der_coor, basis_coor);
          }
        }
        // fill components relative to the rhs value od the interval
        j0 = (interval_coor + 1) * get_dim() * _codom_dim;
        for (std::size_t codom_coor = 0; codom_coor < _codom_dim;
             codom_coor++) {

          for (std::size_t basis_coor = 0; basis_coor < get_dim();
               basis_coor++) {

            result.insert(i0 + codom_coor,
                          j0 + basis_coor * 1.0 + get_dim() * codom_coor) =
                -left_buffer(der_coor, basis_coor);
          }
        }
      }
    }

    // result.makeCompressed();
    continuity_matrix_buff_[_number_of_intervals][_codom_dim][_deriv_order] =
        result;

    continuity_matrix_dynamic_buff_[_number_of_intervals][_codom_dim]
                                   [_deriv_order] = std::move(result);
  }
  /*
    continuity_matrix_dynamic_buff_
        [_number_of_intervals][_codom_dim][_deriv_order] =
            continuity_matrix_buff_[_number_of_intervals][_codom_dim]
                                   [_deriv_order];
                                   */

  Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_dest =
      continuity_matrix_dynamic_buff_[_number_of_intervals][_codom_dim]
                                     [_deriv_order];
  Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_src =
      continuity_matrix_buff_[_number_of_intervals][_codom_dim][_deriv_order];

  // in our case all the rows of the matrix have at more than one non-zero cell
  // By this reason the outer index coincieds with the cols
  for (long k = _codom_dim * (_number_of_intervals - 1);
       k < mat_dest.outerSize(); ++k) {

    std::size_t deg = k / (_codom_dim * (_number_of_intervals - 1));

    Eigen::VectorXd deriv_factor =
        Eigen::pow(2.0 / _interval_lengths.array(), deg);

    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_dest(
        mat_dest, k);
    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_src(mat_src,
                                                                       k);
    for (/*already initialized*/; it_dest; ++it_dest, ++it_src) {
      std::size_t interval = it_dest.col() / (get_dim() * _codom_dim);
      it_dest.valueRef() = it_src.value() * deriv_factor(interval);
      assert(it_dest.col() == it_dest.index());
      assert(it_dest.row() == k);
      /*
      std::cout << "row " << it_dest.row() << " col " << it_dest.col() << "
      index "
                << it_dest.index() << " outer idx " << k << " interval " <<
      interval
                << " deg " << deg << " val " << it_dest.value() << std::endl;
                */
    }
  }

  return continuity_matrix_dynamic_buff_[_number_of_intervals][_codom_dim]
                                        [_deriv_order];
}

const Eigen::SparseMatrix<double, Eigen::RowMajor> &
Basis::gspline_derivative_matrix(
    std::size_t _number_of_intervals, std::size_t _codom_dim,
    std::size_t _deriv_order,
    Eigen::Ref<const Eigen::VectorXd> _interval_lengths) const {

  std::size_t matrix_size = _number_of_intervals * _codom_dim * get_dim();

  if (not(derivative_matrix_buff_.count(_number_of_intervals) and
          derivative_matrix_buff_[_number_of_intervals].count(_codom_dim) and
          derivative_matrix_buff_[_number_of_intervals][_codom_dim].count(
              _deriv_order))) {

    Eigen::SparseMatrix<double, Eigen::RowMajor> result(matrix_size,
                                                        matrix_size);

    const Eigen::MatrixXd &dmat = get_derivative_matrix_block(_deriv_order);
    for (std::size_t block = 0; block < _number_of_intervals * _codom_dim;
         block++) {
      std::size_t i0 = block * get_dim();
      for (std::size_t i = 0; i < get_dim(); i++) {
        for (std::size_t j = 0; j < get_dim(); j++) {
          result.insert(i0 + i, i0 + j) = dmat(i, j);
        }
      }
    }

    derivative_matrix_buff_[_number_of_intervals][_codom_dim][_deriv_order] =
        result;

    derivative_matrix_dynamic_buff_[_number_of_intervals][_codom_dim]
                                   [_deriv_order] = std::move(result);
    // std::chrono::steady_clock::time_point end =
    //    std::chrono::steady_clock::now();
  }

  Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_dest =
      derivative_matrix_dynamic_buff_[_number_of_intervals][_codom_dim]
                                     [_deriv_order];

  Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_src =
      derivative_matrix_buff_[_number_of_intervals][_codom_dim][_deriv_order];

  Eigen::VectorXd deriv_factor =
      Eigen::pow(2.0 / _interval_lengths.array(), _deriv_order);

  for (long k = 0; k < mat_dest.outerSize(); ++k) {

    // std::size_t deg = k / (_codom_dim * (_number_of_intervals - 1));

    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_dest(
        mat_dest, k);

    Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_src(mat_src,
                                                                       k);
    for (/*already initialized*/; it_dest; ++it_dest, ++it_src) {
      std::size_t interval = it_dest.col() / (get_dim() * _codom_dim);
      it_dest.valueRef() = it_src.value() * deriv_factor(interval);
    }
  }

  return derivative_matrix_dynamic_buff_[_number_of_intervals][_codom_dim]
                                        [_deriv_order];
}
} // namespace basis
} // namespace gsplines
