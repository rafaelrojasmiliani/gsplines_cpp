#ifndef INTEGRAL
#define INTEGRAL
#include <gsplines/Functions/FunctionBase.hpp>
#include <utility>
#include <vector>
namespace gsplines {
namespace functional_analysis {

double integral(const functions::FunctionBase &_in, std::size_t _n_glp = 10,
                std::size_t _n_int = 1);

} // namespace functional_analysis
} // namespace gsplines
#endif /* ifndef INTEGR */
