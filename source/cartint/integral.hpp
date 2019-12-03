#ifndef INTEGRAL_HPP
#define INTEGRAL_HPP

#include <cstddef>
#include <cassert>
#include <type_traits>

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <cartint/hermite.h>
#include <cartint/hexpander.h>

namespace cartint {

template<
    std::size_t Dim,
    typename potential_element = struct A(),
    typename = std::enable_if_t<(Dim == 2) || (Dim == 3)>
>
class Integral {
private:
    Eigen::MatrixXi levels;
    Eigen::MatrixXd gaussian_functions;
    
public:
    Integral () {}

    Integral (
            const Eigen::Ref<const Eigen::MatrixXi> angular_moment_marix,
            const Eigen::Ref<const Eigen::MatrixXd> gaussian_matrix
    ) : levels(angular_moment_marix), gaussian_functions(gaussian_matrix) {
        assert(angular_moment_marix.cols() == Dim && "Angular moment matrix dimension is wrong!");
        assert(gaussian_matrix.cols() == Dim && "Gaussian matrix dimension is wrong!");
    } // end constructor

    template<typename = std::enable_if_t<Dim == 2>>
    double overlap_element(const std::size_t& i, const std::size_t& j) {
    }
    
    template<typename = std::enable_if_t<Dim == 2>>
    double coulomb_element(
            const std::size_t& i,
            const std::size_t& j,
            const std::size_t& k,
            const std::size_t& l
    ) {
    }
    
    template<typename = std::enable_if_t<Dim == 2>>
    double laplacian_element(const std::size_t& i, const std::size_t& j) {
    }

    virtual ~Integral () { /* deconstructor */ }
}; // end class Integral

using Integral2D = Integral<2>;
using Integral3D = Integral<3>;

} // end namespace cartint

#endif /* INTEGRAL_HPP */

