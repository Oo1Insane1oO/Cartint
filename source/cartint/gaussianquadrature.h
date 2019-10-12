#ifndef GAUSSIANQUADRATURE_H
#define GAUSSIANQUADRATURE_H

#include <cartint/hermite.h>

#include <boost/math/special_functions/factorials.hpp>
#include <Eigen/Dense>
#include <iostream>

class GaussianQuadrature {
    private:
        GaussianQuadrature() {};
        virtual ~GaussianQuadrature() {};

        static void setWeightsAndPoints(const size_t&, Eigen::VectorXd&,
                Eigen::VectorXd&);
        static void setPoints(const size_t&, Eigen::VectorXd&);
        static void generateRandom(double*);
        static void findRootsHermite(const size_t&, Eigen::VectorXd&);

    public:
        template<typename U, typename F, typename... Args> static inline double
            gaussHermiteQuad(const size_t n, const U& obj, F f, Args... args) {
            /* calculate integral of f using Gauss-Hermite Quadrature */
            Eigen::VectorXd points, weights;
            
            /* set points, weights and normalization constant */
            setWeightsAndPoints(n, weights, points);
            double A = 1./(pow(2,n) * boost::math::factorial<double>(n) *
                    sqrt(M_PI));

            double sum = 0;
            for (unsigned int i = 0; i < points.size(); ++i) {
                double hi = H(points(i), i);
                sum += weights(i) * A * hi*hi * (obj->*f)(points(i), args...);
            } // end fori

            return sum;
        } // end function gaussHermiteQuad
        
        template<typename U, typename F, typename... Args> static inline double
            gaussChebyshevQuad(const unsigned int n, const U& obj, F f, Args...
                    args) {
            /* calculate integral of f using Gauss-Chebyshev Quadrature */
            unsigned int s = (n<2 ? 2 : n); // n > 2 is assumed
            double pis = M_PI/s; // all weights are pi/n
            double sum = 0.0;
            #pragma omp parallel for reduction(+:sum)
            for (unsigned int i = 0; i < s; ++i) {
                sum += (obj->*f)(cos(pis * (i-0.5)), args...);
            } // end fori

            return sum * pis;
        } // end function gaussChebyshevQuad
};

#endif /* GAUSSIANQUADRATURE_H */
