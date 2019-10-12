#include <cartint/gaussianquadrature.h>

#include <random>
#include <cmath>
#include <memory>

void GaussianQuadrature::setWeightsAndPoints(const size_t& size,
        Eigen::VectorXd& weights, Eigen::VectorXd& x) {
    /* set weights and points in hermite quadrature */
    unsigned int s = (size<2 ? 2 : size);
    double a = pow(2, s-1) * boost::math::factorial<double>(s) * sqrt(M_PI) /
        (s*s);
    setPoints(s, x);
    weights = Eigen::VectorXd(s);
    for (unsigned int i = 0; i < s; ++i) {
        double h = H(x(i), s-1);
        weights(i) = a / (h*h);
    } // end fori
} // end function setWeight

void GaussianQuadrature::setPoints(const size_t& n, Eigen::VectorXd& x) {
    /* set points in hermite quadrature */

    // make initial guess of points
    x = Eigen::VectorXd(n);
    double nInv = 1./n;
    std::unique_ptr<double[]> pair{new double[2]{}};
    for (unsigned int i = 0; i < n; ++i) {
        generateRandom(pair.get());
        x(i) = pair[0]*n + pair[1]*nInv;
    } // end fori

    // find roots
    findRootsHermite(n, x);
} // end function setPoints

void GaussianQuadrature::findRootsHermite(const size_t& level, Eigen::VectorXd&
        x) {
    /* find roots of hermite polynomial of order level with Aberth-Erlich
     * method */
    size_t size = x.size();
    double tmpSum = 0;
    double dx = 0.001;
    double dxHalf = dx/2.;
    Eigen::VectorXd prevP = x;

    // run algorithm
    bool diffTest = true;
    while(diffTest) {
        for (unsigned int i = 0; i < size; ++i) {
            for (unsigned int j = 0; j < size; ++j) {
                if(i!=j) {
                    tmpSum += 1/(x(i)-x(j));
                } // end if j not i
            } // end forj
            prevP(i) += 1/(tmpSum - (H(x(i)+dxHalf,level) -
                        H(x(i)-dxHalf,level)) / (dx*H(x(i),level)));
            tmpSum = 0;
        } // end fori
        for (unsigned int k = 0; k < prevP.size(); ++k) {
            if(fabs(x(k)-prevP(k)) <= 1e-14) {
                /* check that all roots are within convergence */
                prevP(k) = (fabs(prevP(k))<=1e-14) ? 0 : prevP(k);
                diffTest = false;
            } else {
                /* break and rerun if convergence is not reached */
                diffTest = true;
                break;
            } // end ifelse
        } // end fork
        x = prevP;
    } // end while diffTest
} // end function findRootsHermite

void GaussianQuadrature::generateRandom(double* pair) {
    /* generate pair of random numbers based on Box-Muller transform */

    // initialize RNG
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> randomReal(-1,1);

    // find suitable pair of numbers
    double u, v, s;
    while (true) {
        u = randomReal(generator);
        v = randomReal(generator);
        s = u*u + v*v;
        if((s>0) && (s<1)) {
            break;
        } // endif s
    } // end while true
    double S = sqrt(-2*log(s)/s);

    // set numbers in pair array
    pair[0] = S*u;
    pair[1] = S*v;
} // end function generateRandom
