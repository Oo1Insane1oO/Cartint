#include <cartint/hexpander.h>

#include <cartint/gaussianquadrature.h>

#include <cmath>
#include <iostream>

Hexpander::Hexpander() {
    /* default constructor */
} // end constructor

Hexpander::~Hexpander() {
} // end deconstructor

void Hexpander::setCoefficients(unsigned int iMax, unsigned int jMax, double a,
        double b, double Qx) {
    /* find coefficients for overlap distributions (product between two
     * hermite-gaussians */
    iMp1 = iMax+1;
    jMp1 = jMax+1;
    nMp1 = iMp1+jMp1+1;

    double p = a + b;
    double q = a*b/p;

    coefficients = Eigen::ArrayXd::Zero(iMp1*jMp1*nMp1);
    coefficients(0) = exp(-q*Qx*Qx);
    for (unsigned int i = 0; i < jMp1; ++i) {
        for (unsigned int t = 0; t < nMp1; t++) {
            if ((i==0) && (t==0)) {
                continue;
            } // end if
            int im = i - 1;
            int tm = t - 1;
            int tp = t + 1;

            double Em = 0;
            if(checkIndices(0, im, tm)) {
                Em = coefficients(cidx(0,im,tm));
            } // end if
            double Ed = 0;
            if(checkIndices(0, im, t)) {
                Ed = coefficients(cidx(0,im,t));
            } // end if
            double Ep = 0;
            if(checkIndices(0, im, tp)) {
                Ep = coefficients(cidx(0,im,tp));
            } // end if
            coefficients(cidx(0,i,t)) = 0.5/p*Em + q*Qx/b*Ed + tp*Ep;
        } // end fort
    } // end fori

    for (unsigned int i = 1; i < iMp1; ++i) {
        for (unsigned int j = 0; j < jMp1; ++j) {
            for (unsigned int t = 0; t < nMp1; ++t) {
                int im = i - 1;
                int tm = t - 1;
                int tp = t + 1;

                double Em = 0;
                if(checkIndices(im, j, tm)) {
                    Em = coefficients(cidx(im,j,tm));
                } // end if
                double Ed = 0;
                if(checkIndices(im, j, t)) {
                    Ed = coefficients(cidx(im,j,t));
                } // end if
                double Ep = 0;
                if(checkIndices(im, j, tp)) {
                    Ep = coefficients(cidx(im,j,tp));
                } // end if
                coefficients(cidx(i,j,t)) = 0.5/p*Em + q*Qx/a*Ed + tp*Ep;
            } // end fort
        } // end forj
    } // end fori
} // end function setCoefficients

void Hexpander::setAuxiliary2D(unsigned int xMax, unsigned int yMax, double a,
        double b, double c, double d, const Eigen::VectorXd& PQ) {
    /* setup integral elements in 2D */
    x2DMp1 = xMax + 1;
    y2DMp1 = yMax + 1;
    n2DMp1 = x2DMp1 + y2DMp1 + 1;

    double p = (a+c)*(b+d) / (a+c+b+d);

    integrals2D = Eigen::ArrayXd::Zero(x2DMp1 * y2DMp1 * n2DMp1);

    double powVal = 1;
    for (unsigned int n = 0; n < n2DMp1; ++n) {
        /* calculate initial integrals */
        integrals2D(I2Didx(0,0,n)) = powVal *
            GaussianQuadrature::gaussChebyshevQuad(50, this,
                    &Hexpander::modifiedIntegrand, n, p*PQ.squaredNorm());
        powVal *= -2*p;
    } // end forn

    for (unsigned int ts = 1; ts < n2DMp1; ++ts) {
        for (unsigned int n = 0; n < n2DMp1-ts; ++n) {
            for (unsigned int i = 0; i < x2DMp1; ++i) {
                for (unsigned int j = 0; j < y2DMp1; ++j) {
                    if ((i+j != ts) || (i+j == 0)) {
                        /* out of bounds */
                        continue;
                    } // end if
                    unsigned int ijMax = (i > j ? i : j);
                    int i2 = i;
                    int j2 = j;
                    int i3 = i;
                    int j3 = j;
                    int factor = i;
                    double XPQ = PQ(0);
                    if (ijMax == i) {
                        i2 = i-2;
                        i3 = i-1;
                        factor = i3;
                        XPQ = PQ(0);
                    } else {
                        j2 = j-2;
                        j3 = j-1;
                        factor = j3;
                        XPQ = PQ(1);
                    } // end ifelse
                    double currValue = 0;
                    if (i2 >= 0 && j2 >= 0) {
                        currValue += factor * integrals2D(I2Didx(i2,j2,n+1));
                    } // end if
                    if (i3 >= 0 && j3 >= 0) {
                        currValue += XPQ * integrals2D(I2Didx(i3,j3,n+1));
                    } // end if
                    integrals2D(I2Didx(i,j,n)) = currValue;
                } // end forj
            } // end fori
        } // end forn
    } // end forts
} // end function setAuxiliary2D

void Hexpander::setAuxiliary3D(unsigned int xMax, unsigned int yMax, unsigned
        int zMax, double a, double b, double c, double d, const
        Eigen::VectorXd& PQ) {
    /* setup integral elements in 3D */
    x3DMp1 = xMax + 1;
    y3DMp1 = yMax + 1;
    z3DMp1 = zMax + 1;
    n3DMp1 = xMax + yMax + zMax + 1;

    double p = (a+c)*(b+d) / (a+c+b+d);

    integrals3D = Eigen::ArrayXd::Zero(x3DMp1 * y3DMp1 * z3DMp1 * n3DMp1);

    double powVal = 1;
    for (unsigned int n = 0; n < n3DMp1; ++n) {
        integrals3D(I3Didx(0,0,0,n)) = powVal *
            GaussianQuadrature::gaussChebyshevQuad(50, this,
                    &Hexpander::boysIntegrand, n, p*PQ.squaredNorm());
        powVal *= -2*p;
    } // end forn

    for (unsigned int ts = 1; ts < n3DMp1; ++ts) {
        for (unsigned int n = 0; n < n3DMp1 - ts; ++n) {
            for (unsigned int i = 0; i < x3DMp1; ++i) {
                for (unsigned int j = 0; j < y3DMp1; ++j) {
                    for (unsigned int k = 0; k < z3DMp1; ++k) {
                        if (i+j+k != ts || i+j+k == 0) {
                            /* out of bounds */
                            continue;
                        } // end if
                        unsigned int ijkMax = (i>j ? (i>k ? i : k) : (j>k ? j :
                                    k));
                        int i2 = i;
                        int j2 = j;
                        int k2 = k;
                        int i3 = i;
                        int j3 = j;
                        int k3 = k;
                        int factor = i;
                        double XPQ = PQ(0);
                        if (ijkMax == i) {
                            i2 = i - 2;
                            i3 = i - 1;
                            factor = i3;
                            XPQ = PQ(0);
                        } else if (ijkMax == j) {
                            j2 = j - 2;
                            j3 = j - 1;
                            factor = j3;
                            XPQ = PQ(1);
                        } else {
                            k2 = k - 2;
                            k3 = k - 1;
                            factor = k3;
                            XPQ = PQ(2);
                        } // end ifeifelse
                        double  currValue = 0;
                        if (i2 >= 0 && j2 >= 0 && k2 >= 0) {
                            currValue += factor *
                                integrals3D(I3Didx(i2,j2,k2,n+1));
                        } // end if
                        if (i3 >= 0 && j3 >= 0 && k3 >= 0) {
                            currValue += XPQ *
                                integrals3D(I3Didx(i3,j3,k3,n+1));
                        } // end if
                        integrals3D(I3Didx(i,j,k,n)) = currValue;
                    } // end fork
                } // end forj
            } // end fori
        } // end forn
    } // end forts
} // en function setAuxiliary3D

bool Hexpander::checkIndices(const int& ia, const int& ib, const int& t) {
    if ((t < 0) || (t > (ia+ib)) || (ia < 0) || (ib < 0)) {
        return false;
    } else {
        return true;
    } // end if
} // end function checkIndices

double Hexpander::boysIntegrand(double u, const unsigned int& n, const double&
        pRR) {
    /* Integrand of Boys function (using Gauss-Chebyshev method) */
    double powu = 1;
    double uu = -u*u;
    for (unsigned int i = 0; i < 2*n; ++i) {
        powu *= u;
    } // end fori
    return powu * sqrt(1 + uu) * exp(uu*pRR);
} // end function boysIntegrand

double Hexpander::modifiedIntegrand(double u, const unsigned int& n, const
        double& pRR) {
    /* integrand assuming Gauss-Chebyshev method is used (the 1/sqrt(1-u^2)
     * part is absorbed) */
    double powu = 1;
    for (unsigned int i = 0; i < 2*n; ++i) {
        powu *= u;
    } // end fori
    return powu * exp(-u*u*pRR); // * 1/sqrt(1-u*u);
} // end function modifiedIntegrand

const double& Hexpander::coeff(const unsigned int& i, const unsigned int& j,
        const unsigned int& t) const {
    /* return coefficient E^{ij}_t */
    return coefficients(cidx(i,j,t));
} // end coeff

const double& Hexpander::auxiliary2D(const unsigned int& n, const unsigned int&
        i, const unsigned int& j) const {
    /* return integral R^{ij}_n */
    return integrals2D(I2Didx(i,j,n));
} // end function integral2D

const double& Hexpander::auxiliary3D(const unsigned int& n, const unsigned int&
        i, const unsigned int& j, const unsigned int& k) const {
    /* return integral R^{ij}_n */
    return integrals3D(I3Didx(i,j,k,n));
} // end function integral2D
