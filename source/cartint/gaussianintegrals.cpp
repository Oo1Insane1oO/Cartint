#include <cartint/gaussianintegrals.h>
#include <cartint/hermite.h>

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>

GaussianIntegrals::GaussianIntegrals(const unsigned int dim, unsigned int cutOff, double scaling) :
    GaussianBasis(cutOff, dim) {
    m_cutOff = cutOff;
    m_dim = dim;
    expScaleFactor = scaling;
    sqrtFactor = sqrt(scaling);
    coeffs = std::make_unique<Hexpander>();
} // end constructor

GaussianIntegrals::~GaussianIntegrals() {
} // end deconstructor

std::string GaussianIntegrals::initializeParameters(double omega) {
    /* set value of oscillator frequency */
    xScale = omega;
    xScaleHalf = 0.5*xScale;
    sqrtScale1 = sqrt(pow(xScale, m_dim));
    sqrtScale = 1./sqrtScale1;
    powScale = pow(xScale, 2*m_dim);
    coulomb2DFactor = pow(M_PI/xScale, 1.5) / sqrt(2);
    coulomb3DFactor = pow(M_PI/xScale, 2.5) / sqrt(2);
    
    nSum = Eigen::ArrayXi::Zero(m_dim);

    // choose coulombElement function for 2D or 3D and set all coefficients and
    // 1D integral elements
    const Eigen::Ref<const Eigen::VectorXi> nvalues =
        GaussianBasis::Cartesian::getn();
    unsigned int pmax = nvalues(nvalues.size()-1);
    unsigned int auxMax = 2*(pmax+pmax);
    
    xScaleSqrtPow = Eigen::ArrayXd::Zero((pmax+1)*4*m_dim);
    for (unsigned int i = 0; i < xScaleSqrtPow.size(); ++i) {
        xScaleSqrtPow(i) = pow(xScale, i/2.);
    } // end fori

    static Eigen::VectorXd centerVec = Eigen::VectorXd::Constant(m_dim, 0.0);
    coeffs->setCoefficients(pmax, pmax, xScaleHalf, xScaleHalf, 0.0);
    if (m_dim == 2) {
        coeffs->setAuxiliary2D(auxMax, auxMax, xScaleHalf, xScaleHalf,
                xScaleHalf, xScaleHalf, centerVec);
        coulombFunc = &GaussianIntegrals::coulomb2D;
    } else {
        coeffs->setAuxiliary3D(auxMax, auxMax, auxMax, xScaleHalf, xScaleHalf,
                xScaleHalf, xScaleHalf, centerVec);
        coulombFunc = &GaussianIntegrals::coulomb3D;
    } // end if

    // calculate and set normalization factors to array
    setNormalizations();
    
    bool isFull = false;
    for (unsigned int i = 0; i < GaussianBasis::getMagic().size(); ++i) {
        if (m_cutOff==GaussianBasis::getMagic(i)) {
            isFull = true;
            break;
        } // end if
    } // end fori

    if (!isFull) {
        /* print message if setup is unsuccessfull */
        std::string possibleN = " ";
        for (unsigned int i = 0; i < GaussianBasis::getMagic().size(); ++i) {
            possibleN += std::to_string(GaussianBasis::getMagic(i)) + " ";
        } // end fori
        return "Slater not full, possible N:" + possibleN;
    } else {
        return "";
    } // end if
} // end function initializeParameter 

void GaussianIntegrals::setNormalizations() {
    /* calculate and set normalization factors for all basis functions */
    normalizationFactors =
        Eigen::ArrayXd::Constant(GaussianBasis::Cartesian::getNumberOfStates(),
                1.0);
    for (unsigned int i = 0; i < GaussianBasis::Cartesian::getNumberOfStates();
            ++i) {
        for (unsigned int d = 0; d < m_dim; ++d) {
            const int& n = GaussianBasis::Cartesian::getn(i,d);
//             normalizationFactors(i) *= 1.0 / ddexpr(n,n,
//                     &GaussianIntegrals::ddexprOverlap);
            normalizationFactors(i) *= 1.0 / (sqrt(M_PI) * pow(2,n) *
                    boost::math::factorial<double>(n));
        } // end ford
    } // end fori
    normalizationFactors = normalizationFactors.cwiseSqrt();
} // end function setNormalizations

const double& GaussianIntegrals::normalizationFactor(const unsigned int& n)
    const {
    /* normalization for Gauss-Hermite of order n */
    return normalizationFactors(n); 
} // end function normalizationFactor

double GaussianIntegrals::overlapd(const unsigned int& n, const unsigned int&
        m) {
    /* calculate and return <g_n|g_m> (overlap) in 1 dimension */
    int s = n + m;
    if ((s<=-1) || (s%2==1)) {
        return 0.0;
    } // end if

    return boost::math::tgamma<double>((s+1)/2.) / sqrt(xScale);
} // end function overlapd

double GaussianIntegrals::ddexpr(const int& ndd, const int& mdd,
        double(GaussianIntegrals::* f)(const int&, const int&)) {
    /* expression for sum over contracted functions */
    double sums = 0.0;
    for (int p = 0; p <= ndd; ++p) {
        for (int q = 0; q <= mdd; ++q) {
            sums += HC::coeff(ndd,p)*HC::coeff(mdd,q) * (this->*f)(p,q);
        } // end forq
    } // end forp

    return sums;
} // end function ddexpr

inline double GaussianIntegrals::ddexprOverlap(const int& p, const int& q) {
    /* expression for 1D overlap element */
    return overlapd(p,q);
} // end function ddexpr1

inline double GaussianIntegrals::ddexprLaplacian(const int& p, const int& q) {
    /* expression for 1D laplacian element */
    return xScale * (q*(q-1)*overlapd(p,q-2) - (2*q+1)*overlapd(p,q) +
            overlapd(p,q+2));
} // end function ddexpr2

inline double GaussianIntegrals::ddexprPotential(const int& p, const int& q) {
    /* expression for 1D potential element */
    return overlapd(p,q+2);
} // end function ddexpr1

inline double GaussianIntegrals::laplacianElement(const unsigned int& i, const
        unsigned int& j) {
    /* calculate and return the laplacian integral element <i|nabla|j> */
    double sums = 0.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        double tmpProdsd = 1.0;
        for (unsigned int dd = 0; dd < m_dim; ++dd) {
            const int& ndd = GaussianBasis::Cartesian::getn(i,dd);
            const int& mdd = GaussianBasis::Cartesian::getn(j,dd);
            if (dd != d) {
                tmpProdsd *= ddexpr(ndd, mdd,
                        &GaussianIntegrals::ddexprOverlap);
            } else {
                tmpProdsd *= ddexpr(ndd, mdd,
                        &GaussianIntegrals::ddexprLaplacian);
            } // end ifelse
        } // end fordd
        sums += tmpProdsd;
    } // end ford

    return sums * normalizationFactor(i) * normalizationFactor(j);
} // end function laplacianElement 

inline double GaussianIntegrals::potentialElement(const unsigned int& i, const
        unsigned int& j) {
    /* calculate and return the HO potential integral element <i|0.5wr^2|j> */
    double sums = 0.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        double tmpProdsd = 1.0;
        for (unsigned int dd = 0; dd < m_dim; ++dd) {
            const int& ndd = GaussianBasis::Cartesian::getn(i,dd);
            const int& mdd = GaussianBasis::Cartesian::getn(j,dd);
            if (dd != d) {
                tmpProdsd *= ddexpr(ndd, mdd,
                        &GaussianIntegrals::ddexprOverlap);
            } else {
                tmpProdsd *= ddexpr(ndd, mdd,
                        &GaussianIntegrals::ddexprPotential);
            } // end ifelse
        } // end fordd
        sums += tmpProdsd;
    } // end ford

    return xScaleHalf*sums * normalizationFactor(i) * normalizationFactor(j);
} // end function potentialElement

inline double GaussianIntegrals::coulombElement2D(const unsigned int& ix, const
        unsigned int& iy, const unsigned int& jx, const unsigned int& jy, const
        unsigned int& kx, const unsigned int& ky, const unsigned int& lx, const
        unsigned int& ly) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> in 2D for level (i,k,j,l) (calculate 1D integral
     * numerically) */
    double sum = 0.0;
    for (unsigned int px = 0; px <= ix+kx; ++px) {
        const double& cpx = coeffs->coeff(ix,kx,px);
        for (unsigned int py = 0; py <= iy+ky; ++py) {
            int pSign = (((px+py)%2==0) ? 1 : -1);
            const double& cpy = coeffs->coeff(iy,ky,py);
            for (unsigned int qx = 0; qx <= jx+lx; ++qx) {
                const double& cqx = coeffs->coeff(jx,lx,qx);
                for (unsigned int qy = 0; qy <= jy+ly; ++qy) {
                    sum += cpx * cpy * cqx * coeffs->coeff(jy,ly,qy) *
                        coeffs->auxiliary2D(0, px+qx, py+qy) * pSign;
                } // end forqy
            } // end forqx
        } // end forpy
    } // end forpx
    return sum;
} // end function coulombElement2D

inline double GaussianIntegrals::coulomb2D(const unsigned int& i, const unsigned
        int& j, const unsigned int& k, const unsigned int& l) {
    /* calculate full coulomb integral in 2D case */

    // grab hermite coefficients
    const int& Ix = GaussianBasis::Cartesian::getn(i,0);
    const int& Iy = GaussianBasis::Cartesian::getn(i,1);
    const int& Jx = GaussianBasis::Cartesian::getn(j,0);
    const int& Jy = GaussianBasis::Cartesian::getn(j,1);
    const int& Kx = GaussianBasis::Cartesian::getn(k,0);
    const int& Ky = GaussianBasis::Cartesian::getn(k,1);
    const int& Lx = GaussianBasis::Cartesian::getn(l,0);
    const int& Ly = GaussianBasis::Cartesian::getn(l,1);

    double sum = 0.0;
    for (int ix = ((Ix+1)%2==0 ? 1 : 0); ix <= Ix; ix+=2)
    for (int iy = ((Iy+1)%2==0 ? 1 : 0); iy <= Iy; iy+=2)
    for (int jx = ((Jx+1)%2==0 ? 1 : 0); jx <= Jx; jx+=2)
    for (int jy = ((Jy+1)%2==0 ? 1 : 0); jy <= Jy; jy+=2)
    for (int kx = ((Kx+1)%2==0 ? 1 : 0); kx <= Kx; kx+=2)
    for (int ky = ((Ky+1)%2==0 ? 1 : 0); ky <= Ky; ky+=2)
    for (int lx = ((Lx+1)%2==0 ? 1 : 0); lx <= Lx; lx+=2)
    for (int ly = ((Ly+1)%2==0 ? 1 : 0); ly <= Ly; ly+=2)
    {
        sum += xScaleSqrtPow(ix+iy + kx+ky + jx+jy + lx+ly) *
            HC::coeff(Ix,ix) * HC::coeff(Iy,iy) * HC::coeff(Kx,kx) *
            HC::coeff(Ky,ky) * HC::coeff(Jx,jx) * HC::coeff(Jy,jy) *
            HC::coeff(Lx,lx) * HC::coeff(Ly,ly) * coulombElement2D(ix,iy,
                    jx,jy, kx,ky, lx,ly);
    } // end for ix,iy,jx,jy,kx,ky,lx,ly

    return sum * coulomb2DFactor;
} // end function coulomb2D

inline double GaussianIntegrals::coulombElement3D(const unsigned int& ix, const
        unsigned int& iy, const unsigned int& iz, const unsigned int& jx, const
        unsigned int& jy, const unsigned int& jz, const unsigned int& kx, const
        unsigned int& ky, const unsigned int& kz, const unsigned int& lx, const
        unsigned int& ly, const unsigned int& lz) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> in 3D for level (i,k,j,l) (calculate 1D integral
     * numerically) */
    double sum = 0.0;
    for (unsigned int px = 0; px <= ix+kx; ++px) {
        const double& cpx = coeffs->coeff(ix,kx,px);
        for (unsigned int py = 0; py <= iy+ky; ++py) {
            const double& cpy =coeffs->coeff(iy,ky,py);
            for (unsigned int pz = 0; pz <= iz+kz; ++pz) {
                const double& cpz = coeffs->coeff(iz,kz,pz);
                int pSign = (((px+py+pz)%2==0) ? 1 : -1);
                for (unsigned int qx = 0; qx <= jx+lx; ++qx) {
                    const double& cqx = coeffs->coeff(jx,lx,qx);
                    for (unsigned int qy = 0; qy <= jy+ly; ++qy) {
                        const double& cqy = coeffs->coeff(jy,ly,qy);
                        for (unsigned int qz = 0; qz <= jz+lz; ++qz) {
                            sum += cpx * cpy * cpz * cqx * cqy *
                                coeffs->coeff(jz,lz,qz) *
                                coeffs->auxiliary3D(0, px+qx, py+qy, pz+qz) *
                                pSign;
                        } // end forqz
                    } // end forqy
                } // end forqx
            } // end forpz
        } // end forpy
    } // end forpx
    return sum;
} // end function coulombElement3D

inline double GaussianIntegrals::coulomb3D(const unsigned int& i, const unsigned
        int& j, const unsigned int& k, const unsigned int& l) {
    /* calculate full coulomb integral in 2D case */
    const int& Ix = GaussianBasis::Cartesian::getn(i,0);
    const int& Iy = GaussianBasis::Cartesian::getn(i,1);
    const int& Iz = GaussianBasis::Cartesian::getn(i,2);
    const int& Jx = GaussianBasis::Cartesian::getn(j,0);
    const int& Jy = GaussianBasis::Cartesian::getn(j,1);
    const int& Jz = GaussianBasis::Cartesian::getn(j,2);
    const int& Kx = GaussianBasis::Cartesian::getn(k,0);
    const int& Ky = GaussianBasis::Cartesian::getn(k,1);
    const int& Kz = GaussianBasis::Cartesian::getn(k,2);
    const int& Lx = GaussianBasis::Cartesian::getn(l,0);
    const int& Ly = GaussianBasis::Cartesian::getn(l,1);
    const int& Lz = GaussianBasis::Cartesian::getn(l,2);
    
    double sum = 0.0;
    for (int ix = ((Ix+1)%2==0 ? 1 : 0); ix <= Ix; ix+=2)
    for (int iy = ((Iy+1)%2==0 ? 1 : 0); iy <= Iy; iy+=2)
    for (int iz = ((Iz+1)%2==0 ? 1 : 0); iz <= Iz; iz+=2)
    for (int jx = ((Jx+1)%2==0 ? 1 : 0); jx <= Jx; jx+=2)
    for (int jy = ((Jy+1)%2==0 ? 1 : 0); jy <= Jy; jy+=2)
    for (int jz = ((Jz+1)%2==0 ? 1 : 0); jz <= Jz; jz+=2)
    for (int kx = ((Kx+1)%2==0 ? 1 : 0); kx <= Kx; kx+=2)
    for (int ky = ((Ky+1)%2==0 ? 1 : 0); ky <= Ky; ky+=2)
    for (int kz = ((Kz+1)%2==0 ? 1 : 0); kz <= Kz; kz+=2)
    for (int lx = ((Lx+1)%2==0 ? 1 : 0); lx <= Lx; lx+=2)
    for (int ly = ((Ly+1)%2==0 ? 1 : 0); ly <= Ly; ly+=2)
    for (int lz = ((Lz+1)%2==0 ? 1 : 0); lz <= Lz; lz+=2)
    {
        sum += xScaleSqrtPow(ix+iy+iz + kx+ky+kz + jx+jy+jz + lx+ly+lz) *
            HC::coeff(Ix,ix) * HC::coeff(Iy,iy) * HC::coeff(Iz,iz) *
            HC::coeff(Kx,kx) * HC::coeff(Ky,ky) * HC::coeff(Kz,kz) *
            HC::coeff(Jx,jx) * HC::coeff(Jy,jy) * HC::coeff(Jz,jz) *
            HC::coeff(Lx,lx) * HC::coeff(Ly,ly) * HC::coeff(Lz,lz) *
            coulombElement3D(ix,iy,iz, jx,jy,jz, kx,ky,kz, lx,ly,lz);
    } // end for ix,iy,iz,jx,jy,jz,kx,ky,kz,lx,ly,lz

    return sum * coulomb3DFactor;
} // end function coulomb3D

double GaussianIntegrals::overlapElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate and return the overlap integral element <i|j> */
    double prod = 1.0;
    for (unsigned int d = 0; d < m_dim; ++d) {
        prod *= ddexpr(GaussianBasis::Cartesian::getn(i,d),
                GaussianBasis::Cartesian::getn(j,d),
                &GaussianIntegrals::ddexprOverlap);
    } // end ford

    return prod * normalizationFactor(i) * normalizationFactor(j);
} // end function overlapElement

double GaussianIntegrals::oneBodyElement(const unsigned int& i, const unsigned
        int& j) {
    /* calculate and return oneBodyElement <i|h|k> = <i|K|j> + <i|P|j>, where K
     * is the kinetic part and P is the potential part */

    return -0.5*laplacianElement(i,j) + potentialElement(i,j);
} // end function oneBodyElements

double GaussianIntegrals::coulombElement(const unsigned int& i, const unsigned
    int& j, const unsigned int& k, const unsigned int& l) {
    /* calculate and return the two-body coulomb integral element
     * <ij|1/r_12|kl> */
    nSum = (GaussianBasis::Cartesian::getnStates(i) +
            GaussianBasis::Cartesian::getnStates(j) +
            GaussianBasis::Cartesian::getnStates(k) +
            GaussianBasis::Cartesian::getnStates(l));

    bool integrandIsEven = false;
    for (unsigned int i = 0; i < nSum.size(); ++i) {
        for (unsigned int j = 0; j < nSum.size(); ++j) {
            if (nSum(i)%2 == nSum(j)%2) {
                integrandIsEven = true;
            } else { 
                integrandIsEven = false;
                break;
            } // end if
        } // end forj
    } // end fori

    if (integrandIsEven) {
        /* make sure integrand is even (odd integrand yields zero) */
        return normalizationFactors(i) * normalizationFactors(j) *
            normalizationFactors(k) * normalizationFactors(l) *
            (this->*coulombFunc)(i,j,k,l);
    } else {
        /* return 0 if integrand is odd (in which case integral is zero) */
        return 0.0;
    } // end ifselse
} // end function coulombElement
