#ifndef GAUSSIANINTEGRALS_H
#define GAUSSIANINTEGRALS_H

#include <cartint/gaussianbasis.h>
#include <cartint/hexpander.h>

#include <string>

class GaussianBasis;

class GaussianIntegrals : protected GaussianBasis {
    private:
        double expScaleFactor, sqrtFactor, F0, F1, xScale, sqrtScale, powScale,
               sqrtScale1, xScaleHalf, coulomb2DFactor, coulomb3DFactor;

        unsigned int m_dim;

        Eigen::ArrayXd normalizationFactors;
        const double& normalizationFactor(const unsigned int&) const;

        Eigen::ArrayXd xScaleSqrtPow;

        Eigen::ArrayXi nSum;

        std::unique_ptr<Hexpander> coeffs;
        
        inline double laplacianElement(const unsigned int&, const unsigned
                int&);
        inline double potentialElement(const unsigned int&, const unsigned
                int&);

        inline double ddexprLaplacian(const int&, const int&);
        inline double ddexprPotential(const int&, const int&);

        inline double overlapd(const unsigned int&, const unsigned int&); 

        double (GaussianIntegrals::*coulombFunc)(const unsigned int&,
                const unsigned int&, const unsigned int&, const unsigned int&);
        inline double coulombElement2D(const unsigned int&, const unsigned
                int&, const unsigned int&, const unsigned int&, const unsigned
                int&, const unsigned int&, const unsigned int&, const unsigned
                int&);
        inline double coulombElement3D(const unsigned int&, const unsigned
                int&, const unsigned int&, const unsigned int&, const unsigned
                int&, const unsigned int&, const unsigned int&, const unsigned
                int&, const unsigned int&, const unsigned int&, const unsigned
                int&, const unsigned int&);
        inline double coulomb2D(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&);
        inline double coulomb3D(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&);

        void setFileName(const double&);

    protected:
        unsigned int m_cutOff;
        
        double ddexpr(const int&, const int&,
                double(GaussianIntegrals::*)(const int&, const int&));
        double ddexprOverlap(const int&, const int&);
        void setNormalizations();

    public:
        GaussianIntegrals(const unsigned int, unsigned int, double=1);
        virtual ~GaussianIntegrals();

        std::string initializeParameters(double);

        double overlapElement(const unsigned int&, const unsigned int&);
        double oneBodyElement(const unsigned int&, const unsigned int&);
        double coulombElement(const unsigned int&, const unsigned int&, const
                unsigned int&, const unsigned int&);
};

#endif /* GAUSSIANINTEGRALS_H */
