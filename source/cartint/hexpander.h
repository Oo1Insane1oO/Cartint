#ifndef HEXPANDER_H
#define HEXPANDER_H

#include <Eigen/Dense>

class Hexpander {
    private:
        Eigen::ArrayXd integrals2D;
        Eigen::ArrayXd coefficients;
        Eigen::ArrayXd integrals3D;
    
        unsigned int iMp1, jMp1, nMp1;
        
        unsigned int x2DMp1, y2DMp1, n2DMp1;
        unsigned int x3DMp1, y3DMp1, z3DMp1, n3DMp1;

        double boys(const unsigned int&, const double&);

        double boysIntegrand(double, const unsigned int&, const double&);
        double modifiedIntegrand(double, const unsigned int&, const double&);

        bool checkIndices(const int&, const int&, const int&);

        inline unsigned int cidx(const unsigned int& i, const unsigned int& j,
                const unsigned int& n) const {
            /* calculate index in coefficients matrix */
            return n + nMp1 * (j + jMp1 * i);
        } // end function cidx
        
        inline unsigned int I2Didx(const unsigned int& i, const unsigned int&
                j, const unsigned int& n) const {
            /* calculate index in integrals2D matrix */
            return n + n2DMp1 * (j + y2DMp1 * i);
        } // end function I2Didx

        inline unsigned int I3Didx(const unsigned int& i, const unsigned int&
                j, const unsigned int& k, const unsigned int& n) const {
            /* calculate index in coefficients matrix */
            return n + n3DMp1 * (k + z3DMp1 * (j + y3DMp1 * i));
        } // end function cidx

    public:
        Hexpander();
        virtual ~Hexpander ();
        
        void setCoefficients(unsigned int, unsigned int, double, double, double);
        void setAuxiliary2D(unsigned int, unsigned int, double, double, double,
                double, const Eigen::VectorXd&);
        void setAuxiliary3D(unsigned int, unsigned int, unsigned int, double,
                double, double, double, const Eigen::VectorXd&);

        const double& coeff(const unsigned int&, const unsigned int&, const
                unsigned int&) const;
        const double& auxiliary2D(const unsigned int&, const unsigned int&,
                const unsigned int&) const;
        const double& auxiliary3D(const unsigned int&, const unsigned int&,
                const unsigned int&, const unsigned int&) const;
};

#endif /* HEXPANDER_H */
