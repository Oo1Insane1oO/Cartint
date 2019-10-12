#ifndef CARTESIAN_H
#define CARTESIAN_H

#include <string>

#include <Eigen/Dense>

#include <cartint/methods.h>

class Cartesian {
    private:
        int s;
        unsigned int m_dim, m_size, m_numStates;
        Eigen::VectorXi n, ms, E, M, m;

        Eigen::VectorXi angularMomenta;

        Eigen::MatrixXi states;

        void addState(unsigned int&, const unsigned int, const unsigned int,
                const unsigned int, const unsigned int, const unsigned int);
        void addState(unsigned int&, const unsigned int, const unsigned int,
                const unsigned int, const unsigned int);
        void addState(unsigned int&, const unsigned int, const unsigned int,
                const unsigned int);
        void findPrincipal(const unsigned int& , Eigen::MatrixXi&);
        void sumn();

        unsigned int calculateDegeneracy(const unsigned int&);

        void setStates(const unsigned int&);
    
    public:
        Cartesian ();
        virtual ~Cartesian ();

        void setup(unsigned int, const unsigned int);

        const Eigen::MatrixXi& getStates() const;
        const Eigen::Ref<const Eigen::VectorXi> getStates(const unsigned int&)
            const;
        const Eigen::Ref<const Eigen::VectorXi> getnStates(const unsigned int&)
            const;
        const int& getn(const unsigned int&, const unsigned int&) const;
        const Eigen::VectorXi &getn() const;
        const int &getn(const unsigned int&) const;
        const Eigen::VectorXi &getE() const;
        const int &getE(const unsigned int&) const;
        const Eigen::VectorXi &getMagic() const;
        const int &getMagic(const unsigned int&) const;
        const int &getSumn(const unsigned int&) const;
        const Eigen::VectorXi& getSumn() const;
        const unsigned int& getSize() const;
        const unsigned int& getNumberOfStates() const;

        void restructureStates();

        void printStates();

        void writeStatesToFile(std::string);
};

#endif /* CARTESIAN_H */
