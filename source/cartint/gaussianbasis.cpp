#include <cartint/gaussianbasis.h>

GaussianBasis::GaussianBasis() : Cartesian::Cartesian() {
    /* default constructor */
} // end constructor

GaussianBasis::GaussianBasis(unsigned int numSingleParticleStates, unsigned int
        dimension): Cartesian::Cartesian() {
    setup(numSingleParticleStates, dimension);
} // end constructor

GaussianBasis::~GaussianBasis() {
} // end deconstructor

unsigned int GaussianBasis::getSize() {
    /* return number of states (basis functions) */
    return Cartesian::getStates().rows();
} // end function getSize

void GaussianBasis::setup(unsigned int numSingleParticleStates, unsigned int
        dimension) {
    /* set number of dimensions */
    m_dim = dimension;
    Cartesian::setup(numSingleParticleStates, m_dim);
    Cartesian::restructureStates();
} // end function setup
