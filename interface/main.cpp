#include <cartint.hpp>
#include <cartint/methods.h>

int main(void) {
    unsigned int dim = 2;
    unsigned int num_basis = 6;
    double scaling = 1.0;
    double omega = 1.0;

    GaussianIntegrals* I = new GaussianIntegrals(dim, num_basis, scaling);
    I->initializeParameters(omega);

    for (unsigned int p = 0; p < num_basis; ++p) {
        for (unsigned int q = 0; q < num_basis; ++q) {
            for (unsigned int r = 0; r < num_basis; ++r) {
                for (unsigned int s = 0; s < num_basis; ++s) {
                    Methods::sepPrint(p, q, r, s, I->coulombElement(p, q, r, s));
                } // end fors
            } // end forr
        } // end forq
    } // end forp

    delete I;
    return 0;
}
