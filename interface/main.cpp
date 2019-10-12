#include <iomanip>
#include <cstdlib>
#include <vector>

#include <cartint.hpp>

int main(int argc, char *argv[]) {
    unsigned int dim = std::atoi(argv[1]);
    unsigned int cut_off = std::atoi(argv[2]);
    double scaling = std::atof(argv[3]);
    double omega = std::atof(argv[4]);

    GaussianIntegrals* I = new GaussianIntegrals(dim, cut_off, scaling);
    I->initializeParameters(omega);

    double num_basis = cut_off / 2;
    for (unsigned int p = 0; p < num_basis; ++p) {
        for (unsigned int q = 0; q < num_basis; ++q) {
                    // std::cout << std::setprecision(15) << p << " " << q << " " << I->oneBodyElement(p, q) << std::endl;
            for (unsigned int r = 0; r < num_basis; ++r) {
                for (unsigned int s = 0; s < num_basis; ++s) {
                    std::cout << std::setprecision(15) << p << " " << q << " " << r << " " << s << " " << I->coulombElement(p, q, r, s) << std::endl;
                } // end fors
            } // end forr
        } // end forq
    } // end forp


    delete I;
    return 0;
}
