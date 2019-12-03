#include<pybind11/pybind11.h>

#include <cartint/integral.hpp>

int main(void) {
    return 0;
}

PYBIND11_MODULE(intobj, m) {
    cartint::Integral2D intObj = cartint::Integral2D();
    m.def("overlap", [&intObj](){return intObj.overlap_element(0,0);});
}
