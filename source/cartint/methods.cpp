#include <cartint/methods.h>
#include <cmath>

void Methods::printProgressBar(std::string &flip, const float 
        &percentage, const int width, const std::string &extra) {
    /* print progressbar with percentage as output */
    std::cout << flip + std::string(width*percentage,'-') +
        std::string(width*(1-percentage),' ') + "] " +
        std::to_string((int)(round(percentage*100.0))) + " % " + extra + "\n";
    std::cout.flush();
} // end function printProgressBar

std::string Methods::stringPos(const unsigned int &rank, const int &displ) {
    /* return a position for progress based on rank and displacement*/
    return "\033[" + std::to_string(rank + displ) + "H\r";
} // end function progressPositions

int Methods::divider(const unsigned int &index, const unsigned int &maxValue,
        const int &divisor) {
    /* divide a range in even whole number sections, return a given value
     * depending on index */
    return static_cast<int>(index < (maxValue%divisor) ?
            ceil((double)maxValue/divisor) : floor((double)maxValue/divisor));
} // end function divider

void Methods::setDisplacement(int recvc[], int displ[], const int sizes[], int
        P) {
    /* function for setting receive count and displacement used in MPI gatherv
     * and scatterv functions */
    for (int p = 0; p < P; ++p) {
        recvc[p] = sizes[p];
        displ[p] = sizes[p] * p;
    } // end forp
} // end function setDisplacement
