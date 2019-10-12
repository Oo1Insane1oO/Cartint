#ifndef METHODS_H
#define METHODS_H

#include <Eigen/Dense>
#include <string>
#include <iostream>
#include <iomanip>

class Methods {
    private:
        Methods ();
        virtual ~Methods ();

    public:
        static void printProgressBar(std::string&, const float&, const int,
                const std::string& extra="");
        static std::string stringPos(const unsigned int&, const int&);
        static int divider(const unsigned int&, const unsigned int&, const
                int&);
        static void setDisplacement(int[], int[], const int[], int);
        
        template <typename T> static inline T refSum(const Eigen::Matrix<T*,
                Eigen::Dynamic, 1>& vec) {
            /* sum values of a vector of pointers(de-reference and sum) */
            T sum = 0;
            for (unsigned int i = 0; i < vec.size(); ++i) {
                sum += *(vec(i));
            } // end fori
            return sum;
        } // end function refSum

        template<typename T> static void setMax(T& a, T b) {
            /* set a to max of a and b */
            if (b > a) {
                a = b;
            } // end if
        } // end function setMax
        
        template<typename T> static void setMin(T& a, T b) {
            /* set a to min of a and b */
            if (b < a) {
                a = b;
            } // end if
        } // end function setMax
        
        struct expand {
            /* struct for expanding function pack */
            template<typename... T> expand(T&&...) {}
        };

        template<typename T, typename... Args> static inline T max(T arg0,
                Args...  args) {
            /* return maximum of args */
            T prev = arg0;

            expand{0, (setMax(prev, args), 0)...};
            return prev;
        } // end function max
        
        template<typename T, typename... Args> static inline T min(T arg0,
                Args...  args) {
            /* return minimum of args */
            T prev = arg0;

            expand{0, (setMin(prev, args), 0)...};
            return prev;
        } // end function max

        template<typename T> static inline void print1(T arg) {
            std::cout << arg << " ";
        } // end function print1

        template<typename... Args> static inline void sepPrint(Args... args) {
            /* print args with spaces */
            expand{0, (print1(args), 0)...};
            std::cout << "\n";
        } // end function sepPrint

        template<typename T> static inline void refSum(Eigen::Array<T,
                Eigen::Dynamic, 1>& buffer, const Eigen::Array<T*,
                Eigen::Dynamic, 1>& refArray) {
            /* override for expansion with more than 1 array */
            for (unsigned int i = 0; i < buffer.size(); ++i) {
                buffer(i) += *(refArray(i));
            } // end fori
        } // end function refSum

        template<typename T, typename ...Args> static inline void
            refSum(Eigen::Array<T, Eigen::Dynamic, 1>& buffer, const Args&...
                    arrays) {
            /* sum arrays Elementwise */
            expand{0, (refSum(buffer, arrays), 0)...};
        } // end function refSum
};

#endif /* METHODS_H */
