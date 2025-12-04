#include <iostream>
#include <vector>

#include "../src/autodiff.hpp"

using namespace ASC_ode;

// ---------------------------------------------------------
// Build Legendre polynomials P0..P5 using recursion
// ---------------------------------------------------------
template <size_t N, typename T = double>
std::vector<AutoDiff<N,T>> Legendre(const AutoDiff<N,T>& x)
{
    std::vector<AutoDiff<N,T>> P(6);

    P[0] = AutoDiff<N,T>(1.0);
    P[1] = x;

    for (int n = 2; n <= 5; n++)
    {
        double a = (2.0*n - 1.0)/n;
        double b = (n - 1.0)/n;
        P[n] = a * x * P[n-1] - b * P[n-2];
    }

    return P;
}

// ---------------------------------------------------------
// MAIN: prints values and derivatives for P0..P5 on [-1,1]
// ---------------------------------------------------------
int main()
{
    const int N = 1; // only x is the variable

    std::cout << "x   P0  dP0  P1  dP1  P2  dP2  P3  dP3  P4  dP4  P5  dP5\n";

    for (int i = 0; i <= 40; i++)
    {
        double xv = -1.0 + 2.0 * i / 40.0;

        AutoDiff<N> x = Variable<0,double>(xv);

        auto P = Legendre<N>(x);

        std::cout << xv << " ";
        for (int n = 0; n <= 5; n++)
        {
            std::cout << P[n].value() << " " << P[n].deriv()[0] << "  ";
        }
        std::cout << "\n";
    }
    return 0;
}
