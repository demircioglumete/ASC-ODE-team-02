#include <iostream>
#include "../src/autodiff.hpp"

using namespace ASC_ode;

const double g = 9.81;
const double l = 1.0;

// -----------------------------------------------------
// Generic pendulum RHS: works for double AND AutoDiff
// -----------------------------------------------------
template <typename T>
std::array<T,2> T_evaluate(const std::array<T,2> &y)
{
    std::array<T,2> f;

    f[0] = y[1];
    f[1] = - (g/l) * sin(y[0]);

    return f;
}

int main()
{
    std::cout << "=== TEST PENDULUM RHS WITH DOUBLE ===\n";

    std::array<double,2> y_double = {1.0, 0.2};

    auto f_val = T_evaluate<double>(y_double);

    std::cout << "f0 = " << f_val[0] << "\n";
    std::cout << "f1 = " << f_val[1] << "\n";


    std::cout << "\n=== TEST PENDULUM RHS WITH AUTODIFF<J> ===\n";

    // 2 variables -> Jacobian is 2×2
    std::array<AutoDiff<2>,2> y_AD;
    y_AD[0] = Variable<0,double>(1.0);  // y0 = alpha
    y_AD[1] = Variable<1,double>(0.2);  // y1 = alpha'

    auto f_AD = T_evaluate<AutoDiff<2>>(y_AD);

    // Print values
    std::cout << "f0 value = " << f_AD[0].value() << "\n";
    std::cout << "f1 value = " << f_AD[1].value() << "\n";

    // Print Jacobian entries df_i/dy_j
    std::cout << "\nJacobian:\n";
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<2;j++)
            std::cout << f_AD[i].deriv()[j] << "  ";
        std::cout << "\n";
    }

    return 0;
}
