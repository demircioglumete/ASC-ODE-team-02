#include <iostream>
#include <memory>
#include <string>

#include "nonlinfunc.hpp"
#include "explicitRK.hpp"

using namespace ASC_ode;

// y = [x, v], y' = [v, -k/m x]
class MassSpring : public NonlinearFunction
{
    double m, k;

public:
    MassSpring(double m_, double k_) : m(m_), k(k_) {}

    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }

    void evaluate(VectorView<double> x, VectorView<double> f) const override
    {
        f(0) = x(1);
        f(1) = -k / m * x(0);
    }

    void evaluateDeriv(VectorView<double>, MatrixView<double>) const override {}
};

int main(int argc, char** argv)
{
    std::string method = "rk4";
    if (argc > 1)
        method = argv[1];

    auto rhs = std::make_shared<MassSpring>(1.0, 1.0);

    
    std::unique_ptr<Matrix<>> A;
    std::unique_ptr<Vector<>> b;
    std::unique_ptr<Vector<>> c;

    if (method == "rk2")
    {
        std::cout << "Using RK2 (Midpoint)\n";

        A = std::make_unique<Matrix<>>(2, 2);
        *A = 0.0;
        (*A)(1, 0) = 0.5;

        b = std::make_unique<Vector<>>(Vector<>({0.0, 1.0}));
        c = std::make_unique<Vector<>>(Vector<>({0.0, 0.5}));
    }
    else
    {
        std::cout << "Using RK4\n";

        A = std::make_unique<Matrix<>>(4, 4);
        *A = 0.0;
        (*A)(1, 0) = 0.5;
        (*A)(2, 1) = 0.5;
        (*A)(3, 2) = 1.0;

        b = std::make_unique<Vector<>>(Vector<>({
            1.0 / 6.0,
            1.0 / 3.0,
            1.0 / 3.0,
            1.0 / 6.0
        }));

        c = std::make_unique<Vector<>>(Vector<>({
            0.0,
            0.5,
            0.5,
            1.0
        }));
    }

    ExplicitRungeKutta stepper(rhs, *A, *b, *c);

    Vector<> y({1.0, 0.0});   // x = 1, v = 0
double t  = 0.0;
double dt = 0.01;

for (int n = 0; n < 10; ++n)
{
    stepper.doStep(dt, y);
    t += dt;

    std::cout << "t = " << t
              << "  x = " << y(0)
              << "  v = " << y(1) << std::endl;
}


    return 0;
}
