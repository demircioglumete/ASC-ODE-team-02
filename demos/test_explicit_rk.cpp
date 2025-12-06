#include <iostream>
#include <fstream>
#include <memory>

#include "nonlinfunc.hpp"
#include "explicitRK.hpp"

using namespace ASC_ode;

// Simple mass-spring system: y = [x, v], y' = [v, -k/m * x]
class MassSpring : public NonlinearFunction
{
private:
    double m_mass;
    double m_stiffness;

public:
    MassSpring(double m, double k) : m_mass(m), m_stiffness(k) {}

    size_t dimX() const override { return 2; }
    size_t dimF() const override { return 2; }

    void evaluate(VectorView<double> x, VectorView<double> f) const override
    {
        f(0) = x(1);
        f(1) = -m_stiffness / m_mass * x(0);
    }

    void evaluateDeriv(VectorView<double>, MatrixView<double>) const override
    {
        // Explicit RK does not need Jacobian
    }
};

int main()
{
    // Define the RHS object for ODE y' = f(y)
    auto rhs = std::make_shared<MassSpring>(1.0, 1.0);

    // RK4 Butcher tableau (classical 4-stage)
    Matrix<> A(4,4);
    A = 0.0;
    A(1,0) = 0.5;
    A(2,1) = 0.5;
    A(3,2) = 1.0;

    Vector<> b = {1.0/6, 1.0/3, 1.0/3, 1.0/6};
    Vector<> c = {0.0, 0.5, 0.5, 1.0};

    // Construct the explicit RK stepper
    ExplicitRungeKutta stepper(rhs, A, b, c);

    // Simulation parameters
    double tend = 10.0;
    int steps = 500;
    double tau = tend / steps;

    Vector<> y = {1.0, 0.0}; // initial condition: x=1, v=0

    std::ofstream out("rk4_output.txt");
    out << 0.0 << " " << y(0) << " " << y(1) << "\n";

    // Time stepping loop
    for(int i = 1; i <= steps; i++)
    {
        stepper.doStep(tau, y);
        double t = i * tau;
        out << t << " " << y(0) << " " << y(1) << "\n";
    }

    std::cout << "Explicit RK4 simulation completed. Output written to rk4_output.txt\n";
    return 0;
}
