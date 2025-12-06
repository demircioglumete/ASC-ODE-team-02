#include <iostream>
#include <fstream>
#include <cmath>

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;


// -----------------------------------------------------------------------------
// RC Circuit in autonomous form
//
// State vector:
//   x1 = U_C(t)  (capacitor voltage)
//   x2 = t       (time variable, added to make the system autonomous)
//
// ODE system:
//   x1' = (1/(R*C)) * ( cos(omega * x2) - x1 )
//   x2' = 1
//
// where omega = 100*pi.
//
// Initial condition:
//   U_C(0) = 0
//   t(0)   = 0
// -----------------------------------------------------------------------------

class RCCircuit : public NonlinearFunction
{
private:
  double R;
  double C;
  double omega;   // driving frequency: 100*pi

public:
  RCCircuit(double R_, double C_)
    : R(R_), C(C_), omega(100.0 * M_PI) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }

  // Evaluate the ODE function f(x)
  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    double UC = x(0);   // capacitor voltage U_C
    double t  = x(1);   // time variable

    // x1' = dU_C/dt
    f(0) = (1.0/(R*C)) * ( std::cos(omega * t) - UC );

    // x2' = dt/dt = 1
    f(1) = 1.0;
  }

  // Jacobian (df/dx)
  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    double t = x(1);

    df = 0.0;

    // Partial derivatives of f0
    df(0,0) = -1.0/(R*C);                                 // ∂f0/∂x1
    df(0,1) = (1.0/(R*C)) * (-omega * std::sin(omega*t)); // ∂f0/∂x2

    // f1 = 1  → derivatives are zero
    df(1,0) = 0.0;
    df(1,1) = 0.0;
  }
};


int main()
{
  // Circuit parameters
 
  double R = 1;
  double C = 1;

   //double R = 100.0;
   //double C = 1e-6;

  // Time interval
  double t0   = 0.0;
  double tend = 0.02;        
  int    steps = 1000;
  double tau   = (tend - t0) / steps;

  // Initial condition: U_C(0) = 0,  t(0) = 0
  Vector<> y = { 0.0, t0 };

  // Create the right-hand side object
  auto rhs = std::make_shared<RCCircuit>(R, C);

  // Choose a time-stepping method:
  // ExplicitEuler  stepper(rhs);
   //ImplicitEuler  stepper(rhs);
   //ImprovedEuler  stepper(rhs);
 CrankNicolson stepper(rhs);

  // Output file: columns = t, U_C, U_0
  std::ofstream outfile("output_test_ode_circuit.txt");
  double omega = 100.0 * M_PI;

  double t  = y(1);
  double UC = y(0);
  double U0 = std::cos(omega * t);

  std::cout  << t << "  " << UC << "  " << U0 << std::endl;
  outfile    << t << "  " << UC << "  " << U0 << std::endl;

  for (int i = 0; i < steps; i++)
  {
    stepper.DoStep(tau, y);

    t  = y(1);
    UC = y(0);
    U0 = std::cos(omega * t);

    std::cout  << t << "  " << UC << "  " << U0 << std::endl;
    outfile    << t << "  " << UC << "  " << U0 << std::endl;
  }

  return 0;
}