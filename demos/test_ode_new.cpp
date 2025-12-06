#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;


// Mass–spring system: x0 = position, x1 = velocity
class MassSpring : public NonlinearFunction
{
private:
  double mass;
  double stiffness;

public:
  MassSpring(double m, double k) : mass(m), stiffness(k) {}

  size_t dimX() const override { return 2; }
  size_t dimF() const override { return 2; }

  void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);                   // dx/dt = v
    f(1) = -stiffness/mass * x(0); // dv/dt = -k/m * x
  }

  void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1.0;
    df(1,0) = -stiffness/mass;
  }
};


int main(int argc, char** argv)
{
  // ------------------------------------------------------------------
  // 1) Read command line arguments
  // ------------------------------------------------------------------
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " METHOD STEPS\n";
    std::cerr << "  METHOD = exp | imp | impr | cn\n";
    std::cerr << "  STEPS  = number of time steps, e.g. 10 50 100\n";
    return 1;
  }

  std::string method = argv[1];        // "exp", "imp", "impr", "cn"
  int steps = std::stoi(argv[2]);      // e.g. 10, 50, 100

  double tend = 4 * M_PI;
  double tau  = tend / steps;

  // ------------------------------------------------------------------
  // 2) Create RHS and choose time stepper
  // ------------------------------------------------------------------
  auto rhs = std::make_shared<MassSpring>(1.0, 1.0);

  std::unique_ptr<TimeStepper> stepper;
  std::string prefix;

  if (method == "exp")
  {
    stepper = std::make_unique<ExplicitEuler>(rhs);
    prefix  = "exp";
  }
  else if (method == "imp")
  {
    stepper = std::make_unique<ImplicitEuler>(rhs);
    prefix  = "imp";
  }
  else if (method == "impr")
  {
    stepper = std::make_unique<ImprovedEuler>(rhs);
    prefix  = "impr";
  }
  else if (method == "cn")
  {
    stepper = std::make_unique<CrankNicolson>(rhs);
    prefix  = "cn";
  }
  else
  {
    std::cerr << "Unknown method: " << method << "\n";
    std::cerr << "Use: exp, imp, impr, cn\n";
    return 1;
  }

  // ------------------------------------------------------------------
  // 3) Build output file name: <prefix>_t<steps>.txt
  // ------------------------------------------------------------------
  std::ostringstream fname;
  fname << prefix << "_t" << steps << ".txt";
  std::string filename = fname.str();

  std::ofstream outfile(filename);
  if (!outfile)
  {
    std::cerr << "Cannot open output file: " << filename << "\n";
    return 1;
  }

  std::cout << "Running method = " << method
            << ", steps = " << steps
            << ", dt = " << tau
            << ", writing to " << filename << std::endl;

  // ------------------------------------------------------------------
  // 4) Initial condition and time stepping
  // ------------------------------------------------------------------
  Vector<> y = { 1.0, 0.0 }; // position = 1, velocity = 0

  double t = 0.0;
  outfile << t << " " << y(0) << " " << y(1) << std::endl;

  for (int i = 0; i < steps; i++)
  {
    stepper->DoStep(tau, y);
    t = (i+1) * tau;
    outfile << t << " " << y(0) << " " << y(1) << std::endl;
  }

  return 0;
}
