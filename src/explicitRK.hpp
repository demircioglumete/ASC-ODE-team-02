#pragma once

#include "timestepper.hpp"

namespace ASC_ode {

class ExplicitRungeKutta : public TimeStepper
{
private:
    Matrix<> A;
    Vector<> b;
    Vector<> c;
    int s;

public:
    ExplicitRungeKutta(std::shared_ptr<NonlinearFunction> rhs,
                       const Matrix<>& A_,
                       const Vector<>& b_,
                       const Vector<>& c_)
        : TimeStepper(rhs), A(A_), b(b_), c(c_), s(c_.size())
    {}

    void doStep(double tau, VectorView<double> y) override
    {
        int n = y.size();

        std::vector<Vector<>> k(s, Vector<>(n));
        Vector<> ytemp(n);

        for (int j = 0; j < s; j++)
        {
            ytemp = y;

            for (int i = 0; i < j; i++)
            {
                ytemp += tau * A(j,i) * k[i];
            }

            m_rhs->evaluate(ytemp, k[j]);   
        }

        for (int j = 0; j < s; j++)
            y += tau * b(j) * k[j];
    }
};

} // namespace ASC_ode
