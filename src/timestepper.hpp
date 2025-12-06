#ifndef TIMERSTEPPER_HPP
#define TIMERSTEPPER_HPP

#include <functional>
#include <exception>

#include "Newton.hpp"


namespace ASC_ode
{
  
  class TimeStepper
  { 
  protected:
    std::shared_ptr<NonlinearFunction> m_rhs;
  public:
    TimeStepper(std::shared_ptr<NonlinearFunction> rhs) : m_rhs(rhs) {}
    virtual ~TimeStepper() = default;
    virtual void DoStep(double tau, VectorView<double> y) = 0;
  };

  class ExplicitEuler : public TimeStepper
  {
    Vector<> m_vecf;
  public:
    ExplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_vecf(rhs->dimF()) {}
    void DoStep(double tau, VectorView<double> y) override
    {
      this->m_rhs->evaluate(y, m_vecf);
      y += tau * m_vecf;
    }
  };

  class ImplicitEuler : public TimeStepper
  {
    std::shared_ptr<NonlinearFunction> m_equ;
    std::shared_ptr<Parameter> m_tau;
    std::shared_ptr<ConstantFunction> m_yold;
  public:
    ImplicitEuler(std::shared_ptr<NonlinearFunction> rhs) 
    : TimeStepper(rhs), m_tau(std::make_shared<Parameter>(0.0)) 
    {
      m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
      auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());
      m_equ = ynew - m_yold - m_tau * m_rhs;
    }

    void DoStep(double tau, VectorView<double> y) override
    {
      m_yold->set(y);
      m_tau->set(tau);
      NewtonSolver(m_equ, y);
    }
  };

class ImprovedEuler : public TimeStepper
  {
    Vector<> m_vecf;      // will hold f(y_n) or f(y_tilde)
    Vector<> m_ytilde;    // temporary y_tilde

  public:
    ImprovedEuler(std::shared_ptr<NonlinearFunction> rhs)
      : TimeStepper(rhs),
        m_vecf(rhs->dimF()),
        m_ytilde(rhs->dimX())
    { }

    void DoStep(double tau, VectorView<> y) override
    {
      // 1) Compute f(y_n)
      this->m_rhs->evaluate(y, m_vecf);

      // 2) y_tilde = y_n + (tau/2) * f(y_n)
      m_ytilde = y;
      m_ytilde += 0.5 * tau * m_vecf;

      // 3) Compute f(y_tilde)
      this->m_rhs->evaluate(m_ytilde, m_vecf);

      // 4) y_{n+1} = y_n + tau * f(y_tilde)
      y += tau * m_vecf;
    }
  };


// Crank Nicolson 

class CrankNicolson : public TimeStepper
{
    // G(y_{n+1}) = y_{n+1} - y_n - τ * ( 0.5 f(y_n) + 0.5 f(y_{n+1}) )
    // Newton ile G(y_{n+1}) = 0 çözüyoruz.

    std::shared_ptr<NonlinearFunction> m_equ;   // G
    std::shared_ptr<Parameter>         m_tau;   // τ
    std::shared_ptr<ConstantFunction>  m_yold;  // y_n
    std::shared_ptr<ConstantFunction>  m_fold;  // f(y_n)

public:
    CrankNicolson(std::shared_ptr<NonlinearFunction> rhs)
      : TimeStepper(rhs),
        m_tau(std::make_shared<Parameter>(0.0))
    {
        // y_n (sabit vektör)
        m_yold = std::make_shared<ConstantFunction>(rhs->dimX());
        // f(y_n) (sabit vektör)
        m_fold = std::make_shared<ConstantFunction>(rhs->dimF());

        // y_{n+1} kimliği
        auto ynew = std::make_shared<IdentityFunction>(rhs->dimX());

        // 0.5 f(y_{n+1}) + 0.5 f(y_n)
        std::shared_ptr<NonlinearFunction> f_comb =
            0.5 * m_rhs
          + 0.5 * std::static_pointer_cast<NonlinearFunction>(m_fold);

        // G(y_{n+1}) = y_{n+1} - y_n - τ * (0.5 f(y_n) + 0.5 f(y_{n+1}))
        m_equ = ynew - m_yold - m_tau * f_comb;
    }

    void DoStep(double tau, VectorView<> y) override
    {
        // y_n kaydet
        m_yold->set(y);

        // f(y_n) hesapla ve sabit fonksiyona koy
        Vector<> f_n(m_rhs->dimF());
        m_rhs->evaluate(y, f_n);
        m_fold->set(f_n);

        // τ güncelle
        m_tau->set(tau);

        // Newton ile G(y_{n+1}) = 0 çöz → sonuç doğrudan y'ye yazılıyor
        NewtonSolver(m_equ, y);
    }
};

}


#endif
