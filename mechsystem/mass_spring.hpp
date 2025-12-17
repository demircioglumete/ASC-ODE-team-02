#ifndef MASS_SPRING_HPP
#define MASS_SPRING_HPP

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;

#include <vector.hpp>
using namespace nanoblas;


template <int D>
class Mass
{
public:
  double mass;
  Vec<D> pos;
  Vec<D> vel = 0.0;
  Vec<D> acc = 0.0;
};


template <int D>
class Fix
{
public:
  Vec<D> pos;
};


class Connector
{
public:
  enum CONTYPE { FIX=1, MASS=2 };
  CONTYPE type;
  size_t nr;
};

std::ostream & operator<< (std::ostream & ost, const Connector & con)
{
  ost << "type = " << int(con.type) << ", nr = " << con.nr;
  return ost;
}

class Spring
{
public:
  double length;  
  double stiffness;
  std::array<Connector,2> connectors;
};

/// Distance constraint between two connectors: |p2 - p1| = length
struct DistanceConstraint
{
  Connector connectors[2];  // same idea as in Spring
  double    length;         // prescribed distance
};

template <int D>
class MassSpringSystem
{
  std::vector<Fix<D>> m_fixes;
  std::vector<Mass<D>> m_masses;
  std::vector<Spring> m_springs;
  std::vector<DistanceConstraint> m_constraints;
  Vec<D> m_gravity=0.0;
public:

// read-only access (comme pour masses(), fixes(), springs())
auto const & constraints() const { return m_constraints; }

void addDistanceConstraint(Connector c1, Connector c2, double length)
{
  DistanceConstraint dc;
  dc.connectors[0] = c1;
  dc.connectors[1] = c2;
  dc.length        = length;
  m_constraints.push_back(dc);
}

  void setGravity (Vec<D> gravity) { m_gravity = gravity; }
  Vec<D> getGravity() const { return m_gravity; }

  Connector addFix (Fix<D> p)
  {
    m_fixes.push_back(p);
    return { Connector::FIX, m_fixes.size()-1 };
  }

  Connector addMass (Mass<D> m)
  {
    m_masses.push_back (m);
    return { Connector::MASS, m_masses.size()-1 };
  }
  
  size_t addSpring (Spring s) 
  {
    m_springs.push_back (s); 
    return m_springs.size()-1;
  }

  auto & fixes() { return m_fixes; } 
  auto & masses() { return m_masses; } 
  auto & springs() { return m_springs; }

  void getState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        valmat.row(i) = m_masses[i].pos;
        dvalmat.row(i) = m_masses[i].vel;
        ddvalmat.row(i) = m_masses[i].acc;
      }
  }

  void setState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        m_masses[i].pos = valmat.row(i);
        m_masses[i].vel = dvalmat.row(i);
        m_masses[i].acc = ddvalmat.row(i);
      }
  }
};

template <int D>
std::ostream & operator<< (std::ostream & ost, MassSpringSystem<D> & mss)
{
  ost << "fixes:" << std::endl;
  for (auto f : mss.fixes())
    ost << f.pos << std::endl;

  ost << "masses: " << std::endl;
  for (auto m : mss.masses())
    ost << "m = " << m.mass << ", pos = " << m.pos << std::endl;

  ost << "springs: " << std::endl;
  for (auto sp : mss.springs())
    ost << "length = " << sp.length << ", stiffness = " << sp.stiffness
        << ", C1 = " << sp.connectors[0] << ", C2 = " << sp.connectors[1] << std::endl;
  return ost;
}


template <int D>
class MSS_Function : public NonlinearFunction
{
  MassSpringSystem<D> & mss;

public:
  MSS_Function (MassSpringSystem<D> & _mss)
    : mss(_mss) { }

  virtual size_t dimX() const override { return D * mss.masses().size(); }
  virtual size_t dimF() const override { return D * mss.masses().size(); }

 
  virtual void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f = 0.0;

    const size_t nm = mss.masses().size();
    auto xmat = x.asMatrix(nm, D); 
    auto fmat = f.asMatrix(nm, D);  

    for (size_t i = 0; i < nm; i++)
      fmat.row(i) = mss.masses()[i].mass * mss.getGravity();

    for (const auto &spring : mss.springs())
    {
      auto c1 = spring.connectors[0];
      auto c2 = spring.connectors[1];

      Vec<D> p1, p2;

      if (c1.type == Connector::FIX)
        p1 = mss.fixes()[c1.nr].pos;
      else
        p1 = xmat.row(c1.nr);

      if (c2.type == Connector::FIX)
        p2 = mss.fixes()[c2.nr].pos;
      else
        p2 = xmat.row(c2.nr);

      Vec<D> diff = p2 - p1;
      const double dist = norm(diff);
      if (dist <= 1e-12) continue;

      const double force = spring.stiffness * (dist - spring.length);
      Vec<D> dir12 = (1.0 / dist) * diff;   

      if (c1.type == Connector::MASS)
        fmat.row(c1.nr) += force * dir12;   

      if (c2.type == Connector::MASS)
        fmat.row(c2.nr) -= force * dir12;   
    }

    for (size_t i = 0; i < nm; i++)
      fmat.row(i) *= 1.0 / mss.masses()[i].mass;

   

    for (const auto &con : mss.constraints())
    {
      const auto c1 = con.connectors[0];
      const auto c2 = con.connectors[1];

      const bool isMass1 = (c1.type == Connector::MASS);
      const bool isMass2 = (c2.type == Connector::MASS);

      // positions
      Vec<D> p1, p2;
      if (c1.type == Connector::FIX)
        p1 = mss.fixes()[c1.nr].pos;
      else
        p1 = xmat.row(c1.nr);

      if (c2.type == Connector::FIX)
        p2 = mss.fixes()[c2.nr].pos;
      else
        p2 = xmat.row(c2.nr);

      Vec<D> diff = p2 - p1;
      const double dist = norm(diff);
      if (dist <= 1e-12) continue; 

      Vec<D> dir = (1.0 / dist) * diff;

      double invm1 = 0.0;
      double invm2 = 0.0;
      if (isMass1) invm1 = 1.0 / mss.masses()[c1.nr].mass;
      if (isMass2) invm2 = 1.0 / mss.masses()[c2.nr].mass;

      const double denom = dist * (invm1 + invm2);
      if (denom <= 1e-12) continue;

      Vec<D> a1 = 0.0, a2 = 0.0;
      if (isMass1) a1 = fmat.row(c1.nr);
      if (isMass2) a2 = fmat.row(c2.nr);

      double num = 0.0;
      for (int d = 0; d < D; ++d)
        num += diff(d) * (a2(d) - a1(d));

      const double lambda = - num / denom;

      if (isMass1)
        fmat.row(c1.nr) += (-lambda * invm1) * dir;

      if (isMass2)
        fmat.row(c2.nr) += (+lambda * invm2) * dir;
    }
  }
  
  virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    // TODO: exact differentiation
    double eps = 1e-8;
    Vector<> xl(dimX()), xr(dimX()), fl(dimF()), fr(dimF());
    for (size_t i = 0; i < dimX(); i++)
      {
        xl = x;
        xl(i) -= eps;
        xr = x;
        xr(i) += eps;
        evaluate (xl, fl);
        evaluate (xr, fr);
        df.col(i) = 1/(2*eps) * (fr-fl);
      }
  }
  
};

#endif
