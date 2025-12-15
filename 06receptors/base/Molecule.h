#include <list>
#include <stdio.h>
#include "Vect.h"

class Molecule
{
public:
  Molecule();
  Molecule(Vect3D pos_, double time);
  Molecule(Vect3D pos_, double time, int type_);

  Vect3D v;
  Vect3D f;
  double r;
  Vect3D pos;
  double age;
  double time;
  bool MustDelete;
  double Diff;

  int type; // type = 0 for molecules and type = 1 for treatment
  void Move(double dt, gsl_rng *rng);
  void Step(double dt);
  double NormFunction(double x, double delta, double mu);
};

void Molecule::Move(double dt, gsl_rng *rng)
{
  double r1 = 0;
  double r2 = 0;
  double r3 = 0;

  if (type == 1) Diff = 5*diff_egf;
  r1 = gsl_ran_gaussian_ziggurat(rng, 1.0) * sqrt(2 * Diff * dt);
  r2 = gsl_ran_gaussian_ziggurat(rng, 1.0) * sqrt(2 * Diff * dt);
  r3 = gsl_ran_gaussian_ziggurat(rng, 1.0) * sqrt(2 * Diff * dt);

  double dref = sqrt(SQ(pos.x + r1) + SQ(pos.y + r2) + SQ(pos.z + r3));
  double d1 = r - sqrt(SQ(pos.x) + SQ(pos.y) + SQ(pos.z));

  if (sqrt(SQ(pos.x + r1) + SQ(pos.y + r2) + SQ(pos.z + r3)) < 50)
  {
    pos.x += r1;
    pos.y += r2;
    pos.z += r3;
  }

  else
  { // reflective BC
    pos.x += ((50 - d1) / (dref - d1)) * r1 - ((dref - r) / (dref - d1)) * r1;
    pos.y += ((50 - d1) / (dref - d1)) * r2 - ((dref - r) / (dref - d1)) * r2;
    pos.z += ((50 - d1) / (dref - d1)) * r3 - ((dref - r) / (dref - d1)) * r3;
  }
}

void Molecule::Step(double dt)
{
  age += dt;
  if (age >= time)
    MustDelete = true;
}

double Molecule::NormFunction(double x, double delta, double mu)
{
  double f = 1 / (delta * sqrt(2 * M_PI)) * exp(-(x - mu) * (x - mu) / (2 * delta * delta));
  return f;
}

Molecule::Molecule(Vect3D pos_, double time_)
{
  age = 0;
  pos = pos_;
  time = time_;
  Diff = diff_egf;
  r = 0.02;
  MustDelete = false;
  type = 0;
}

Molecule::Molecule(Vect3D pos_, double time_, int type_)
{
  age = 0;
  pos = pos_;
  time = time_;
  Diff = diff_egf;
  r = 0.02;
  MustDelete = false;
  type = type_;
}

Molecule::Molecule()
{
  age = 0;
  time = 5000;
  r = 0.01;
  type = 0;
}
