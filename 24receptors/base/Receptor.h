#include <list>
#include <stdio.h>
#include "Vect.h"

#define SQ(x) ((x) * (x))
enum status{inactive = 0, active = 1, blocked = 2};
using namespace std;

class Receptor
{ 
 public:
  Receptor();
  Receptor(Vect3D posEx_, Vect3D posIn_, double r0_, enum status cur_status_, double timer_);
 
  Vect3D posIn;
  double age;
  Vect3D posEx;
  double r;
  enum status cur_status;
  double timer;
  bool MustDelete;
  void StepR(double dt);
};

void Receptor::StepR(double dt)
{
	if (cur_status != inactive) {
    age += dt;
  }
	
  if (age >= timer) {
    cur_status = inactive; 
    age = 0;  
     
  }
}

Receptor::Receptor(Vect3D posEx_, Vect3D posIn_, double r0_, enum status cur_status_, double timer_)
{
  posEx = posEx_, posIn = posIn_;
  r=r0_; 
  cur_status = cur_status_;
  timer = timer_;
  //timer = 0;
  MustDelete = false;
  age = 0;

}






