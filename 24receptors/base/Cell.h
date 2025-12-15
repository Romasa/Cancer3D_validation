#include <list>
#include <stdio.h>
#include "Vect.h"
#include "Receptor.h"
#define allreceptors(rec) for(list <Receptor>::iterator rec = RR.begin(); rec != RR.end();rec++)
#define SQ(x) ((x) * (x))


bool Delete1(Receptor& rec) {return (rec.MustDelete == true);}


class Cell
{ 
 public:
  Cell();
  Cell(Vect3D pos_, double r0_, int type_, double TimeG01_, double cellCycle_, int erk_, gsl_rng *rng);
 
  double **node_x;
  double **node_y;
  double *node_z;
  double age;
  Vect3D pos;
  double r, R;
  int type;
  double TimeG01, cellCycle;
  bool MustDelete;
  bool MustDivide;
  bool MustDiff;
  int erk, erk_div=0;
  double erktrans;
  Vect3D v;
  Vect3D f;
  double Nrec;
  void Move(double dt, gsl_rng *rng);
  void Step(double dt, gsl_rng *rng);
  double NormFunction(double x, double delta, double mu);
  list <Receptor> RR;
};


void Cell::Move(double dt, gsl_rng *rng)
{ 
  
  if (f.x == 0) v.x = 0;
  else v.x = (1*v.x + dt*(f.x))/(1+dt*1e4);
  if (f.y == 0) v.y = 0;
  else v.y = (1*v.y + dt*(f.y))/(1+dt*1e4);
  if (f.z == 0) v.z = 0;
  else v.z = (1*v.z + dt*(f.z))/(1+dt*1e4);

  pos+=v*dt;
}

void Cell::Step(double dt, gsl_rng *rng)
{
	age += dt;
	int active_rec_count = 0;

	allreceptors(rec){
		if (rec->cur_status != inactive) {
			rec->StepR(dt);
		}
		if(rec->cur_status == active) active_rec_count++;
	}
	
//	erk = 0;
	
	double dErk_dt = K1*active_rec_count - K2*erk;	// N is the number of active receptors, calculating differential change
	erk = erk + dErk_dt*dt; 
	
	if (age >= TimeG01 && age <= TimeG01 + dt){
		int nTF = erk;
		
		if (nTF >= 50000) {
			std::cout << "cell division " << std::endl; 
			MustDiff = true;
			erk_div=erk;
		}
	} 

	if (age >= TimeG01 && age < cellCycle) {
		if (MustDiff == true){
			r = R + 0.414*R*(age-TimeG01)/(cellCycle - TimeG01);

			allreceptors(rec){
				double d = sqrt(SQ(rec->posEx.x - pos.x) + SQ(rec->posEx.y - pos.y) + SQ(rec->posEx.z - pos.z));
				
				//cout << "receptors moving! "  << std::endl;
				rec->posEx.x += R*0.112e-4*(rec->posEx.x - pos.x)/d;
				rec->posEx.y += R*0.112e-4*(rec->posEx.y - pos.y)/d;
				rec->posEx.z += R*0.112e-4* (rec->posEx.z - pos.z)/d;

				rec->posIn.x += R*0.112e-4*(rec->posIn.x - pos.x)/d;
				rec->posIn.y += R*0.112e-4*(rec->posIn.y - pos.y)/d;
				rec->posIn.z += R*0.112e-4* (rec->posIn.z - pos.z)/d;
				
			}
		}
	}

	if (age >= cellCycle) {
		if (MustDiff == true) MustDivide = true;
		
		if (MustDiff == false){
			MustDelete = true;
			allreceptors(rec) rec->MustDelete = true;
		}
	}
}


Cell::Cell(Vect3D pos_, double r0_, int type_, double TimeG01_, double cellCycle_, int erk_, gsl_rng *rng)
{
  pos=pos_; r=r0_; type = type_;
  node_x = new double*[32];
  node_y = new double*[32];
  node_z = new double[32];
  
  age = 0; 
  erktrans = 0;
  R = r0_;
  erk = erk_;
  TimeG01 = TimeG01_;
  cellCycle = cellCycle_;
  cout << "new cell "  << pos.x << " " << pos.y << std::endl;
  MustDelete = false;
  MustDivide = false;
  MustDiff = false;


  for (int x=0; x<32; x++){
    node_x[x] = new double[32];
    node_y[x] = new double[32];
  }
  for (int n=0; n<32; n++){
   for (int m=0; m<32; m++){
   node_x[n][m] = pos.x + r * sin(2. * M_PI * (double) n / 32) * sin(M_PI * (double) m/32);
   node_y[n][m] = pos.y + r * cos(2. * M_PI * (double) n / 32)* sin(M_PI * (double) m/32);
   node_z[m] = pos.z + r * cos(M_PI * (double) m/32);
   }
  }

	//set receptors
  double k=0;
  int total_rec = NRx*NRy;
	for (int i = 1; i <= NRx; i++){
    
		for (int j = 1; j <= NRy; j ++){

      double phi = acos(1-2*k/total_rec);
			double theta = M_PI * (1 + sqrt(5))* k;
			/*RR.push_back(Receptor(
      Vect3D(pos.x + (r + 0.08)*sin(2. * M_PI * (double) i / NRx) * sin(M_PI * (double) j/NRy),pos.y + (r + 0.08)*cos(2. * M_PI * (double) i / NRx) * sin(M_PI * (double) j / NRy),pos.z + (r + 0.08)*cos(M_PI * (double) j/NRy)), 
			Vect3D(pos.x + (r - 0.08)*sin(2. * M_PI * (double) i / NRx) * sin(M_PI * (double) j/NRy),pos.y + (r - 0.08)*cos(2. * M_PI * (double) i / NRx) * sin(M_PI * (double) j /NRy),pos.z + (r - 0.08)*cos(M_PI * (double) j/NRy)), 
								0.05, inactive, 0));
      */
		RR.push_back(Receptor(
		    Vect3D(
		      pos.x + (r + 0.08)*cos(theta) * sin(phi),
		      pos.y + (r + 0.08)*sin(theta) * sin(phi),
		      pos.z + (r + 0.08)*cos(phi)), 
		    Vect3D(
		      pos.x + (r - 0.08)*cos(theta) * sin(phi),
		      pos.y + (r - 0.08)*sin(theta) * sin(phi),
		      pos.z + (r - 0.08)*cos(phi)), 
		    0.05, inactive, 1));

    /*RR.push_back(Receptor( 
      Vect3D(pos.x + (r + 0.08)*sin(phi) * cos(theta),pos.y + (r + 0.08)*sin(phi) * sin(theta),pos.z + (r + 0.08)*cos(phi)), 
			Vect3D(pos.x + (r - 0.08)*sin(phi) * cos(theta),pos.y + (r - 0.08)*sin(phi) * sin(theta),pos.z + (r - 0.08)*cos(phi)), 
								0.05, inactive, 0));
                */
               k = k + 1;
    }
	}
}

Cell::Cell()
{
  pos = 0;
  r = 0;

}



double Cell::NormFunction(double x, double delta, double mu)
{
 double f=1/(delta*sqrt(2*M_PI))*exp(-(x-mu)*(x-mu)/(2*delta*delta));
 return f;
}

