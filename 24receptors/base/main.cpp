#include <stdio.h>
#include <locale.h>
#include<chrono>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <assert.h>
#include <string.h>

#include <list>
#include "Vect.h"

// The GSL library
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_errno.h"

#ifdef __MACH__
#include <mach/mach_time.h>
#endif

#include "parameters.h"
#include "fitting.h"
#include "Cell.h"
#include "Molecule.h"
#include <vector>         // vector containers
#include <cmath>          // mathematical library
#include <fstream>        // file streams
#include <sstream>        // string streams
#include <cstdlib>        // standard library
#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code

#ifdef allcells
#undef allcells
#endif
#define allcells(c) for (list<Cell>::iterator c = CC.begin(); c != CC.end(); c++)

#ifdef allmol
#undef allmol
#endif
#define allmol(m) for (list<Molecule>::iterator m = MM.begin(); m != MM.end(); m++)

bool Delete(Cell &c) { return (c.MustDelete == true); }
bool DeleteM(Molecule &m) { return (m.MustDelete == true); }

using namespace std;
list<Cell> CC;
list<Molecule> MM;

int numb = 0;
int numb1 = 0;
int numb2 = 0;
int numb3 = 0;
int numb4 = 0;
int numb5 = 0;
int numb6 = 0;
int numb7 = 0;
int numb8 = 0;
int numb9 = 0;
int numb10 = 0;
int numbMol = 0;
int numbRec = 0;
int numbRecIn = 0;
int numbRecBl = 0;
int numbDrugMol = 0;
double allerk = 0;
int allgencells = 0;
double conc = 0; // drug concentration in plasma at any time dt 
int ti = 0; // time to represent the start of new day ti % 60*24/dt
int daycount = 0; // number of days since the first administration of drug
int latestDrugMolCount = 0;

// Declaring functions
void initialize(gsl_rng *rng);
void set_particle();
void write_particle1_vtk(int, int);
void write_particle_vtk(int, int);
void write_molecule_vtk(int, int);
void write_treatment_molecule_vtk(int, int);
void write_receptor_vtk(int, int);
void write_receptorbl_vtk(int, int);
void write_receptorIn_vtk(int, int);
void countCells(int);
void generateMolecules(gsl_rng *rng);
void generateTreatmentMolecules(gsl_rng *rng);
void introduceMolecules(gsl_rng *rng);
void introduceTreatmentMolecules(gsl_rng *rng);
void introduceDrug(int t, gsl_rng *rng);
void computeDivision(int t, gsl_rng *rng);
void molecularDegradation(double);
void computeForces();
void MolecularMotion(double, gsl_rng *rng);
double Fn(double h, double h0, bool &success);
void computeMotion(double, gsl_rng *rng);
void write_data(int);
void write_conc(int);
void write_cell_division_data(int time, int erk_div, int erk);
void computeReaction();
void activeRec(gsl_rng *rng);
void intraCell(double, gsl_rng *rng);
int countTreatmentMolecules();

// Main function
int main()
{

  auto start = std::chrono::high_resolution_clock::now();
  srand(time(NULL));

  gsl_rng *rng;
#ifdef __MACH__
  uint64_t seed;
  seed = mach_absolute_time();
#elif __linux__
  timespec t;

  clock_gettime(CLOCK_REALTIME, &t);

  unsigned long seed = (unsigned long)(1e9 * t.tv_sec + t.tv_nsec);
  srand48(seed);
#endif
  rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, seed);

  initialize(rng);

  generateMolecules(rng);

  countCells(0); // Counting total number of cells, receptors and molecules(both drug and EGF)

  cout << "Parameters loaded! ready to begin" << endl;
  cout << "=================================" << endl;

  cout << endl;
  cout << "starting simulation" << endl;

  for (int t = 1; t <= t_num; ++t)
  { // loop of computation

    countCells(t);
    MolecularMotion(dt, rng); // Moving molecules

    if (t % t_egf == 0)
    {
      introduceMolecules(rng);
    }
    
    if (t % tn == 0 && t >= START_DAY) // Introducing treatment molecules after n simulation days
    //if(t % tn == 0 && numb1 >= 50)
    {
      //if (numb1 == 50) cout << "Treatment starts .... \n";
      introduceDrug(t, rng);
    }
    activeRec(rng);
    molecularDegradation(dt);

    computeForces();
    computeMotion(dt, rng);

    intraCell(dt, rng);
    computeDivision(t, rng);


    if (t % t_disk == 0)
    { // write only each t_disk

      write_particle_vtk(t, 1);
      write_particle1_vtk(t, 1);
      write_molecule_vtk(t, numbMol);
      write_treatment_molecule_vtk(t, numbDrugMol);
      write_receptor_vtk(t, numbRec);
      write_receptorbl_vtk(t, numbRecBl);
      write_receptorIn_vtk(t, numbRecIn);
    }

    /// Report end of time step

    if (t % t_info == 0)
    {
      write_data(t);    
      cout << "completed time step " << t << " in [1, " << t_num << "]" << endl;
    }

    if( t % t_conc == 0) // write concentration data every 30 minutes
    {
      write_conc(t);
    }
    if(numb >= 200){
      auto end = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
      std::cout << "Runtime: " << (duration.count())/1000.0 << " seconds" << std::endl;
      std::cout << "Total number of cells generated " << allgencells << endl;
      exit(0);
    }
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Runtime: " << (duration.count())/1000.0 << " seconds" << std::endl;
  std::cout << "Total number of cells generated " << allgencells << endl;
}

void initialize(gsl_rng *rng)
{

  int ignore; // ignore return value of system calls
  // ignore = system("mkdir -p vtk_fluid"); // create folder if not existing
  ignore = system("mkdir -p vtk_particle");  // create folder if not existing
  ignore = system("rm -f data.txt");         // delete file if existing
  ignore = system("rm -f time-conc.txt"); // delete time-conc.txt file if existing
  ignore = system("rm -f vtk_particle/*.*"); // delete files if existing
  cout << ignore << endl;
  /// Set cells

  double Time1 = 0;
  // Time1 = 600 + 300*(double) rand()/RAND_MAX;
  Time1 = 720;

  CC.push_back(Cell(Vect3D(0, 0, 0), CellRad, 1, Time1, 2 * Time1, 0, rng));
  allgencells++;
  return;
}

void generateMolecules(gsl_rng *rng)
{
  int n = 0;

  // cubic generation
  while (n < NUM_MOLECULES)
  {
    double c1 = 0;
    double c2 = 0;
    double c3 = 0;
    double times = 0;

    c1 = -50 + 100 * (double)rand() / RAND_MAX;
    c2 = -50 + 100 * (double)rand() / RAND_MAX;
    c3 = -50 + 100 * (double)rand() / RAND_MAX;

    times = gsl_ran_exponential(rng, 1 / d_egf);
    //cout << times << endl;

    if (sqrt(SQ(c1) + SQ(c2) + SQ(c3)) < MAXR)
    {
      MM.push_back(Molecule(Vect3D(c1, c2, c3), times));
      n++;
    }
  }
}

void generateTreatmentMolecules(gsl_rng *rng)
{

  cout << "generating " << NUM_TREATMENT_MOLECULES << " molecules" <<endl;
  int n = 0;
  // cubic generation
  while (n < NUM_TREATMENT_MOLECULES)
  {
    double c1 = 0;
    double c2 = 0;
    double c3 = 0;
    double times = 0;

    c1 = -50 + 100 * (double)rand() / RAND_MAX;
    c2 = -50 + 100 * (double)rand() / RAND_MAX;
    c3 = -50 + 100 * (double)rand() / RAND_MAX;

    times = gsl_ran_exponential(rng, 1 / d_drg);

    if (sqrt(SQ(c1) + SQ(c2) + SQ(c3)) < MAXR)
    {
      MM.push_back(Molecule(Vect3D(c1, c2, c3), times, 1));
      n++;
    }
  }
}

void introduceMolecules(gsl_rng *rng)
{

  int n = 0;

  // spheric generation
  /*
    while (n < 20){
      double rad = 0;
            double phi = 0;
      double psi = 0;
      double times = 0;

      rad = 140 +(double) rand()/RAND_MAX*150;
      phi = (double) rand()/RAND_MAX*2*M_PI;
      psi = (double) rand()/RAND_MAX*M_PI;

      double node1_x = 0 + rad * sin(phi) * sin(psi);
                  double node1_y = 0 + rad * cos(phi)* sin(psi);
                  double node1_z = 0 + rad * cos(psi);
      times = 50 + 20*(double) rand()/RAND_MAX;

                  MM.push_back(Molecule(Vect3D(node1_x,node1_y,node1_z), times));

      n++;


    }
   */
  // cubic generation
  while (n < NUM_MOLECULES_R)
  {
    double c1 = 0;
    double c2 = 0;
    double c3 = 0;
    double times = 0;

    c1 = -50+ 100 * (double)rand() / RAND_MAX;
    c2 = -50 + 100 * (double)rand() / RAND_MAX;
    c3 = -50 + 100 * (double)rand() / RAND_MAX;

    times = gsl_ran_exponential(rng, 1 / d_egf);
    //cout << "Times = " <<times <<endl;
    if (sqrt(SQ(c1) + SQ(c2) + SQ(c3)) < MAXR && sqrt(SQ(c1) + SQ(c2) + SQ(c3)) >= MAXR / 3)
    {
      MM.push_back(Molecule(Vect3D(c1, c2, c3), times));
      n++;
    }
  }
}

void introduceTreatmentMolecules(gsl_rng *rng)
{
  int n = 0;
  int calculatedConc = (int) (conc*MUL_FACTOR);
  int diff = calculatedConc - latestDrugMolCount;
  
  if(diff > 0)
  {
    n = 0;
    while (n < diff)
    {
      double c1 = 0;
      double c2 = 0;
      double c3 = 0;
      double times = 0;

      c1 = -50 + 100 * (double)rand() / RAND_MAX;
      c2 = -50 + 100 * (double)rand() / RAND_MAX;
      c3 = -50 + 100 * (double)rand() / RAND_MAX;

      times = gsl_ran_exponential(rng, 1 / d_drg);
      if (sqrt(SQ(c1) + SQ(c2) + SQ(c3)) < MAXR && sqrt(SQ(c1) + SQ(c2) + SQ(c3)) >= MAXR / 3)
      {
        MM.push_back(Molecule(Vect3D(c1, c2, c3), times, 1));
        n++;
      }
    }
    latestDrugMolCount += diff;
    //cout << "Added " << diff << " drug molecules. Total Drug count = "<< latestDrugMolCount << endl;;
    
  }
  /*
  else if (diff < 0)
  {
    cout << "Deleted " << diff << " drug molecules. Total Drug count = "<< latestDrugMolCount +diff<< endl;
    n = -1*diff;
    allmol(m)
    {
      if(n <= 0) break; //Stop looping over if required number of drug molecules are marked to be deleted
    
      m->MustDelete;
      n--;
    }
    latestDrugMolCount += diff;
  }
  */
  
  //MM.remove_if(DeleteM);
}

void introduceDrug(int t, gsl_rng *rng)
{
  ti = t%((int) (60*24/dt));
  double dose = DOSE;
  if( t > (int) (10*24*60/dt)) dose = DOSE/2;

  if(ti < Tmax*60/dt) 
  {
    conc = conc + dt*(Ka*DOSE - Ke*conc);
  }
  else
    conc = conc + dt*(-Ke*conc);

  introduceTreatmentMolecules(rng);
}

void MolecularMotion(double dt, gsl_rng *rng)
{
  // random motion
  allmol(m)
  {
      m->Move(dt, rng);
    // prevent infiltration into cells
    allcells(c)
    {
      if (sqrt(SQ(m->pos.x - c->pos.x) + SQ(m->pos.y - c->pos.y) + SQ(m->pos.z - c->pos.z)) <= c->r + m->r)
      {
        double d = sqrt(SQ(m->pos.x - c->pos.x) + SQ(m->pos.y - c->pos.y) + SQ(m->pos.z - c->pos.z));
        m->pos.x += 0.1 * (m->pos.x - c->pos.x) / d;
        m->pos.y += 0.1 * (m->pos.y - c->pos.y) / d;
        m->pos.z += 0.1 * (m->pos.z - c->pos.z) / d;
      }
    }
  }
}

void molecularDegradation(double dt)
{
  allmol(m) m->Step(dt);
  allmol(m){ 
    if(m->MustDelete == true && m->type == 1) 
    {
      latestDrugMolCount--;
    }
  }
  MM.remove_if(DeleteM);
}

void computeDivision(int t, gsl_rng *rng)
{
  allcells(c)
  {

    if (c->MustDivide == true)
    {
      int erk1 = (int) c->erk / 2;
      c->MustDelete = true;
      c->RR.remove_if(Delete1);
      double dir = 0;
      double dir2 = 0;

      dir = (double)rand() / RAND_MAX * 2 * M_PI;
      dir2 = (double)rand() / RAND_MAX * M_PI;

      Vect3D dpos = Vect3D(CellRad * cos(dir) * cos(dir2), CellRad * sin(dir) * cos(dir2), CellRad * sin(dir2));

      double TimeG = PG1;
      TimeG += gsl_ran_gaussian_ziggurat(rng, PG1sd);
      

      CC.push_back(Cell(Vect3D(c->pos.x - dpos.x, c->pos.y - dpos.y, c->pos.z - dpos.z), CellRad, 1, TimeG, 2 * TimeG, erk1, rng));

      TimeG = PG1;
      TimeG += gsl_ran_gaussian_ziggurat(rng, PG1sd);

      CC.push_back(Cell(Vect3D(c->pos.x + dpos.x, c->pos.y + dpos.y, c->pos.z + dpos.z), CellRad, 1, TimeG, 2 * TimeG, erk1, rng));
      allgencells++;
    }
  }
  CC.remove_if(Delete);
}

void computeMotion(double dt, gsl_rng *rng)
{
  allcells(c)
  {

    c->Move(dt, rng);
  }
}

void intraCell(double dt, gsl_rng *rng)
{
  allcells(c)
  {
    c->Step(dt, rng);
  }
}

double Fn(double h, double h0, bool& success) {
    if ( h >= h0-2*NucRad || h < 0) success = true; 
    else  
      {	success = true; 
	return 0;
      }  // 0- if ()-false
    return (h<h0)?K * ((h0-h) / (h-h0+2*NucRad)): 0 ; }

void computeForces(){
  allcells(c)   {c->f=0;} 
	allcells(c){
	  bool success;
	  allcells(cp){
      if (cp != c){
        double d = sqrt((cp->pos.x - c->pos.x)*(cp->pos.x - c->pos.x)+(cp->pos.y - c->pos.y)*(cp->pos.y - c->pos.y) + (cp->pos.z - c->pos.z)*(cp->pos.z - c->pos.z));
        double F = Fn(d, c->r + cp->r, success);
        if (d <= c->r + cp->r){
          c->f.x -= (F)*(cp->pos.x-c->pos.x)/d;
          c->f.y -=  (F)*(cp->pos.y-c->pos.y)/d;
          c->f.z -=  (F)*(cp->pos.z-c->pos.z)/d;

          cp->f.x +=  (F)*(cp->pos.x-c->pos.x)/d; 
          cp->f.y +=  (F)*(cp->pos.y-c->pos.y)/d;
          cp->f.z +=  (F)*(cp->pos.z-c->pos.z)/d;
        }
      }
	  }
  }
}

void countCells(int time)
{
  numb = 0;
  numbMol = 0;
  numbDrugMol = 0;
  numbRec = 0;
  numbRecIn = 0;
  numbRecBl = 0;
  numb1 = 0;
  numb2 = 0;
  numb3 = 0;
  numb4 = 0;
  numb5 = 0;
  numb6 = 0;
  numb7 = 0;
  numb8 = 0;
  numb9 = 0;
  numb10 = 0;
  allerk = 0;
  allcells(c)
  {
    numb++;
    allerk += c->erktrans;
    if (c->type == 1)
      numb1++;
    if (c->type == 2)
      numb2++;
    if (c->type == 3)
      numb3++;
    if (c->type == 4)
      numb4++;
    if (c->type == 5)
      numb5++;
    if (c->type == 6)
      numb6++;
    if (c->type == 7)
      numb7++;
    if (c->type == 8)
      numb8++;
    if (c->type == 9)
      numb9++;
    if (c->type == 10)
      numb10++;

    for (list<Receptor>::iterator rec = c->RR.begin(); rec != c->RR.end(); rec++)
    {
      if (rec->cur_status == active)
        numbRec++;
      else if (rec->cur_status == blocked)
        numbRecBl++;
      else {
        numbRecIn++;
      }
    }
  }

  // Counting all molecules including drug and EGF
  allmol(m) (m->type == 0)? numbMol++ : numbDrugMol++; 
}

// In this function, EGF molecule activates Receptor, Drug molecule blocks receptor
void activeRec(gsl_rng *rng)
{
  allcells(c)
  { 
    for (list<Receptor>::iterator rec = c->RR.begin(); rec != c->RR.end(); rec++)
    {
      allmol(m)
      {
        double d = sqrt(SQ(rec->posEx.x - m->pos.x) + SQ(rec->posEx.y - m->pos.y) + SQ(rec->posEx.z - m->pos.z));

        if(d <= 0.5) // increased EGF affinity
        { 
          if(m->type == 0 && rec->cur_status == inactive)//EFG Molecule -> activates receptor
          {
            rec->cur_status = active;
            cout << endl << "receptor activation" << endl << endl;
            rec->timer = gsl_ran_exponential(rng, 1 / d_rec);
            m->MustDelete = true;
          }
          //TKI-inhibitor blocking respective EGFR
          else if((m->type == 1 && DRUG_NAME == 'E' && rec->cur_status != blocked) ||
            (m->type == 1 && DRUG_NAME == 'L' && rec->cur_status == inactive) ||
            (m->type == 1 && DRUG_NAME == 'G' && rec->cur_status == active))//Treatment molecule -> blocks receptor
          {
            rec->cur_status = blocked;
            m->MustDelete = true;
            m->age = m->time;
            cout << endl <<endl<< "A receptor is blocked" << endl << endl;
            rec->timer = gsl_ran_exponential(rng, 1 / (0.0000000000000001*d_rec)); 
          }

        } 

      }     
    }
  }
  allmol(m){ 
    if(m->MustDelete == true && m->type == 1) 
    {
      latestDrugMolCount--;
    }
  }
  MM.remove_if(DeleteM);
}

int countTreatmentMolecules()
{
  numbDrugMol = 0;

  for (list<Molecule>::iterator m = MM.begin(); m != MM.end(); m++)
  {
    if (m->type == 1)
      numbDrugMol++;
  }
  return numbDrugMol;
}


void write_data(int time)
{
  string output_filename("data.txt");
  ofstream output_file;

  output_file.open(output_filename.c_str(), fstream::app);

  output_file << time  << " "; // time step
  output_file << numb1 << " ";     // number of cells
  output_file << numbMol << " "; // number of molecules
  output_file << numbDrugMol << " ";           // number of drug molecules
  output_file << allerk << " ";               // count of erk inside the cell
  output_file << numbRec << " "; // Number of Active receptors
  output_file << numbRecIn << " "; // Number of Inactive receptors
  output_file << numbRecBl << "\n"; // Number of Blocked Receptors
  
  output_file.close();
/*
  if(numb1 == 0)
  {
    cout << "All cells died. Aborting!";
    exit(0);
  }
*/
  return;
}

void write_conc(int time)
{
  string output_filename("time-conc.txt");
  ofstream output_file;

  output_file.open(output_filename.c_str(), fstream::app);
  output_file << (int)time*dt/60 << " "; 
  output_file << conc << " ";
  output_file << latestDrugMolCount << " ";
  output_file << numbDrugMol << "\n";
  output_file.close();
  return;

}

void write_particle_vtk(int time, int Type)
{

  /// Create filename
  double numbz = 0;
  if (Type == 1)
    numbz = numb1;

  stringstream output_filename;
  output_filename << "vtk_particle/particle_t" << Type << " " << time * 1/t_disk << ".vtk";
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.str().c_str());

  /// Write VTK header

  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions
  if (numbz == 0)
  {
    output_file << "POINTS " << 1 << " float\n";

    output_file << 0 << " " << 0 << " " << 0 << "\n";
  }
  else
  {
    output_file << "POINTS " << numbz << " float\n";
    allcells(c)
    {
      if (c->type == Type)
      {

        output_file << c->pos.x << " " << c->pos.y << " " << c->pos.z << "\n";
      }
    }
  }

  /// Close file

  output_file.close();

  return;
}

void write_particle1_vtk(int time, int Type)
{

  /// Create filename
  double numbz = 0;
  if (Type == 1)
    numbz = numb1;
  if (Type == 2)
    numbz = numb2;
  if (Type == 3)
    numbz = numb3;
  if (Type == 4)
    numbz = numb4;
  if (Type == 5)
    numbz = numb5;
  if (Type == 6)
    numbz = numb6;
  if (Type == 7)
    numbz = numb7;
  if (Type == 8)
    numbz = numb8;
  if (Type == 9)
    numbz = numb9;
  if (Type == 10)
    numbz = numb10;
  stringstream output_filename;
  output_filename << "vtk_particle/particle1_t" << Type << " " << time *1/t_disk << ".vtk";
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.str().c_str());

  /// Write VTK header

  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions

  output_file << "POINTS " << SQ(32) * numbz << " float\n";
  allcells(c)
  {
    if (c->type == Type)
    {

      for (int n = 0; n < 32; n++)
      {
        for (int m = 0; m < 32; m++)
        {
          c->node_x[n][m] = c->pos.x + c->r * sin(2. * M_PI * (double)n / 32) * sin(M_PI * (double)m / 32);
          c->node_y[n][m] = c->pos.y + c->r * cos(2. * M_PI * (double)n / 32) * sin(M_PI * (double)m / 32);
          c->node_z[m] = c->pos.z + c->r * cos(M_PI * (double)m / 32);
        }
      }

      for (int m = 0; m < 32; ++m)
      {
        for (int n = 0; n < 32; ++n)
        {
          output_file << c->node_x[n][m] << " " << c->node_y[n][m] << " " << c->node_z[m] << "\n";
        }
      }
    }
  }

  output_file << "POLYGONS " << 992 * numbz << " " << 992 * numbz * 5 << "\n";
  for (int p = 0; p < numbz; p++)
  {
    for (int n = 0; n < 31; ++n)
    {
      for (int j = 0; j < 32; j++)
      {
        output_file << 4 << " " << j + 32 * n + p * 1024 << " " << j + 32 * n + 32 + p * 1024 << " " << j + 32 * n + 33 + p * 1024 << " " << j + 32 * n + 1 + p * 1024 << "\n";
      }
    }
  }
  /// Close file

  output_file.close();

  return;
}

void write_molecule_vtk(int time, int numbMol)
{
  int numbz = numbMol;
  stringstream output_filename;
  output_filename << "vtk_particle/molecule_t"     << " " << time * 1/t_disk << ".vtk";
  ofstream output_file;

  /// Open file
  output_file.open(output_filename.str().c_str());

  /// Write VTK header
  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions
  if (numbz == 0)
  {
    output_file << "POINTS " << 1 << " float\n";
    output_file << 0 << " " << 0 << " " << 0 << "\n";
  }
  else
  {
    output_file << "POINTS " << numbz << " float\n";
    allmol(m)
    {
      if (m->type == 0)
        output_file << m->pos.x << " " << m->pos.y << " " << m->pos.z << "\n";
    }
  }

  /// Close file
  output_file.close();

  return;
}

void write_treatment_molecule_vtk(int time, int numbDrugMol)
{
  /// Create filename
  int numbz = numbDrugMol;
  stringstream output_filename;
  output_filename << "vtk_particle/treatment_molecule_t"   << " " << time *1/t_disk << ".vtk";
  ofstream output_file;

  /// Open file
  output_file.open(output_filename.str().c_str());

  /// Write VTK header
  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions
  output_file << "POINTS " << numbz << " float\n";
  allmol(m)
  {
    if (m->type == 1)
      output_file << m->pos.x << " " << m->pos.y << " " << m->pos.z << "\n";
  }

 /// Close file
  output_file.close();

  return;
}

void write_receptor_vtk(int time, int numbRec)
{

  /// Create filename

  int numbz = numbRec;

  // int numbz = 20;
  stringstream output_filename;
  output_filename << "vtk_particle/receptor_t" << " " << time*1/t_disk << ".vtk";
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.str().c_str());

  /// Write VTK header

  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions
  if (numbz == 0)
  {
    output_file << "POINTS " << 1 << " float\n";

    output_file << 100 << " " << 100 << " " << 100 << "\n";
  }
  else
  {
    output_file << "POINTS " << 2 * numbz << " float\n";
    allcells(c)
    {
      for (list<Receptor>::iterator rec = c->RR.begin(); rec != c->RR.end(); rec++)
      {
        if (rec->cur_status == active)
        {
          output_file << rec->posIn.x << " " << rec->posIn.y << " " << rec->posIn.z << "\n";
          output_file << rec->posEx.x << " " << rec->posEx.y << " " << rec->posEx.z << "\n";
        }
      }
    }
  }

  /// Close file

  output_file.close();

  return;
}

void write_receptorbl_vtk(int time, int numbRec)
{
  int numbz = numbRec;

  stringstream output_filename;
  output_filename << "vtk_particle/receptorbl_t" << " " << time * 1/t_disk << ".vtk";
  ofstream output_file;

  output_file.open(output_filename.str().c_str());

  /// Write VTK header
  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions
  if (numbz == 0)
  {
    output_file << "POINTS " << 1 << " float\n";

    output_file << 100 << " " << 100 << " " << 100 << "\n";
  }
  else
  {
    output_file << "POINTS " << 2 * numbz << " float\n";
    allcells(c)
    {
      for (list<Receptor>::iterator rec = c->RR.begin(); rec != c->RR.end(); rec++)
      {
        if (rec->cur_status == blocked)
        {
          output_file << rec->posIn.x << " " << rec->posIn.y << " " << rec->posIn.z << "\n";
          output_file << rec->posEx.x << " " << rec->posEx.y << " " << rec->posEx.z << "\n";
        }
      }
    }
  }

  /// Close file

  output_file.close();

  return;
}

void write_receptorIn_vtk(int time, int numbRecIn)
{

  /// Create filename

  int numbz = numbRecIn;

  // int numbz = 20;
  stringstream output_filename;
  output_filename << "vtk_particle/receptorIn_t"   << " " << time * 1/t_disk<< ".vtk";
  ofstream output_file;

  /// Open file

  output_file.open(output_filename.str().c_str());

  /// Write VTK header

  output_file << "# vtk DataFile Version 3.0\n";
  output_file << "particle_state\n";
  output_file << "ASCII\n";
  output_file << "DATASET POLYDATA\n";

  /// Write node positions
  if (numbz == 0)
  {
    output_file << "POINTS " << 1 << " float\n";
    output_file << 100 << " " << 100 << " " << 100 << "\n";
  }
  else
  {
    output_file << "POINTS " << 2 * numbz << " float\n";
    allcells(c)
    {
      for (list<Receptor>::iterator rec = c->RR.begin(); rec != c->RR.end(); rec++)
      {
        
        if (rec->cur_status == 0)
        {
          output_file << rec->posIn.x << " " << rec->posIn.y << " " << rec->posIn.z << "\n";
          output_file << rec->posEx.x << " " << rec->posEx.y << " " << rec->posEx.z << "\n";
        }
      }
    }
  }

  /// Close file

  output_file.close();

  return;
}
