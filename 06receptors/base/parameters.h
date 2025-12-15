// ----------------------------------------------------------------------------
// Simulation parameters
// ----------------------------------------------------------------------------
// Time step
double dt = 0.02;

// Total simulation time (in minutes)
const int t_num = 15*24*60/dt; 

// Time to write results
const int t_disk = 720/dt; 

// Time to return timestep info
const int t_info = 100;  

// Time to write concentration results
const int t_conc = 60/dt; // 1 hour of simulation clock


// ----------------------------------------------------------------------------
// Population-level parameters 
// ----------------------------------------------------------------------------
int NUM_TREATMENT_MOLECULES = 50;
int tn = 50;

// Friction coefficient with the medium
const double K = 50000; 

// EGF degradation OR Dissociation constant K_D
const double d_egf = 0.0007;
//const double d_egf = 0.476; // ng/ml
// Drug degradation
const double d_drg = 0.1; 

// EGF diffusion coefficient
const double diff_egf = 0.5;

// EGFR-ligand dissociation constant
//const double d_egfr = 0.01;


// Outer boundary. All molecules must exist within in a sphere of radius MAXR.
const double MAXR = 50.0; // in micrometer

// ----------------------------------------------------------------------------
// Cell characteristics and regulation
// ----------------------------------------------------------------------------
// Average phase G1 of the cell cycle
const double PG1 = 750;
const double PG1sd = 60;

// Radius of the cell and the nucleus
const double CellRad = 4;
const double NucRad = 2;

// rest length of cell-cell interaction
double const s = CellRad*2;

// maximum interaction length
double const ra = CellRad + s;

// Reaction radii
const double radd = 0.5;

// Number of receptor clusters (horizontal x vertical)
const int NRx = 2;
const int NRy = 3;


// Diffusion constants of intracellular proteins
const double diffconst = 0.1;


// Intracellular degradation rates
const double d_rec = 0.000001;
const double d_ras = 0.01;
const double d_raf = 0.1;
const double d_erk = 0.5;
const double d_TF = 0.3;

// Number of molecules in cells
const double n_ras = 150;
const double n_raf = 50; 
const double n_erk = 200;
const double n_TF = 50;

// Differential Equation constants
const double K1 = 2285.3*60;
const double K2 = 0.0017293*60;



//Drug specific parameters
const int START_DAY = (int) (99*24*60/dt); // First day to introduce drug/treatment

//For Erlotinib
const char DRUG_NAME = 'E';
//const double Ka = 0.0037578449; // Drug absorption constant
const double Ka = 0.003/2;
//const double Ke = 0.236995758; // Drug excretion constant ln(2)/half-life of erlotinib
const double Ke = 0.037;
const int STEADY_STATE = 7; // Number of days required to achieve steady-state since first dose
const float BIOAVAIL = 0.6; // Bioavailability of orally administered drug
const float DOSE = 100; // in milligrams
const int Tmax = 2; // hours to achieve the maximum concentration
const int MUL_FACTOR = 1600;


    //For Gefitinib {{Need to work on it }}
    /*
    const char DRUG_NAME = 'G';
    const double Ka = 0.00248413; // Drug absorption constant
    const double Ke = 0.23104906018664842; // Drug excretion constant ln(2)/half-life of erlotinib
    const int STEADY_STATE = 7; // Number of days required to achieve steady-state since first dose
    const float BIOAVAIL = 0.6; // Bioavailability of orally administered drug
    const float DOSE = 250; // in milligrams
    const int Tmax = 4; // hours to achieve the maximum concentration
    const int MUL_FACTOR = 1000;
*/

  //For Lapatinib {{Need to work on it }}
    /*
    const char DRUG_NAME = 'G';
    const double Ka = 0.00248413; // Drug absorption constant
    const double Ke = 0.23104906018664842; // Drug excretion constant ln(2)/half-life of erlotinib
    const int STEADY_STATE = 7; // Number of days required to achieve steady-state since first dose
    const float BIOAVAIL = 0.6; // Bioavailability of orally administered drug
    const float DOSE = 250; // in milligrams
    const int Tmax = 4; // hours to achieve the maximum concentration
    const int MUL_FACTOR = 1000;
*/