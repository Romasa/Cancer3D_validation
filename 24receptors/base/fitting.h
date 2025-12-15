/* EGF concentrations to be consdiered 
    c = [0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.5] ng/mL
*/
/************  EGF concentration parameter fitting  ************/

// For EGF concentration = 0.05 ng/ml, # of molecules should be around 2463 => 25
 /* const int t_egf = 2700;
  int NUM_MOLECULES = 25;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 1;
*/

// For EGF concentration = 0.1 ng/ml, # of molecules should be around 4927 => 49
/*  const int t_egf = 1470;
  int NUM_MOLECULES = 49;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 1;
*/

// For EGF concentration = 0.2 ng/ml, # of molecules should be around 9855 => 99
/*  const int t_egf = 1415;
  int NUM_MOLECULES = 99;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 2;
*/

// For EGF concentration = 0.4 ng/ml, # of molecules should be around 19707/100 => 197
 const int t_egf = 1440;
int NUM_MOLECULES = 197; // number of molecules to be introduced at bulk in the start for once
int NUM_MOLECULES_R = 4;  // number of molecules to be introduced every t_egf time period


// For EGF concentration = 0.6 ng/ml, # of molecules should be around 29560 => 295
/*  const int t_egf = 1200;
  int NUM_MOLECULES = 295;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 5; 
*/

// For EGF concentration = 0.8 ng/ml, # of molecules should be around 39414 => 394
/*  const int t_egf = 1075;
  int NUM_MOLECULES = 394;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 6; 
*/
  // For EGF concentration = 1.0 ng/ml, # of molecules should be around 49267 => 493
/*  const int t_egf = 1015;
  int NUM_MOLECULES = 493;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 7;
*/

// For EGF concentration = 1.2 ng/ml, # of molecules should be around 59120 => 591
/*  const int t_egf = 965;
  int NUM_MOLECULES = 591;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 8;
*/

// For EGF concentration = 2.0 ng/ml, # of molecules should be around 98535 => 985
/* const int t_egf = 500;
  int NUM_MOLECULES = 50;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 7;
*/

// For EGF concentration = 2.5 ng/ml, # of molecules should be around 123168 => 1232
/* const int t_egf = 510;
  int NUM_MOLECULES = 50;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 9;
*/

// For EGF concentration = 3.0 ng/ml, # of molecules should be around 147802 => 1478
 /* const int t_egf = 490;
  int NUM_MOLECULES = 50;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 10;
*/

// For EGF concentration = 3.5 ng/ml, # of molecules should be around 172436 => 1724
 /* const int t_egf = 450;
  int NUM_MOLECULES = 100;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 11;
*/

// For EGF concentration = 4.0 ng/ml, # of molecules should be around 197070 => 1970
/*  const int t_egf = 440;
  int NUM_MOLECULES = 150;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 12;
*/

// For EGF concentration = 4.5 ng/ml, # of molecules should be around 221703 => 2217
 /* const int t_egf = 430;
  int NUM_MOLECULES = 150;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 13;
*/

// For EGF concentration = 5.0 ng/ml, # of molecules should be around 246337 => 2463
/*  const int t_egf = 425;
  int NUM_MOLECULES = 150;
  // Number of molecules introduced every tn timestep
  int NUM_MOLECULES_R = 15;
*/
/************  End of EGF concentration  ************/