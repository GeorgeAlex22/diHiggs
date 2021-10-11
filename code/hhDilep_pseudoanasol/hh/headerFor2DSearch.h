//////////////////////////////////////////////////////////////////////////////////////////////
#define MTOPWORLD        172                 // Measured World Average (PDG) for Mtop 
#define MWWORLD           80                 // Measured World Average (PDG) for MW
#define NUMOFCOMBINATIONS  2                 // Number of combinations for topology
#define STEP               5                 // STEP for both Mtop, MW
#define NSMEARING        100                 // number of smearing
#define NMT0               0                 // Start value for Mtop 150
#define NMT             2000                 // Stop  value for Mtop  201 2001
#define NMW0               0                 // start value for MW 70
#define NMW             2000                 // Stop  value for MW 91
#define MTOPGTMW        true                 // if true then MTop should be greater than MW 
#define ECM            13000                 // Collission energy
#define MTBINS (NMT-NMT0)/STEP               // bins for MT axis
#define MWBINS (NMW-NMW0)/STEP               // bins for MW axis
#define MASSBQUARK      4.8                  // mass b quark
#define MASSLEPTON      0.00051              // mass lepton  0.10566; 
#define MASSNEUTRINO    0                    // mass neutrino       
#define NUMOFSOLUTIONS  4*NUMOFCOMBINATIONS  // number of solutions is 4 times combinations
#define STOREINITIALEVENT     false           // store initial event nice plots?
#define STORESMEAREDEVENT     false           // store smeared event nice plots?
//////////////////////////// DEBUGGING /////////////////////////////////////////////////////////
//#define NMT0                 121           // Start value for Mtop 150
//#define NMT                  121           // Stop  value for Mtop  201 2001
//#define NMW0                 594           // start value for MW 70
//#define NMW                  594           // Stop  value for MW 91
#define TESTSMEARING       false             // suppress smearing
#define PRINT_NEUTRINOS    false             // debugging - print neutrino solutions
#define PRINT_X1X2MTT      false             // debugging - print x1, x2, mttbar
#define PRINT_MTTSOLPDF    false             // debugging - print mtt, Solvability, PDF weight
////////////////////////////////////////////////////////////////////////////////////////////////
// Changes from CMS to Delphes and the opposite
#define delphesRootTree true
#define extendedRootTree false
#define superExtendedRootTree  false

