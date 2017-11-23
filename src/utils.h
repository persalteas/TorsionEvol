#include <iostream>
#include <vector>
#include "IniReader.h"
#include <boost/algorithm/string.hpp>
using namespace std;

/* A struct to store parameters of the .ini file. These are default values.
   We need to get the real values from the params.ini file. */
typedef struct {
    // INPUTS
    string GFF, TSS, TTS, BARR_FIX;
    // GLOBAL
    double  J_0 = 0.1;                          // The basal local flux of supercoiling (in each delta_x)
    int     D = 800;                            // The effective diffusivity of supercoiling along DNA (bp²/s)
    double  m = 2.20;                           // for k calculation
    double  sigma_t = -0.04, epsilon = 0.01;    // from Meyer article
    // SIMULATION
    double  SIGMA_0 = -0.02;
    int     DELTA_X = 60;                       // Spatial discretisation in bp
    int     DELTA_T = 2;                        // Time step corresponding to 25~30 nt/s
    int     RNAPS_NB = 8;                       // Number of RNAP
    int     ITERATIONS_NB = 10000;              // Number of iterations (s)
    int     OUTPUT_STEP = 1000;                 // Time interval at which a graphical and text output is given
    double  GYRASE_CONC = 0.0;                  // Gyrase concentration (microM) 0.1
    double  TOPO_CONC = 0.0;                    // Topoisomerase I concentration (microM) 0.041
    double  GYRASE_EFFICIENCY = 0.09;           // Gyrase e fficiency
    double  TOPO_EFFICIENCY = 0.31;             // Guess what.
    double  GYRASE_CTE = 0.001;                 // Gyrase constant
    double  TOPO_CTE = 0.0005;                  // Topo constant
    int     k_GYRASE = 50;
    double  x0_GYRASE = 0.016;
    int     k_TOPO = 80;
    double  x0_TOPO = -0.04;
} Params;

typedef struct {
    const char*   prot_name;
    uint    prot_pos;
} prot_t;

typedef struct {
    uint    TUindex;
    char    TUorient;
    uint    TTS_pos;
    double  TTS_proba_off;
} TTS_t;

typedef struct {
    uint    TUindex;
    char    TUorient;
    uint    TSS_pos;
    double  TSS_strength;
} TSS_t;

Params *readIni(const char *cfgFile);
vector<prot_t>* readProt(string protFile);
vector<TSS_t>* readTSS(string TSSFile);
vector<TTS_t>* readTTS(string TTSFile);