#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <vector>
#include <valarray>
#include <cctype>
#include <cmath>
#include <map>
#include <algorithm>
#include "IniReader.h"
#include <boost/algorithm/string.hpp>
#include <boost/multi_array.hpp>
using namespace std;

/* Vocabulary:
    TU = Transcription Unit
    TSS = Transcription Starting Site
    TTS = Transcription Termination Site
    RNAP = RNA Polymerase

    prot.dat =  proteins that make topological barriers that prevent supercoiling to spread.
                Supercoiling is constant between 2 barriers.
                RNAPs are barriers too.
    TSS.dat =   transcription starting sites.
                n° of the TU  |  strand   |   position   |   strength
    TTS.dat =   transcription termination sites.
                n° of the TU  |  strand   |   position   |   strength
    tousgenesidentiques.gff = annotations of TUs
*/


/* ===========================================================
                  New types (and file types)
   =========================================================== */
typedef unsigned int DNApos;
typedef boost::multi_array<int, 3> array3;
typedef array3::index arr_index;

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
    string  prot_name;
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

typedef struct {
    string  seqname;    // name of the chromosome or scaffold
    string  source;     // name of the program that generated this feature, or the data source (database or project name) 
    string  feature;    // feature type name, e.g. Gene, Variation, Similarity
    uint    start;      // Start position of the feature, with sequence numbering starting at 1
    uint    end;        // End position of the feature, with sequence numbering starting at 1
    double  score;
    char    strand;     // defined as + (forward) or - (reverse)
    int     frame;      // One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on
    string  attribute;  // A semicolon-separated list of tag-value pairs, providing additional information about each feature
} GFF_t;                // (Format like GFF file format in "Ensembl" databases)

typedef vector<prot_t> prot_file;
typedef vector<TSS_t> TSS_file;
typedef vector<TTS_t> TTS_file;
typedef vector<GFF_t> GFF_file;


/* ===========================================================
                    Declarations of functions
   =========================================================== */

Params*     readIni(const char *cfgFile);
void        readProt(prot_file& data, string protFile);
void        readTSS(TSS_file& data, string TSSFile);
void        readTTS(TTS_file& data, string TTSFile);
void        readGFF(GFF_file& data, string GFFFile);

ostream &operator<<(ostream &stream, TSS_t const &s);
ostream &operator<<(ostream &stream, TTS_t const &s);
ostream &operator<<(ostream &stream, GFF_t const &s);
ostream &operator<<(ostream &stream, prot_t const &s);


uint	get_genome_size(GFF_file& gff_df);
map< uint , vector<uint> > get_TU_tts_pos(TSS_file& tss, TTS_file& tts);
map< uint , vector<uint> > get_TU_tts(TSS_file& tss);
double  f_prob_unhooked_rate(double sum_Kon, int DELTA_T, size_t RNAPs_unhooked_nbr);
void    random_choice(vector<int>& result, const vector<int>& array, uint n, const vector<double>& probs);
void    calc_sigma(vector<double>& Barr_sigma, Params* params);

/* ===========================================================
                     Definition of templates
   =========================================================== */

/* Prints a vector's content, elements separated by spaces */
template<typename file_type>
void    display_vector(file_type& v){
    for (size_t n = 0; n < v.size(); n++)
    	cout << v[n] << " ";
  	cout << endl;
}


/* prints a valarray content */
template<typename T>
void display_array(T& v)
{
    cout << "{ ";
    for (size_t i = 0; i < v.size(); i++)
        cout << v[i] << " ";
    cout << "}" << endl;
}


/* sums the elements of a vector. */
template<typename T>
T	vector_sum(vector<T>& V) {
	T sum = 0; for (size_t i = 0; i < V.size(); i++) sum += V[i];
	return sum;
}



/* returns the index in Barr_pos where to place a TSS or a new barrier */
template<typename T>
void searchsorted(vector<uint>& result, vector<uint>& Barr_pos, vector<T>&TSS_pos)
{
	uint index;
	for (auto tss : TSS_pos)
	{
		index = find_if(Barr_pos.begin(), Barr_pos.end(), 
						[tss](T barr)->bool{ return barr>tss; }) - Barr_pos.begin();
		result.push_back(index);
	}
}



#endif //UTILS_H_