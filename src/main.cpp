#include <cstdlib>
#include <ctime>
#include <boost/filesystem.hpp>
#include "utils.h"
#include "Transcript.h"
#include "Individual.h"

static uint POP_SIZE = 1;

int main(int argc, char** argv) {
    if (argc==1) std::cout << "torsionEvol path/to/params.ini path/to/working/folder path/to/environment.dat[ IndelPoissonMean InvP ]" << std::endl;
    if (argc==1) return EXIT_FAILURE;
    if (argc>5) Individual::set_mutation( atof(argv[4]), atof(argv[5]));
    
    srand(time(NULL));

    // define the input/output directories
    Params* params = readIni(argv[1]);
    boost::filesystem::create_directories(argv[2]);    // output
    argv[1][strlen(argv[1])-10] = 0;     //removes "params.ini" in string "path/to/params.ini"
    const char* pth = argv[1];            // by setting "p" to 0 (end of string)

    // read files. The "file" types are vectors of structs
    prot_file prot;
    TTS_file tts;
    TSS_file tss;
    GFF_file gff_df_raw;
    vector<double> env;
    readProt(prot, pth + params->BARR_FIX);
    readTTS(tts, pth + params->TTS); 
    readTSS(tss, pth + params->TSS); 
    readGFF(gff_df_raw, pth + params->GFF);
    readEnv(env, argv[3]);
    Individual::set_simulation_params(params);
    Individual::set_target_envir(env);

    // =========================== Transcripts definition ==============================
    // Map of transciption units with the list of tts belonging to TU. [ TU ]  = (tts1, tts2, ... )
    map< uint , vector<uint> > TU_tts = get_TU_tts(tss);
    vector<Transcript> tr;
    for (TSS_t this_TSS : tss)
    {
        // Get the list of TTS that are in the same TU of this tss_id
        // TU_tts ex : { <0: {0, 1}>, <1: {2}> }  ==>  0 -> [0, 1]
        vector<uint> this_TU_tts = TU_tts[this_TSS.TUindex];

        float proba_rest = 1.0;
        uint k = 0; // Start scanning TTS.dat lines from the TUindex of this_TU_tts[0]
        while (proba_rest > 0 and k<this_TU_tts.size())
        {
            TTS_t this_TTS = tts[this_TU_tts[k]];
            tr.push_back(Transcript(this_TSS.TUindex,
                                    this_TSS.TSS_pos, 
                                    this_TTS.TTS_pos, 
                                    (this_TSS.TUorient=='+')-(this_TSS.TUorient=='-'), 
                                    this_TSS.TSS_strength*this_TTS.TTS_proba_off*proba_rest,
                                    params->DELTA_X)
                        );
            proba_rest *= (1.0 - this_TTS.TTS_proba_off);
            k++;
        }
    }

    // ====================== Topological barriers ============================
    std::vector<DNApos> Barr_fix; // Get the fixed topo barriers in a vector
    std::transform(prot.begin(), prot.end(), back_inserter(Barr_fix), 
                    [params](prot_t const& x) { return int(x.prot_pos/params->DELTA_X); });    

    // a bit of control on what happens
    std::cout << std::endl << "Starting with individuals with the following genome:" << std::endl;
    for (Transcript t : tr) 
    {
        std::cout << "TU n°" << t.TUindex_ <<  ": ";
        std::cout << t.TSS_ << '-' << t.TTS_ << " (length = " << t.size_ << "), strand ";
        std::cout << t.s_ << " with rate " << t.r_ << std::endl;
    }
    std::cout << std::endl;

    // Creation of the population:
    std::cout << "Creating " << POP_SIZE << " individuals..." << std::endl;
    vector<Individual*> population;
    population.reserve(POP_SIZE);
    uint genome_size = get_genome_size(gff_df_raw);
    for (uint i=0; i<POP_SIZE; ++i)
        population.push_back(new Individual(genome_size, tr, Barr_fix));

    for (Individual* indiv : population)
        indiv->estimate_exression();

    // Cleaning
    std::cout << std::endl <<  "Simulation completed. Deleting individuals..." << std::endl;
    for (Individual* indiv : population)
        delete indiv;
    delete params;

    // Exiting
    return(EXIT_SUCCESS);
}
