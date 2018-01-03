#include <cstdlib>
#include <ctime>
#include <boost/filesystem.hpp>
#include "utils.h"
#include "Transcript.h"
#include "RNAP.h"

using namespace std;

static float RNAPs_genSC = 0.1;

int main(int argc, char** argv) {
    if (!argc) cout << "torsionEvol path/to/params.ini path/to/working/folder";
    Params* params = readIni(argv[1]);
    srand(time(NULL));

    // define the input/output directories
    boost::filesystem::create_directories(argv[2]);    // output
    argv[1][strlen(argv[1])-10] = 0;     //removes "params.ini" in string "path/to/params.ini"
    const char* pth = argv[1];            // by setting "p" to 0 (end of string)

    // read files. The "file" types are vectors of structs
    prot_file prot;
    TTS_file tts;
    TSS_file tss;
    GFF_file gff_df_raw;
    readProt(prot, pth + params->BARR_FIX);
    readTTS(tts, pth + params->TTS); // TTS means "Transcription Termination Site"
    readTSS(tss, pth + params->TSS); // TSS means "Transcription Start Site"
    readGFF(gff_df_raw, pth + params->GFF);

    uint genome_size = get_genome_size(gff_df_raw);
    int genome = int(genome_size/params->DELTA_X);
    uint Niter = params->ITERATIONS_NB/params->DELTA_T;

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

    // a bit of control on what happens
    cout << endl << "Starting with individuals with the following genome:" << endl;
    for (Transcript t : tr) 
    {
        cout << "TU n°" << t.TUindex_ <<  ": ";
        cout << t.TSS_ << '-' << t.TTS_ << " (length = " << t.size_ << "), strand ";
        cout << t.s_ << " with rate " << t.r_ << endl;
    }
    cout << endl;

    vector<uint> ts_remain_all;
    transform(tr.begin(), tr.end(), back_inserter(ts_remain_all), 
                    [](Transcript const& x) { return x.size_n_; });


    // ====================== Topological barriers ============================
    vector<uint> Barr_fix; // Get the fixed topo barriers in a vector
    transform(prot.begin(), prot.end(), back_inserter(Barr_fix), 
                    [params](prot_t const& x) { return int(x.prot_pos/params->DELTA_X); });    
    vector<uint> Barr_pos (Barr_fix.size(), 0);
    copy(Barr_fix.begin(), Barr_fix.end(), Barr_pos.begin());
    vector<int> Dom_size (Barr_pos.size(), 0);
    adjacent_difference(Barr_pos.begin(), Barr_pos.end(), Dom_size.begin());
    Dom_size.erase(Dom_size.begin()); // because unlike np.ediff1d, adjacent_diff keeps the first element.
    Dom_size.push_back(Barr_fix[0]+(genome-Barr_fix.back()));
    vector<int> Barr_type (Barr_fix.size(), 0);
    vector<double> Barr_sigma (Barr_fix.size(), params->SIGMA_0);
    // For later efficiency while inserting RNAPs in the barriers:
    Barr_pos.reserve(Barr_fix.size()+params->RNAPS_NB);
    Barr_type.reserve(Barr_fix.size()+params->RNAPS_NB);
    Barr_sigma.reserve(Barr_fix.size()+params->RNAPS_NB);
    Dom_size.reserve(Barr_fix.size()+params->RNAPS_NB);

    // ===================== Polymerases initialization ======================

    uint N_RNAPs_unhooked(params->RNAPS_NB);    // Available polymerases unhooked.
    vector<RNAP> RNAPs_hooked;                    // Vector of hooked polymerases, empty for now
    RNAPs_hooked.reserve(params->RNAPS_NB);

    // ===================== Simulation ======================================
    vector<uint> TSS_pos_idx, rm_RNAPs_idx;
    vector<int> tss_and_unhooked_RNAPs, picked_tr;
    vector<double> init_rates, prob_init_rates, prob_unhooked_rates, all_prob;
    valarray<bool> isfinished;
    double sum_init_rates, prob_unhooked_rate;
    size_t index, j;
    DNApos pos;
    init_rates.reserve(tr.size());
    prob_init_rates.reserve(tr.size());
    all_prob.resize(tr.size() + N_RNAPs_unhooked);
    picked_tr.reserve(N_RNAPs_unhooked);
    rm_RNAPs_idx.reserve(N_RNAPs_unhooked);
    for (size_t t=0; t<Niter; ++t)
    {
        // cout << "=========== t = " << t << "===========" << endl;
        TSS_pos_idx.clear();
        init_rates.clear();
        prob_init_rates.clear(); 
        // get init_rate (raw)
        for (auto it = tr.begin(); it!=tr.end(); it++)
            init_rates.push_back(it->f_init_rate(params->sigma_t, params->epsilon, params->m, Barr_pos, Barr_sigma));
        // compute the probabilities of expression of each tr
        sum_init_rates = vector_sum(init_rates);
        for (auto it = tr.begin(); it!=tr.end(); it++) 
            prob_init_rates.push_back(it->f_prob_init_rate(sum_init_rates, params->DELTA_T));
        // cout << params->RNAPS_NB - N_RNAPs_unhooked << " RNAPs working. ";
        // if there are some unhooked RNAPs, hook them
        if (N_RNAPs_unhooked) 
        {
            //get the probability for a RNAP to stay unhooked
            prob_unhooked_rate = f_prob_unhooked_rate(sum_init_rates, params->DELTA_T, N_RNAPs_unhooked);
            prob_unhooked_rates.assign(N_RNAPs_unhooked, prob_unhooked_rate);
            // create an array that will contains [ 1 2 ... nTR , -1 -1 -1 ... ]
            tss_and_unhooked_RNAPs.assign(tr.size() + N_RNAPs_unhooked, -1);
            iota(tss_and_unhooked_RNAPs.begin(), tss_and_unhooked_RNAPs.begin() + tr.size(), 0);
            // concatenation of the probabilities arrays
            copy(prob_init_rates.begin(), prob_init_rates.end(), all_prob.begin());
            copy(prob_unhooked_rates.begin(), prob_unhooked_rates.end(), all_prob.begin()+tr.size());
            // pick up transcripts randomly
            random_choice(picked_tr, tss_and_unhooked_RNAPs, N_RNAPs_unhooked, all_prob);
            // cout << "tir des nouveaux tr: ";
            // if all RNAPs are hooked, assign transcripts to RNAPs
            for (int i : picked_tr)
            {
                // cout << i << ' ';
                if (i != -1)
                {
                    // cout << "(brin" << tr[i].s_ << ") ";
                    // Add a hooked polymerase in the vector, working with the chosen transcript.
                    RNAPs_hooked.push_back(RNAP(tr[i]));
                    N_RNAPs_unhooked--;
                    // Add the RNAP in the list of topological barriers:
                    pos = tr[i].start_;
                    index = find_if(Barr_pos.begin(), Barr_pos.end(), 
                            [pos](DNApos barr)->bool{ return barr>pos; }) - Barr_pos.begin();
                    Barr_pos.insert(Barr_pos.begin()+index, pos);
                    Barr_type.insert(Barr_type.begin()+index, tr[i].s_);
                    Barr_sigma.insert(Barr_sigma.begin()+index, Barr_sigma[index-1]);
                }
            }
            // cout << endl;
            // Update the topological domains 
            Dom_size.resize(Barr_pos.size());
            adjacent_difference(Barr_pos.begin(), Barr_pos.end(), Dom_size.begin());
            Dom_size.erase(Dom_size.begin()); // because unlike np.ediff1d, adjacent_diff keeps the first element.
            Dom_size.push_back(Barr_pos[0]+(genome-Barr_pos.back()));
        }
        // Move each polymerase on its transcript
        for (auto it=RNAPs_hooked.begin(); it!=RNAPs_hooked.end(); it++)
            it->move();
        for (uint i=0; i<Barr_pos.size(); i++)
            Barr_pos[i] += Barr_type[i];
        // cout << "progress: ";
        // for (RNAP pol : RNAPs_hooked)
        //     cout << abs(int(pol.last_pos_) - int(pol.pos_)) << " ";
        // cout << endl;
        // cout << "position of barriers: ";
        // for (DNApos pos : Barr_pos)
        //     cout << pos << " ";
        // cout << endl;
        // Check if RNAPs have finished a transcript
        rm_RNAPs_idx.clear();
        for (uint i = 0; i<RNAPs_hooked.size(); i++)
        {
            if (RNAPs_hooked[i].hasfinished())
            {
                // cout << "a pol has finished !" << endl;
                tr[RNAPs_hooked[i].tr_id_].expr_count_++;
                rm_RNAPs_idx.push_back(i);
            }
        }
        // Unhook the polymerase(s)
        j = 0;
        for (auto it = rm_RNAPs_idx.begin(); it != rm_RNAPs_idx.end(); it++)
        {
            pos = RNAPs_hooked[(*it) - j].pos_;
            index = find_if(Barr_pos.begin(), Barr_pos.end(), 
                            [pos](DNApos barr)->bool{ return barr>pos; }) - Barr_pos.begin();
            Barr_sigma[index-2] = (Dom_size[index-2]*Barr_sigma[index-2] + Dom_size[index-1]*Barr_sigma[index-1]) / (Dom_size[index-2] + Dom_size[index-1]);
            Barr_pos.erase(Barr_pos.begin() + index - 1);
            Barr_type.erase(Barr_type.begin() + index - 1);
            Barr_sigma.erase(Barr_sigma.begin() + index - 1);
            Dom_size.erase(Dom_size.begin() + index - 1);
            RNAPs_hooked.erase(RNAPs_hooked.begin() + (*it) - j);
            N_RNAPs_unhooked++;
            j++;
        }
        // if (rm_RNAPs_idx.size())
        // {
        //     cout << "removed " << rm_RNAPs_idx.size() << " polymerases." << endl;
        //     cout << "new barriers: ";
        //     for (DNApos pos : Barr_pos)
        //         cout << pos << " ";
        //     cout << endl;
        //     cout << "hooked pols: " << RNAPs_hooked.size() << endl;
        // }
        // Update the Dom_size:
        adjacent_difference(Barr_pos.begin(), Barr_pos.end(), Dom_size.begin());
        Dom_size.erase(Dom_size.begin());
        Dom_size.push_back(Barr_pos[0]+(genome-Barr_pos.back()));
        // Update Sigma:
        int a,b,d;
        for (uint i = 0; i<Barr_pos.size(); i++)
        {
            a = Barr_type[i];
            if (i+1==Barr_pos.size()) b = Barr_type[0];
            else b = Barr_type[i+1];
            d = Dom_size[i];
            switch (a)
            {
                case 0:
                    switch (b) 
                    {
                        case 1: Barr_sigma[i] = Barr_sigma[i] * (d-1.0)/d - RNAPs_genSC/d; break;
                        case -1: Barr_sigma[i] = Barr_sigma[i] * (d+1.0)/d + RNAPs_genSC/d; 
                    }
                    break;
                case -1:
                    switch (b) 
                    {
                        case 1: Barr_sigma[i] = Barr_sigma[i] * (d-2.0)/d - 2*RNAPs_genSC/d; break;
                        case 0: Barr_sigma[i] = Barr_sigma[i] * (d-1.0)/d - RNAPs_genSC/d; 
                    }
                    break;
                case 1:
                    switch (b) 
                    {
                        case -1: Barr_sigma[i] = Barr_sigma[i] * (d+2.0)/d + 2*RNAPs_genSC/d; break;
                        case 0: Barr_sigma[i] = Barr_sigma[i] * (d+1.0)/d + RNAPs_genSC/d; 
                    }
            }
        }
        calc_sigma(Barr_sigma, params);
        // cout << "torsion: ";
        // display_vector(Barr_sigma);
    }

    cout << "Simulation completed successfully !! \nNumber of transcripts : "<< endl;
    for (uint i=0; i<tr.size(); i++)
        cout << "Transcript " << i << " : " << tr[i].expr_count_ << endl; 

    // delete loaded param files
    delete params;
    return(EXIT_SUCCESS);
}
