#include <cstdlib>
#include "utils.h"
#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>

using namespace std;

typedef boost::multi_array<int, 3> array3;
typedef array3::index arr_index;

template<typename file_type>
void    display_vector(file_type& v){
    for (size_t n = 0; n < v.size(); n++)
    	cout << v[n] << " ";
  	cout << endl;
}

template<typename T>
T 	vector_sum(vector<T>& V) {
	T sum = 0; for (size_t i = 0; i < V.size(); i++) sum += V[i];
	return sum;
}

int main(int argc, char** argv) {
	if (!argc) cout << "torsionEvol path/to/params.ini";
	Params* params = readIni(argv[1]);

	// define the input/output directories
	boost::filesystem::create_directories(argv[2]);	// output
	argv[1][strlen(argv[1])-10] = 0; 	//removes "params.ini" in string "path/to/params.ini"
	const char* pth = argv[1];			// by setting "p" to 0 (end of string)

	// read files. The "file" types are vectors of structs
	prot_file prot;
	TTS_file tts;
	TSS_file tss;
	GFF_file gff_df_raw;
	readProt(prot, pth + params->BARR_FIX);
	readTTS(tts, pth + params->TTS); // TTS means "Transcription Termination Site"
	readTSS(tss, pth + params->TSS); // TSS means "Transcription Start Site"
	readGFF(gff_df_raw, pth + params->GFF);


	// extract TSS_pos in a vector (transcriptions starting sites)
	vector<uint> TSS_pos;
	std::transform(tss.begin(), tss.end(), std::back_inserter(TSS_pos), 
					[params](TSS_t const& x) { return x.TSS_pos/params->DELTA_X; });
	
	// extract Kon in a vector (probability that a polymerase binds to the promoter ?)
	vector<double> Kon;
	std::transform(tss.begin(), tss.end(), std::back_inserter(Kon), 
					[](TSS_t const& x) { return x.TSS_strength; });				

	// extract Poff in a vector (probability that the polymerase drops at the end of transcript)
	vector<double> Poff;
	std::transform(tts.begin(), tts.end(), std::back_inserter(Poff), 
					[](TTS_t const& x) { return x.TTS_proba_off; });	

	// Genome size
	uint genome_size = get_genome_size(gff_df_raw);

	// Map of transciption units with the list of tts belonging to TU. [ TU ]  = (tts1, tts2, ... )
	map< uint , vector<uint> > TU_tts = get_TU_tts(tss, tts);

	// The RNAPs id
	vector<uint> RNAPs_id (params->RNAPS_NB);
	for (int i = 0; i < params->RNAPS_NB; i++) { RNAPs_id[i] = i; }

	// The position of RNAPs (np.zeros)
	vector<uint> RNAPs_pos (params->RNAPS_NB, 0);

	// RNAPs_last_pos (np.zeros)
	vector<uint> RNAPs_last_pos (params->RNAPS_NB, 0);

	// Strands orientation
	vector<int> strands;
	std::transform(gff_df_raw.begin(), gff_df_raw.end(), std::back_inserter(strands), 
					[](GFF_t const& x) { return 1*(x.strand == '+') - 1*(x.strand == '-'); });	

	// Get TR info
	vector<uint> this_TU_tts;
    vector<uint> tr_id;
    vector<uint> tr_start;
    vector<uint> tr_end;
    vector<int> tr_strand;
	vector<double> tr_rate;
	double sum_Kon = vector_sum(Kon);
	int j = 0; // trancript id indice

	// for (TSS_t line: *tss) {
	for (size_t i = 0; i<tss.size(); i++) {
		TSS_t* transcript_start = &(tss[i]);  //pointer on the ith TSS_t (we should do this with an iterator instead)
	    
        uint TU_id = transcript_start->TUindex; // get the TU of this tss

        // Get the list of TTS that are in the same TU of this tss_id (i)
        // TU_tts ex : { <0: {1150, 2350}>, <1: {6250}> }  ==>  0 -> [1150, 2350]
        vector<uint> this_TU_tts = TU_tts[TU_id];


        if (transcript_start->TUorient == '+')  // go right
		{
        	// tr_rate ======>  [ 0.1875   0.03125  0.00625  0.025    0.0625   0.25     0.125    0.3125 ]
            uint k = TU_id; // TTS id index : k start from the first position of each TU
            float proba_rest = 1.0;
            while (proba_rest > 0 ) {
				TTS_t* transcript_stop = &(tts[k]); //pointer on the kth TTS_t (we should do this with an iterator)
                if (transcript_start->TSS_pos < transcript_stop->TTS_pos) {
                    tr_id.push_back(j);
                    tr_strand.push_back(1);
                    tr_start.push_back(transcript_start->TSS_pos);
                    tr_end.push_back(transcript_stop->TTS_pos);
                    // the probability to choose a specific transcript
                    tr_rate.push_back(Kon[i] * (Poff[k] * proba_rest));
                    proba_rest = (1 - Poff[k]) * proba_rest;
					j++;
				}
				k++;
			}
		} 
		else // go left
		{
            uint k = 0; //TU_id #0 #len(this_TU_tts)# TTS id index ### [0 1 2 3 4 5]
            float proba_rest = 1.0;
            while (proba_rest>0 and k<this_TU_tts.size()) { // >= 0 : #
				TTS_t* transcript_stop = &(tts[k]);  //pointer on the kth TTS_t (we should do this with an iterator)	
                if (transcript_stop->TTS_pos < transcript_start->TSS_pos) {
                    tr_id.push_back(j);
                    tr_strand.push_back(-1);
                    tr_start.push_back(transcript_start->TSS_pos);
                    tr_end.push_back(this_TU_tts[k]);
                    // the probability to choose a specific transcript
                    tr_rate.push_back(Kon[TU_id] * (Poff[k] * proba_rest));
                    proba_rest = (1 - Poff[k]) * proba_rest;
					j++;
				}
				k++;
			}
		}
	}
	// compute the size of transcripts
	vector<uint> tr_size (tr_id.size());
	std::transform(tr_end.begin(), tr_end.end(), tr_start.begin(), tr_size.begin(), minus<int>()); // size = (end - start)
	std::transform(tr_size.begin(), tr_size.end(), tr_size.begin(), ::abs);  // size =  |size| 
	vector<uint> ts_beg_all_trs (tr_id.size(), 0);
	vector<uint> ts_remain_all (tr_id.size());
	std::copy(tr_size.begin(), tr_size.end(), ts_remain_all.begin());
	
	//Divide positions by DELTA_X (yes these are ints. yes this is dirty.)
	cout << "scaling..." << endl;
	std::transform(tr_start.begin(), tr_start.end(), tr_start.begin(), bind2nd(divides<int>(), params->DELTA_X));
	std::transform(tr_end.begin(), tr_end.end(), tr_end.begin(), bind2nd(divides<int>(), params->DELTA_X));
	std::transform(tr_size.begin(), tr_size.end(), tr_size.begin(), bind2nd(divides<int>(), params->DELTA_X));
	std::transform(ts_remain_all.begin(), ts_remain_all.end(), ts_remain_all.begin(), bind2nd(divides<int>(), params->DELTA_X));

	// a bit of control on what happens
	for (size_t i = 0; i < tr_size.size(); i++) 
	{
		cout << "TU n°" << tr_id[i] <<  ": ";
		cout << tr_start[i] << '-' << tr_end[i] << " (length = " << tr_size[i] << "), strand ";
		cout << tr_strand[i] << " with rate " << tr_rate[i] << endl;
	}


	// The number of times transcripts has been transcribed
	vector<uint> tr_nbr (tr_id.size(), 0);
	
	int genome = int(genome_size/params->DELTA_X);

	// Get the topo barriers
	vector<uint> Barr_fix;
	std::transform(prot.begin(), prot.end(), std::back_inserter(Barr_fix), 
					[params](prot_t const& x) { return int(x.prot_pos/params->DELTA_X); });	

	// just for the echo we can assign it directely
    vector<uint> Barr_pos (Barr_fix.size(), 0);
	std::copy(Barr_fix.begin(), Barr_fix.end(), Barr_pos.begin());
	vector<int> Dom_size (Barr_pos.size(), 0);
	std::adjacent_difference(Barr_pos.begin(), Barr_pos.end(), Dom_size.begin());
	Dom_size.erase(Dom_size.begin()); // because unlike np.ediff1d, std::adjacent_diff keeps the first element.
	Dom_size.push_back(genome-Barr_fix[-1]+Barr_fix[0]); // !! change Barr_fix to Barr_pos case : O | |	

	vector<int> Barr_type (Barr_fix.size(), 0);
	vector<double> Barr_sigma (Barr_fix.size(), params->SIGMA_0);

	// here we need to make an Barr_ts_remain
	// to track the position of each RNAPol
    // each position in Barr_ts_remain is associated with the same position in Barr_pos
    vector<uint> Barr_ts_remain (Barr_fix.size()); 
	// The Barr_ts_remain of fixed barr is NaN, but in C++ we did not initialize the values !

	// ######### Variables used to get the coverage ##########

	vector<uint> id_shift_fwd(genome-1);
	iota(id_shift_fwd.begin(), id_shift_fwd.end(), 1); //range(1,genome)
    id_shift_fwd.push_back(0);

	vector<uint> id_shift_bwd(genome-1);
	iota(id_shift_bwd.begin(), id_shift_bwd.end(), 0); //range(1,genome)
    id_shift_bwd.insert(id_shift_bwd.begin(), genome-1);

	vector<uint> cov_bp(genome);
	uint val = 0;
	auto first = cov_bp.begin();
	while (first != cov_bp.end()) {
		*first = val;
		++first;
		val += params->DELTA_X;
	}

	// cout << "======== vectors: ==========" << endl;
	// cout << "Barr_fix: " << endl;
	// display_vector(Barr_fix);
	// cout << "Barr_pos: " << endl;
	// display_vector(Barr_pos);
	// cout << "Barr_type: " << endl;
	// display_vector(Barr_type);
	// cout << "Barr_sigma: " << endl;
	// display_vector(Barr_sigma);
	// cout << "Barr_ts_remain: " << endl;
	// display_vector(Barr_ts_remain);
	// cout << "Dom_size: " << endl;
	// display_vector(Dom_size);
	// cout << "id_shift_fwd: " << endl;
	// display_vector(id_shift_fwd);
	// cout << "id_shift_bwd: " << endl;
	// display_vector(id_shift_bwd);
	// cout << "genome: " << endl;
	// cout << genome << ' ' << genome_size << " " << params->DELTA_X << endl;


	//###########################################################
    //#                 initiation of values                    #
    //###########################################################

	// save the time when RNApoly is starting trasncribing a specific transcript
	map<uint, vector<int> > tr_times;

	// array where will save all RNAPs info (Create a 3D array)
	uint Niter = params->ITERATIONS_NB/params->DELTA_T;
	array3 save_RNAPs_infos(boost::extents[params->RNAPS_NB][2][Niter]);
	for(arr_index i = 0; i != params->RNAPS_NB; ++i) 
		for(arr_index j = 0; j != 2; ++j)
			for(arr_index k = 0; k != Niter; ++k)
				save_RNAPs_infos[i][j][k] = -1; 

	// the same for transcript info
	array3 save_tr_info (boost::extents[tr_id.size()][2][Niter]);
	for(arr_index i = 0; i != int(tr_id.size()); ++i) 
		for(arr_index j = 0; j != 2; ++j)
			for(arr_index k = 0; k != Niter; ++k)
				save_tr_info[i][j][k] = -1; 

	// in those variables, we will save/append info in each time step to save them as --> all_res ;-)
    vector<double> save_Barr_sigma;
    vector<double> save_Dom_size;
    vector<double> save_mean_sig_wholeGenome;

    // ########### Go !

	vector<uint> RNAPs_unhooked_id (RNAPs_id.size(), 0);
	std::copy(RNAPs_id.begin(), RNAPs_id.end(), RNAPs_unhooked_id.begin());

	vector<int> RNAPs_strand (params->RNAPS_NB, -1);
	vector<int> ts_beg (params->RNAPS_NB, -1);
	vector<int> ts_remain (params->RNAPS_NB, -1);
	// RNAPs_tr will contain the id of the picked transcript
    vector<int> RNAPs_tr (params->RNAPS_NB, -1);
	// get the TSSs ids
	vector<uint> tss_id;
	transform(tss.begin(), tss.end(), std::back_inserter(tss_id), 
				[params](TSS_t const& x) { return x.TUindex; });

	// in the case of RNAP_NBR = 0
    //RNAPs_hooked_id = []  (not pertinent in C++)

	// for (size_t t=0; t<Niter; ++t)
	for (size_t t=0; t<2; ++t)
	{
		// we need to know each TSS belong to which Domaine
        vector<uint> TSS_pos_idx;
        vector<double> init_rate, sigma_tr_start, prob_init_rate;
		searchsorted(TSS_pos_idx, Barr_pos, TSS_pos);
        // after knowing the domaine of each TSS we can get sigma
		for (auto i : TSS_pos_idx)
			sigma_tr_start.push_back(Barr_sigma[i-1]);
        // get the initiation rates
		f_init_rate(init_rate, tr_rate, sigma_tr_start, params->sigma_t, params->epsilon, params->m);
        double sum_init_rate = vector_sum(init_rate);

        f_prob_init_rate(prob_init_rate, init_rate, sum_init_rate, params->DELTA_T);
		cout << "prob init rate: " << endl;
		display_vector(prob_init_rate);

		if (RNAPs_unhooked_id.size())
		{
			//get the unhooked rates
            double prob_unhooked_rate = f_prob_unhooked_rate(sum_init_rate, params->DELTA_T, RNAPs_unhooked_id.size());
            vector<double> prob_unhooked_rates (RNAPs_unhooked_id.size(), prob_unhooked_rate);
			

			vector<double> all_prob; // concatenation
			all_prob.reserve( prob_init_rate.size() + prob_unhooked_rates.size() ); // preallocate memory
			all_prob.insert( all_prob.end(), prob_init_rate.begin(), prob_init_rate.end() );
			all_prob.insert( all_prob.end(), prob_unhooked_rates.begin(), prob_unhooked_rates.end() );
            // create the numpy array that will contains [ nTSS , Unhooked RNAPS ]
            vector<int> tss_and_unhooked_RNAPs (tss_id.size()+RNAPs_unhooked_id.size(), -1);
			std::copy(tss_id.begin(), tss_id.end(), tss_and_unhooked_RNAPs.begin());

            // pick up
            //picked_tr = np.random.choice(tss_and_unhooked_RNAPs, len(RNAPs_unhooked_id), replace=False, p=all_prob) // RNAPs_unhooked_id

            // # T# get the unhooked rates
            // prob_unhooked_rate = f_prob_unhooked_rate(sum_init_rate, DELTA_T, len(RNAPs_unhooked_id))
            
		}
	}




	// delete loaded param files
	delete params;
	return(EXIT_SUCCESS);
}
