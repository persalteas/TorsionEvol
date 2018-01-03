#include <cstdlib>
#include <ctime>
#include <boost/filesystem.hpp>
#include "utils.h"
#include "Transcript.h"

using namespace std;

template<typename T>
T 	vector_sum(vector<T>& V) {
	T sum = 0; for (size_t i = 0; i < V.size(); i++) sum += V[i];
	return sum;
}

int main(int argc, char** argv) {
	if (!argc) cout << "torsionEvol path/to/params.ini";
	Params* params = readIni(argv[1]);
	srand(time(NULL));

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




	// Genome size
	uint genome_size = get_genome_size(gff_df_raw);
	// extract TSS_pos in a vector (transcriptions starting sites)
	vector<uint> TSS_pos;
	std::transform(tss.begin(), tss.end(), std::back_inserter(TSS_pos), 
					[params](TSS_t const& x) { return x.TSS_pos/params->DELTA_X; });
	// extract Kon in a vector (probability that a polymerase binds to the promoter = TSS strength)
	vector<double> Kon;
	std::transform(tss.begin(), tss.end(), std::back_inserter(Kon), 
					[](TSS_t const& x) { return x.TSS_strength; });				
	// extract Poff in a vector (probability that the polymerase drops at the TTS = TTS strength)
	vector<double> Poff;
	std::transform(tts.begin(), tts.end(), std::back_inserter(Poff), 
					[](TTS_t const& x) { return x.TTS_proba_off; });
	// Get the topo barriers
	vector<uint> Barr_fix;
	std::transform(prot.begin(), prot.end(), std::back_inserter(Barr_fix), 
					[params](prot_t const& x) { return int(x.prot_pos/params->DELTA_X); });	
	// Map of transciption units with the list of tts belonging to TU. [ TU ]  = (tts1, tts2, ... )
	map< uint , vector<uint> > TU_tts = get_TU_tts(tss);
	// Strands orientation
	vector<int> strands;
	std::transform(gff_df_raw.begin(), gff_df_raw.end(), std::back_inserter(strands), 
					[](GFF_t const& x) { return 1*(x.strand == '+') - 1*(x.strand == '-'); });	

	// Get possible transcripts
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
	for (Transcript t : tr) 
	{
		cout << "TU n°" << t.TUindex_ <<  ": ";
		cout << t.start_ << '-' << t.end_ << " (length = " << t.size_n_ << "*dx), strand ";
		cout << t.s_ << " with rate " << t.r_ << endl;
	}

	// vector<uint> ts_beg_all_trs (tr.size(), 0);
	vector<uint> ts_remain_all;
	std::transform(tr.begin(), tr.end(), std::back_inserter(ts_remain_all), 
					[](Transcript const& x) { return x.size_n_; });

	// The number of times transcripts has been transcribed
	vector<uint> tr_nbr (tr.size(), 0);
	int genome = int(genome_size/params->DELTA_X);

	// Topological barriers
    vector<uint> Barr_pos (Barr_fix.size(), 0);
	std::copy(Barr_fix.begin(), Barr_fix.end(), Barr_pos.begin());
	vector<int> Dom_size (Barr_pos.size(), 0);
	std::adjacent_difference(Barr_pos.begin(), Barr_pos.end(), Dom_size.begin());
	Dom_size.erase(Dom_size.begin()); // because unlike np.ediff1d, std::adjacent_diff keeps the first element.
	Dom_size.push_back(Barr_fix[0]+(genome-Barr_fix.back())); // !! change Barr_fix to Barr_pos case : O | |	
	vector<int> Barr_type (Barr_fix.size(), 0);
	vector<double> Barr_sigma (Barr_fix.size(), params->SIGMA_0);

	// here we need to make an Barr_ts_remain
	// to track the position of each RNAPol
    // each position in Barr_ts_remain is associated with the same position in Barr_pos
    vector<uint> Barr_ts_remain (Barr_fix.size()); 
	// The Barr_ts_remain of fixed barr is NaN, but in C++ we did not initialize the values !

	// For later efficiency while inserting RNAPs in the barriers:
	Barr_pos.reserve(Barr_fix.size()+params->RNAPS_NB);
	Barr_type.reserve(Barr_fix.size()+params->RNAPS_NB);
	Barr_sigma.reserve(Barr_fix.size()+params->RNAPS_NB);
	Barr_ts_remain.reserve(Barr_fix.size()+params->RNAPS_NB);
	Dom_size.reserve(Barr_fix.size()+params->RNAPS_NB);

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

	//###########################################################
    //#                 initiation of values                    #
    //###########################################################

	// save the time when RNApoly is starting trasncribing a specific transcript
	// map<uint, vector<int> > tr_times;

	// array where will save all RNAPs info (Create a 3D array)
	uint Niter = params->ITERATIONS_NB/params->DELTA_T;
	// array3 save_RNAPs_infos(boost::extents[params->RNAPS_NB][2][Niter]);
	// for(arr_index i = 0; i != params->RNAPS_NB; ++i) 
	// 	for(arr_index j = 0; j != 2; ++j)
	// 		for(arr_index k = 0; k != Niter; ++k)
	// 			save_RNAPs_infos[i][j][k] = -1; 

	// // the same for transcript info
	// array3 save_tr_info (boost::extents[tr_id.size()][2][Niter]);
	// for(arr_index i = 0; i != int(tr_id.size()); ++i) 
	// 	for(arr_index j = 0; j != 2; ++j)
	// 		for(arr_index k = 0; k != Niter; ++k)
	// 			save_tr_info[i][j][k] = -1; 

	// // in those variables, we will save/append info in each time step to save them as --> all_res ;-)
    // vector<double> save_Barr_sigma;
    // vector<double> save_Dom_size;
    // vector<double> save_mean_sig_wholeGenome;

    // ########### Go !
	cout << endl;
	uint N_RNAPs_unhooked = params->RNAPS_NB;
	vector<int> RNAPs_pos (params->RNAPS_NB, -1);
	vector<int> RNAPs_last_pos (params->RNAPS_NB, -1);
    vector<int> RNAPs_tr (params->RNAPS_NB, -1); // RNAPs_tr will contain the id of the picked transcript
	vector<int> RNAPs_strand (params->RNAPS_NB, 0);
	vector<int> ts_beg (params->RNAPS_NB, -1);
	vector<int> ts_remain (params->RNAPS_NB, -1);

	Niter = 2; // debug
	for (size_t t=0; t<Niter; ++t)
	{
		cout << "-------------------------" << endl;
		vector<uint> TSS_pos_idx;
		vector<double> init_rates, prob_init_rates;
		for (auto it = tr.begin(); it!=tr.end(); it++)
		{
			// update init_rate
			it->f_init_rate(params->sigma_t, params->epsilon, params->m, Barr_pos, Barr_sigma);
			init_rates.push_back(it->init_rate_);
		}
		double sum_init_rates = vector_sum(init_rates);
		for (auto it = tr.begin(); it!=tr.end(); it++) 
		{
			// compute the expression probabilities
			it->f_prob_init_rate(sum_init_rates, params->DELTA_T);
			prob_init_rates.push_back(it->prob_init_rate_);
		}
		
		if (N_RNAPs_unhooked) // There still are some unfixed RNAPs
		{
			//get the probability for a RNAP to remain unhooked
			double prob_unhooked_rate = f_prob_unhooked_rate(sum_init_rates, params->DELTA_T, N_RNAPs_unhooked);
			vector<double> prob_unhooked_rates(N_RNAPs_unhooked, prob_unhooked_rate);
			// create an array that will contains [ 1 2 ... nTSS , -1 -1 -1 ... ]
			vector<int> tss_and_unhooked_RNAPs (tr.size() + N_RNAPs_unhooked, -1);
			iota(tss_and_unhooked_RNAPs.begin(), tss_and_unhooked_RNAPs.begin()+tr.size(), 0);
			// concatenation of the probabilities arrays
			vector<double> all_prob; 
			all_prob.reserve( prob_init_rates.size() + prob_unhooked_rates.size() ); // preallocate memory
			all_prob.insert( all_prob.end(), prob_init_rates.begin(), prob_init_rates.end() );
			all_prob.insert( all_prob.end(), prob_unhooked_rates.begin(), prob_unhooked_rates.end() );

			// pick up
			vector<int> picked_tr;
			random_choice(picked_tr, tss_and_unhooked_RNAPs, N_RNAPs_unhooked, all_prob);
			cout << "picked_tr: ";
			display_vector(picked_tr);
			cout << "ts_remain_all: ";
			display_vector(ts_remain_all);

			cout << endl << "placing RNAPs..." << endl << endl;
			cout << "RNAPs_tr: ";
			display_vector(RNAPs_tr);
			cout << "RNAPs_strand: ";
			display_vector(RNAPs_strand);
			cout << "RNAPs_pos: ";
			display_vector(RNAPs_pos);
			cout << "ts_remain: ";
			display_vector(ts_remain);
			cout << endl;
			
			if (std::find(picked_tr.begin(), picked_tr.end(), -1)==picked_tr.end()) // if all RNAPs are hooked
			{
				std::copy(picked_tr.begin(), picked_tr.end(), RNAPs_tr.begin());
				for (size_t i=0; i<picked_tr.size(); i++)
				{
					RNAPs_strand[i] = tr[picked_tr[i]].s_;
					RNAPs_pos[i] = tr[picked_tr[i]].start_;
					RNAPs_last_pos[i] = tr[picked_tr[i]].end_;
					ts_beg[i] = 0;
					ts_remain[i] = ts_remain_all[picked_tr[i]];
				}

			}
			cout << "RNAPs_tr: ";
			display_vector(RNAPs_tr);
			cout << "RNAPs_strand: ";
			display_vector(RNAPs_strand);
			cout << "RNAPs_pos: ";
			display_vector(RNAPs_pos);
			cout << "ts_remain: ";
			display_vector(ts_remain);

			// Add the RNAPs in the list of topological barriers:
			cout << endl << "Adding RNAPs to topological barriers..." << endl << endl;
			cout << "Barr_pos: ";
			display_vector(Barr_pos);
			cout << "Barr_type: ";
			display_vector(Barr_type);
			cout << "Dom_size: ";
			display_vector(Dom_size);
			cout << "Barr_sigma: ";
			display_vector(Barr_sigma);
			cout << endl;
			vector<uint> Barr_pos_RNAPs_idx;
			searchsorted<int>(Barr_pos_RNAPs_idx, Barr_pos, RNAPs_pos);
			size_t index;
			DNApos pos;
			for (size_t i=0; i<RNAPs_pos.size(); i++)
			{
				pos = RNAPs_pos[i];
				index = find_if(Barr_pos.begin(), Barr_pos.end(), 
						[pos](DNApos barr)->bool{ return barr>pos; }) - Barr_pos.begin();
				Barr_pos.insert(Barr_pos.begin()+index, pos);
				Barr_type.insert(Barr_type.begin()+index, RNAPs_strand[i]);
				Barr_sigma.insert(Barr_sigma.begin()+index, Barr_sigma[index-1]);
			}
			std::adjacent_difference(Barr_pos.begin(), Barr_pos.end(), Dom_size.begin());
			Dom_size.erase(Dom_size.begin()); // because unlike np.ediff1d, std::adjacent_diff keeps the first element.
			Dom_size.push_back(Barr_pos[0]+(genome-Barr_pos.back()));

			cout << "Barr_pos: ";
			display_vector(Barr_pos);
			cout << "Barr_type: ";
			display_vector(Barr_type);
			cout << "Dom_size: ";
			display_vector(Dom_size);
			cout << "Barr_sigma: ";
			display_vector(Barr_sigma);			
		}

		// Only to remove elements while this code is not finished 
		// and avoid seg faults:
		Barr_pos.resize(Barr_fix.size());
		Barr_type.resize(Barr_fix.size());
		Barr_sigma.resize(Barr_fix.size());
		Dom_size.resize(Barr_fix.size());
	}




	// delete loaded param files
	delete params;
	return(EXIT_SUCCESS);
}
