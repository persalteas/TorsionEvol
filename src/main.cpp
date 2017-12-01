#include <cstdlib>
#include "utils.h"
#include <boost/filesystem.hpp>

using namespace std;

template<typename file_type>
void    display_vector(file_type& v){
    for (size_t n = 0; n < v.size(); n++)
    	cout << v[n] << endl;
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

	


	// delete loaded param files
	delete params;
	return(EXIT_SUCCESS);
}
