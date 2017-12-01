#include "utils.h"

Params *readIni(const char *cfgFile){
    Params* conf = new Params;
    parseIniFile(cfgFile);
	conf->GFF = getOptionToString("GFF");
	conf->TSS = getOptionToString("TSS");
	conf->TTS = getOptionToString("TTS");
	conf->BARR_FIX = getOptionToString("BARR_FIX");
	conf->J_0 = getOptionToDouble("J_0");
	conf->D = getOptionToInt("D");
	conf->m = getOptionToDouble("m");
	conf->sigma_t = getOptionToDouble("sigma_t");
	conf->SIGMA_0 = getOptionToDouble("SIGMA_0");
	conf->DELTA_X = getOptionToInt("DELTA_X");
	conf->DELTA_T = getOptionToInt("DELTA_Y");
	conf->RNAPS_NB = getOptionToInt("RNAPS_NB");
	conf->ITERATIONS_NB = getOptionToInt("ITERATIONS_NB");
	conf->OUTPUT_STEP = getOptionToInt("OUTPUT_STEP");
	conf->GYRASE_CONC = getOptionToDouble("GYRASE_CONC");
	conf->TOPO_CONC = getOptionToDouble("TOPO_CONC");
	conf->GYRASE_EFFICIENCY = getOptionToDouble("GYRASE_EFFICIENCY");
	conf->TOPO_EFFICIENCY = getOptionToDouble("TOPO_EFFICIENCY");
	conf->GYRASE_CTE = getOptionToDouble("GYRASE_CTE");
	conf->TOPO_CTE = getOptionToDouble("TOPO_CTE");
	conf->k_GYRASE = getOptionToInt("k_GYRASE");
	conf->x0_GYRASE = getOptionToDouble("x0_GYRASE");
	conf->k_TOPO = getOptionToInt("k_TOPO");
	conf->x0_TOPO = getOptionToDouble("x0_TOPO");
	cleanupIniReader();
	cout << cfgFile << " parsed successfully." << endl; //Should return nothing as the config items have been cleaned

    return conf;
}

void readProt(prot_file& data, string protFile){
	ifstream file(protFile.c_str());
	string str; 
	getline(file, str); // headers
    while (getline(file, str))
    {
        vector<string> v;
		boost::split(v, str, ::isspace);
		if (v.size() == 2) {
			prot_t line = { v[0], static_cast<uint>(atoi(v[1].c_str())) };
			data.push_back( line );
		}
    }
	display_vector(data);
	cout << endl << protFile << " parsed successfully." << endl;
}

void readTSS(TSS_file& data, string TSSFile){
	ifstream file(TSSFile.c_str());
	string str; 
	getline(file, str); // headers
    while (getline(file, str))
    {
        vector<string> v;
		boost::split(v, str, ::isspace);
		if (v.size() == 4) {
			TSS_t line = {  static_cast<uint>(atoi(v[0].c_str())), 
							*(v[1].c_str()), 
							static_cast<uint>(atoi(v[2].c_str())), 
							atof(v[3].c_str())
						};
			data.push_back( line );
		}
	}
    // display_vector_star(data);
	cout << TSSFile << " parsed successfully." << endl;
}

void readTTS(TTS_file& data, string TTSFile){
	ifstream file(TTSFile.c_str());
	string str; 
	getline(file, str); // headers
    while (getline(file, str))
    {
        vector<string> v;
		boost::split(v, str, ::isspace);
		if (v.size() == 4) {
			TTS_t line = { 	static_cast<uint>(atoi(v[0].c_str())), 
							*(v[1].c_str()), 
							static_cast<uint>(atoi(v[2].c_str())), 
							atof(v[3].c_str())
						};
			data.push_back( line );
		}
    }
	cout << TTSFile << " parsed successfully." << endl;
}

void readGFF(GFF_file& data, string GFFFile){
	ifstream file(GFFFile.c_str());
	string str; 
	do { 
		getline(file, str); 
	} while (str[0] == '#');
    do {
        vector<string> v;
		boost::split(v, str, ::isspace);
		if (v.size() == 9) {
			GFF_t line = { 	v[0],
							v[1],
							v[2],
							static_cast<uint>(atoi(v[3].c_str())),
							static_cast<uint>(atoi(v[4].c_str())),
							atof(v[5].c_str()),
							*(v[6].c_str()),
							atoi(v[7].c_str()),
							v[8]
						};
			data.push_back( line );
		}
    } while (getline(file, str));
	cout << GFFFile << " parsed successfully." << endl;
}

ostream &operator<<(ostream &stream, TSS_t const &s) { 
    return stream << "{ " << s.TUindex << " " << s.TUorient << " " << s.TSS_pos << " " << s.TSS_strength << " }";
}

ostream &operator<<(ostream &stream, prot_t const &s) { 
    return stream << "{ " << s.prot_name << " " << s.prot_pos << " }";
}

ostream &operator<<(ostream &stream, TTS_t const &s) { 
    return stream << "{ " << s.TUindex << " " << s.TUorient << " " << s.TTS_pos << " " << s.TTS_proba_off << " }";
}

ostream &operator<<(ostream &stream, GFF_t const &s) { 
    return stream 	<< "{ " 
					<< s.seqname << " " 
					<< s.source << " " 
					<< s.feature << " " 
					<< s.start << " " 
					<< s.end << " " 
					<< s.score << " " 
					<< s.strand << " " 
					<< s.frame << " " 
					<< s.attribute << " }";
}

template<typename file_type>
void    display_vector(file_type& v){
    for (size_t n = 0; n < v.size(); n++)
    	cout << v[n] << endl;
  	cout << endl;
}

uint	get_genome_size(GFF_file& gff_df) {
	// This is dirty. Guess genome size from the first annotation in GFF.
	GFF_t full_genome = gff_df[0];
	return full_genome.end - full_genome.start + 1;
}

map< uint , vector<uint> > get_TU_tts(TSS_file& tss, TTS_file& tts) {
	vector<uint> TU_values;
	std::transform(tss.begin(), tss.end(), std::back_inserter(TU_values), 
					[](TSS_t const& x) { return x.TUindex; });
	vector<uint> TTS_pos;
	std::transform(tts.begin(), tts.end(), std::back_inserter(TTS_pos), 
					[](TTS_t const& x) { return x.TTS_pos; });

	map< uint , vector<uint> > TU_tts;
	for ( 	size_t i = 0, TU_index_val = TU_values[0] ; 
			i < TU_values.size() ;
			i++, TU_index_val = TU_values[i]
		) 
	{
		TU_tts[TU_index_val].push_back(TTS_pos[i]);
	}

	// display
	// for (size_t i = 0; i < TU_values.size(); i++) 
	// {
	// 	cout << "TU n°" << TU_values[i] << ", tts at ";
	// 	display_vector(TU_tts[TU_values[i]]);
	// }

	return TU_tts;
}