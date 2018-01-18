#include "utils.h"
#include "IniReader.h"
#include <algorithm>
#include <cmath>

Params *readIni(const char *cfgFile) {
  Params *conf = new Params;
  parseIniFile(cfgFile); // from IniReader.h
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
  conf->DELTA_T = getOptionToInt("DELTA_T");
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
  cout << cfgFile << " parsed successfully."
       << endl; // Should return nothing as the config items have been cleaned

  return conf;
}

void readProt(prot_file &data, string protFile) {
  ifstream file(protFile.c_str());
  string str;
  getline(file, str); // headers
  while (getline(file, str)) {
    vector<string> v;
    boost::split(v, str, boost::is_any_of("\t"));
    if (v.size() == 2) {
      prot_t line = {v[0], static_cast<uint>(atoi(v[1].c_str()))};
      data.push_back(line);
    }
  }
  cout << protFile << " parsed successfully." << endl;
}

void readTSS(TSS_file &data, string TSSFile) {
  ifstream file(TSSFile.c_str());
  string str;
  getline(file, str); // headers
  while (getline(file, str)) {
    vector<string> v;
    boost::split(v, str, ::isspace);
    if (v.size() == 4) {
      TSS_t line = {static_cast<uint>(atoi(v[0].c_str())), *(v[1].c_str()),
                    static_cast<uint>(atoi(v[2].c_str())), atof(v[3].c_str())};
      data.push_back(line);
    }
  }
  // display_vector_star(data);
  cout << TSSFile << " parsed successfully." << endl;
}

void readTTS(TTS_file &data, string TTSFile) {
  ifstream file(TTSFile.c_str());
  string str;
  getline(file, str); // headers
  while (getline(file, str)) {
    vector<string> v;
    boost::split(v, str, ::isspace);
    if (v.size() == 4) {
      TTS_t line = {static_cast<uint>(atoi(v[0].c_str())), *(v[1].c_str()),
                    static_cast<uint>(atoi(v[2].c_str())), atof(v[3].c_str())};
      data.push_back(line);
    }
  }
  cout << TTSFile << " parsed successfully." << endl;
}

void readGFF(GFF_file &data, string GFFFile) {
  ifstream file(GFFFile.c_str());
  string str;
  do {
    getline(file, str);
  } while (str[0] == '#');
  do {
    vector<string> v;
    boost::split(v, str, ::isspace);
    if (v.size() == 9) {
      GFF_t line = {v[0],
                    v[1],
                    v[2],
                    static_cast<uint>(atoi(v[3].c_str())),
                    static_cast<uint>(atoi(v[4].c_str())),
                    atof(v[5].c_str()),
                    *(v[6].c_str()),
                    atoi(v[7].c_str()),
                    v[8]};
      data.push_back(line);
    }
  } while (getline(file, str));
  cout << GFFFile << " parsed successfully." << endl;
}

void readEnv(vector<double> &env, const char *Envfile) {
  ifstream file(Envfile);
  string str;
  env.clear();
  do {
    vector<string> v;
    boost::split(v, str, ::isspace);
    if (v.size() == 2)
      env.push_back(atof(v[1].c_str()));
  } while (getline(file, str));
  cout << Envfile << " parsed successfully." << endl;
}

ostream &operator<<(ostream &stream, TSS_t const &s) {
  return stream << "{ " << s.TUindex << " " << s.TUorient << " " << s.TSS_pos
                << " " << s.TSS_strength << " }";
}

ostream &operator<<(ostream &stream, prot_t const &s) {
  return stream << "{ " << s.prot_name << " " << s.prot_pos << " }";
}

ostream &operator<<(ostream &stream, TTS_t const &s) {
  return stream << "{ " << s.TUindex << " " << s.TUorient << " " << s.TTS_pos
                << " " << s.TTS_proba_off << " }";
}

ostream &operator<<(ostream &stream, GFF_t const &s) {
  return stream << "{ " << s.seqname << " " << s.source << " " << s.feature
                << " " << s.start << " " << s.end << " " << s.score << " "
                << s.strand << " " << s.frame << " " << s.attribute << " }";
}

/* Guess genome size from the first annotation in GFF. This is dirty. */
uint get_genome_size(GFF_file &gff_df) {
  GFF_t full_genome = gff_df[0];
  return full_genome.end - full_genome.start + 1;
}

/* Get the transciption unit with the list of tts indexes belonging to TU. */
map<uint, vector<uint>> get_TU_tts(TSS_file &tss) {
  vector<uint> TU_values;
  std::transform(tss.begin(), tss.end(), std::back_inserter(TU_values),
                 [](TSS_t const &x) { return x.TUindex; });
  map<uint, vector<uint>> TU_tts;
  for (size_t i = 0, TU_index_val = TU_values[0]; i < TU_values.size();
       i++, TU_index_val = TU_values[i])
    TU_tts[TU_index_val].push_back(i);
  return TU_tts;
}

double f_prob_unhooked_rate(double sum_Kon, int DELTA_T,
                            size_t RNAPs_unhooked_nbr) {
  // np.exp(-sum_Kon*DELTA_T)/RNAPs_unhooked_nbr
  return (exp(-sum_Kon * double(DELTA_T)) / RNAPs_unhooked_nbr);
}

void random_choice(vector<int> &result, const vector<int> &array, uint n,
                   const vector<double> &proba) {
  vector<double> cum_probs(proba.size(), 0);
  vector<double> probs(proba);
  vector<int> available(array);
  result.clear();
  for (uint i = 0; i < n; ++i) {
    // Compute cumulated probs
    for (uint i = 0; i < probs.size() - 1; i++)
      cum_probs[i + 1] = accumulate(probs.begin(), probs.begin() + i + 2, 0.0);

    // Get an available index at random
    double k = rand() / double(RAND_MAX);
    uint index = find_if(cum_probs.begin(), cum_probs.end(),
                         [k](double p) -> bool { return p > k; }) -
                 cum_probs.begin();
    if (index >= available.size()) {
      cout << "probs: ";
      display_vector(proba);
      cout << "k = " << k << endl;
    }
    result.push_back(available[index]);
    double removed_p = probs[index];
    for (auto it = probs.begin(); it != probs.end(); it++)
      (*it) /= (1.0 - removed_p);

    // Remove the index from availables
    available.erase(available.begin() + index);
    probs.erase(probs.begin() + index);
    cum_probs.clear();
  }
}

void calc_sigma(vector<double> &Barr_sigma, Params *params) {
  for (auto sig = Barr_sigma.begin(); sig != Barr_sigma.end(); sig++) {
    *sig += params->DELTA_T *
            (params->TOPO_CONC * params->TOPO_CTE /
                 (1 + exp(params->k_TOPO * ((*sig) - params->x0_TOPO))) -
             params->GYRASE_CONC * params->GYRASE_CTE /
                 (1 + exp(params->k_GYRASE * ((*sig) - params->x0_GYRASE))));
  }
}