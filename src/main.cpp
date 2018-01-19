#include "Individual.h"
#include "Transcript.h"
#include "utils.h"
#include <boost/filesystem.hpp>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>
#include <thread>

static uint POP_SIZE = 10;
static uint GEN_MAX = 500;

vector<Individual *> &natural_selection(vector<Individual *> &population) {
  // sorts individuals by cost in increasing order
  sort(population.begin(), population.end(), [](Individual *a, Individual *b) {
    return a->get_fitness() < b->get_fitness();
  });

  // we keep the begining of the vector
  population.erase(population.begin() + population.size() / 2,
                   population.end());
  // cout << "\tkeeping the " << population.size() << " best individuals." <<
  // endl;

  return population;
}

int main(int argc, char **argv) {
  if (argc == 1)
    cout << "torsionEvol   path/to/params.ini   name_of_the_run   "
            "path/to/environment.dat   POP_SIZE   GEN_MAX  [ IndelPoissonMean "
            "InvP ]"
         << endl
         << "torsionEvol   path/to/params.ini   name_of_the_run   "
            "path/to/environment.dat   POP_SIZE   GEN_MAX  [ nThreads ]"
         << endl;
  if (argc == 1)
    return EXIT_FAILURE;
  POP_SIZE = atoi(argv[4]);
  GEN_MAX = atoi(argv[5]);
  if (argc > 7)
    Individual::set_mutation(atof(argv[4]), atof(argv[5]));

  int nthreads = thread::hardware_concurrency();
  int time_start = time(NULL);
  if (argc == 7)
    nthreads = atoi(argv[6]);
  cout << "number of threads: " << nthreads << endl;

  srand(time(NULL));

  // define the input/output directories
  Params *params = readIni(argv[1]);
  boost::filesystem::create_directories(argv[2]); // output folder
  ofstream scriptR;
  scriptR.open(string(argv[2]) + "/plot_results.R");
  scriptR << "setwd('" << boost::filesystem::initial_path().string() << '/'
          << argv[2] << "')" << endl;
  scriptR << "fit = as.matrix(read.table('fitnesses.txt', sep = ''))" << endl
          << "m = c()\nmini = c()" << endl;
  scriptR << "for (i in 1:length(fit[, 1])) {" << endl;
  scriptR << "    m[i] = mean(fit[i, ]) \nmini[i] = fit[i, 1]" << endl;
  scriptR << "}\nplot(m, col = 3, xlab = 'iterations', ylab = 'cost',\n";
  scriptR << "  main = 'Evolution and convergence', type = 'l', ylim = c(0, 1))"
          << endl;
  scriptR << "lines(mini, col = 2)" << endl
          << "legend('bottomleft', c('mean', 'mini'), lwd "
             "= 1, col = c(3, 2)) "
          << endl;
  scriptR.close();
  ofstream fitnesses;
  fitnesses.open(string(argv[2]) + "/fitnesses.txt");
  argv[1][strlen(argv[1]) - 10] = 0; // removes params.ini in path/to/params.ini
  const char *pth = argv[1];         // by setting "p" to 0 (end of string)

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

  // ================= Transcripts definition ================================

  // Map of transciption units with the list of tts belonging to TU. [ TU ]  =
  // (tts1, tts2, ... )
  map<uint, vector<uint>> TU_tts = get_TU_tts(tss);
  vector<Transcript> tr;
  for (TSS_t this_TSS : tss) {
    // Get the list of TTS that are in the same TU of this tss_id
    // TU_tts ex : { <0: {0, 1}>, <1: {2}> }  ==>  0 -> [0, 1]
    vector<uint> this_TU_tts = TU_tts[this_TSS.TUindex];

    float proba_rest = 1.0;
    uint k = 0;
    // Start scanning TTS.dat lines from the TUindex of this_TU_tts[0]
    while (proba_rest > 0 and k < this_TU_tts.size()) {
      TTS_t this_TTS = tts[this_TU_tts[k]];
      tr.push_back(Transcript(
          this_TSS.TUindex, this_TSS.TSS_pos, this_TTS.TTS_pos,
          (this_TSS.TUorient == '+') - (this_TSS.TUorient == '-'),
          this_TSS.TSS_strength * this_TTS.TTS_proba_off * proba_rest,
          params->DELTA_X));
      proba_rest *= (1.0 - this_TTS.TTS_proba_off);
      k++;
    }
  }

  // ====================== Topological barriers ============================

  vector<DNApos> Barr_fix; // Get the fixed topo barriers in a vector
  transform(prot.begin(), prot.end(), back_inserter(Barr_fix),
            [params](prot_t const &x) {
              return int(1 + x.prot_pos / params->DELTA_X);
            });

  // a bit of control on what happens
  cout << endl
       << "Starting with " << POP_SIZE
       << " individuals with the following genome:" << endl;
  for (Transcript t : tr) {
    cout << "TU n°" << t.TUindex_ << ": ";
    cout << t.TSS_ << '-' << t.TTS_ << " (length = " << t.size_ << "), strand ";
    cout << t.s_ << " with rate " << t.r_ << endl;
  }
  cout << endl;

  // Creation of the population:
  vector<Individual *> population;
  population.reserve(2 * POP_SIZE);
  uint genome_size = get_genome_size(gff_df_raw);
  for (uint i = 0; i < POP_SIZE; ++i)
    population.push_back(new Individual(genome_size, tr, Barr_fix));

  // =================== Here starts the genetic algorithm =====================
  for (uint generation_counter = 0; generation_counter < GEN_MAX;
       ++generation_counter) {
    printf("=========== GENERATION %d =============", 1 + generation_counter);

    // Mutate individuals
    // printf("\n mutation:\n");
    for (uint i = 0; i < POP_SIZE; i++) {
      Individual *new_guy = new Individual(*(population[i]));
      new_guy->mutate();
      population.push_back(new_guy);
    }

// Attribute a fitness to individuals (a cost, not really a fitness)
// printf("\n fitness:\n");
#pragma omp parallel for schedule(static) num_threads(nthreads)
    for (uint i = 0; i < population.size(); i++) {
      population[i]->update_fitness();
    }

    // Select the most adapted ones (fixed-size population)
    // printf("\n selection:\n");
    population = natural_selection(population);
    // printf("\n");

    // Display the costs
    // printf("\n printing:\n");
    for (auto indiv = population.begin(); indiv < population.end(); indiv++)
      fitnesses << " " << (*indiv)->get_fitness();
    fitnesses << endl;

    printf("\n");
  }

  int time_end = time(NULL);
  int time_full = time_end - time_start;
  cout << endl
       << "Simulation completed in " << time_full / 60 << "min "
       << time_full % 60 << "s. Deleting individuals..." << endl;

  // Cleaning
  for (Individual *indiv : population)
    delete indiv;
  delete params;
  fitnesses.close();

  // Plot with R
  string Rcmd("Rscript " + string(argv[2]) + "/plot_results.R");
  int result = system(Rcmd.c_str());

  // Exiting
  return (EXIT_SUCCESS + result);
}
