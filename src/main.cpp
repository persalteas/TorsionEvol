#include "Individual.h"
#include "Transcript.h"
#include "utils.h"
#include <boost/filesystem.hpp>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>

static uint POP_SIZE = 10;
static uint GEN_MAX = 10;

vector<Individual *> &natural_selection(vector<Individual *> &population) {
  // sorts individuals by cost in increasing order
  std::sort(population.begin(), population.end(),
            [](Individual *a, Individual *b) {
              return a->get_fitness() < b->get_fitness();
            });

  // we keep the begining of the vector
  for (auto it = population.begin() + population.size() / 2;
       it < population.end(); it++) {
    delete (*it);
  }
  population.erase(population.begin() + population.size() / 2,
                   population.end());
  std::cout << "\tkeeping the " << population.size() << " best individuals."
            << std::endl;

  return population;
}

int main(int argc, char **argv) {
  if (argc == 1)
    std::cout << "torsionEvol path/to/params.ini path/to/working/folder "
                 "path/to/environment.dat[ IndelPoissonMean InvP ]"
              << std::endl;
  if (argc == 1)
    return EXIT_FAILURE;
  if (argc > 5)
    Individual::set_mutation(atof(argv[4]), atof(argv[5]));

  srand(time(NULL));

  // define the input/output directories
  Params *params = readIni(argv[1]);
  boost::filesystem::create_directories(argv[2]); // output folder
  ofstream scriptR;
  scriptR.open(std::string(argv[2]) + "/plot_fit.R");
  scriptR << "setwd('" << boost::filesystem::initial_path().string() << '/'
          << argv[2] << "')" << std::endl;
  scriptR << "fit = as.matrix(read.table('fitnesses.txt', sep = ''))"
          << std::endl
          << "m = c()\nmini = c()" << std::endl;
  scriptR << "for (i in 1:length(fit[, 1])) {" << std::endl;
  scriptR << "    m[i] = mean(fit[i, ]) \nmini[i] = fit[i, 1]" << std::endl;
  scriptR << "}\nplot(m, col = 3, xlab = 'iterations', ylab = 'cost',\n";
  scriptR << "  main = 'Evolution and convergence', type = 'l', ylim = c(0, 1))"
          << std::endl;
  scriptR << "lines(mini, col = 2)" << std::endl
          << "legend('bottomleft', c('mean', 'mini'), lwd "
             "= 1, col = c(3, 2)) "
          << std::endl;
  scriptR.close();
  ofstream fitnesses;
  fitnesses.open(std::string(argv[2]) + "/fitnesses.txt");
  argv[1][strlen(argv[1]) - 10] =
      0;                     // removes "params.ini" in "path/to/params.ini"
  const char *pth = argv[1]; // by setting "p" to 0 (end of string)

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
    uint k =
        0; // Start scanning TTS.dat lines from the TUindex of this_TU_tts[0]
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

  std::vector<DNApos> Barr_fix; // Get the fixed topo barriers in a vector
  std::transform(prot.begin(), prot.end(), back_inserter(Barr_fix),
                 [params](prot_t const &x) {
                   return int(1 + x.prot_pos / params->DELTA_X);
                 });

  // a bit of control on what happens
  std::cout << std::endl
            << "Starting with individuals with the following genome:"
            << std::endl;
  for (Transcript t : tr) {
    std::cout << "TU n°" << t.TUindex_ << ": ";
    std::cout << t.TSS_ << '-' << t.TTS_ << " (length = " << t.size_
              << "), strand ";
    std::cout << t.s_ << " with rate " << t.r_ << std::endl;
  }
  std::cout << std::endl;

  // Creation of the population:
  std::cout << "Creating " << POP_SIZE << " individuals...";
  vector<Individual *> population;
  population.reserve(2 * POP_SIZE);
  uint genome_size = get_genome_size(gff_df_raw);
  for (uint i = 0; i < POP_SIZE; ++i)
    population.push_back(new Individual(genome_size, tr, Barr_fix));
  std::cout << " ok" << std::endl;

  // =================== Here starts the genetic algorithm =====================
  uint generation_counter = 0;
  while (generation_counter < GEN_MAX) {
    printf("=========== GENERATION %d: =============", 1 + generation_counter);

    // Mutate individuals
    printf("\n mutation:\n");
    for (auto indiv = population.begin(); indiv < population.begin() + POP_SIZE;
         indiv++) {
      Individual *new_guy = new Individual(**indiv);
      new_guy->mutate();
      population.push_back(new_guy);
    }

    // Attribute a fitness to individuals (a cost, not really a fitness)
    printf("\n fitness:\n");
    std::for_each(population.begin(), population.end(),
                  std::mem_fun(&Individual::update_fitness));

    // Select the most adapted ones (fixed-size population)
    printf("\n selection:\n");
    population = natural_selection(population);
    printf("\n");

    // Display the costs
    printf("\n printing:\n");
    for (auto indiv = population.begin(); indiv < population.end(); indiv++)
      fitnesses << " " << (*indiv)->get_fitness();
    fitnesses << std::endl;

    printf("\n");
    generation_counter++;
  }

  // Cleaning
  std::cout << std::endl
            << "Simulation completed. Deleting individuals..." << std::endl;
  for (Individual *indiv : population)
    delete indiv;
  delete params;
  fitnesses.close();

  // Exiting
  return (EXIT_SUCCESS);
}
