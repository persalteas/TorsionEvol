#include "Individual.h"
#include "RNAP.h"
#include <random>

double Individual::_mean = 5,
       Individual::_p = 0.5; // modifiable with command line arguments 3 and 4
Params *Individual::_params = nullptr; // set later to params.ini
float Individual::_RNAPs_genSC = 0.1;  // not modifiable (yet)
vector<double> Individual::_target_envir;
Random Individual::_rand_generator = Random();

Individual::Individual(const Individual &indiv) {
  _genome_size = indiv._genome_size;
  _genes = new std::vector<Transcript>(*(indiv._genes));
  _barr_fix = new std::vector<DNApos>(*(indiv._barr_fix));
  _fitness = indiv._fitness;
}

Individual::Individual(Individual &&indiv) {
  _genome_size = indiv._genome_size;
  _genes = indiv._genes; // undefined behaviour when indiv will be deleted ?
  _barr_fix =
      indiv._barr_fix; // undefined behaviour when indiv will be deleted ?
  _fitness = indiv._fitness;
}

Individual::Individual(unsigned genome_size, vector<Transcript> &tr,
                       vector<DNApos> &Barr_fix)
    : _genome_size(genome_size) {
  _genes = new std::vector<Transcript>(tr);
  _barr_fix = new std::vector<DNApos>(Barr_fix);
  _fitness = 1.0;
}

Individual &Individual::operator=(const Individual &indiv) {
  if (this != &indiv) {
    _genome_size = indiv._genome_size;
    delete _genes;
    delete _barr_fix;
    _genes = new std::vector<Transcript>(*(indiv._genes));
    _barr_fix = new std::vector<DNApos>(*(indiv._barr_fix));
    _fitness = indiv._fitness;
  }
  return *this;
}

Individual &Individual::operator=(Individual &&indiv) {
  if (this != &indiv) {
    delete _genes;
    delete _barr_fix;
    _genome_size = indiv._genome_size;
    _genes = indiv._genes; // undefined behaviour when indiv will be deleted ?
    _barr_fix =
        indiv._barr_fix; // undefined behaviour when indiv will be deleted ?
    _fitness = indiv._fitness;
  }
  return *this;
}

void Individual::set_mutation(double indel_nb, double inv_p) {
  _mean = indel_nb;
  _p = inv_p;
  _rand_generator = Random(_p, _mean);
}

void Individual::set_simulation_params(Params *params) { _params = params; }

void Individual::set_target_envir(vector<double> &env) { _target_envir = env; }

void Individual::mutate(void) {
  std::vector<DNApos> genes_pos;
  std::vector<size_t> Dom_size(2 * _genes->size(), 0);

  // get the start of each gene on a segment [0, _genome_size]
  for (Transcript &tr : (*_genes)) {
    if (tr.s_ == 1) { // (+) strand
      genes_pos.push_back(tr.TSS_);
      genes_pos.push_back(tr.TTS_ + 1);
    } else { // (-) strand
      genes_pos.push_back(tr.TTS_);
      genes_pos.push_back(tr.TSS_ + 1);
    }
  }
  // size of domains
  adjacent_difference(genes_pos.begin(), genes_pos.end(), Dom_size.begin());
  // (beware, adjacent_difference keeps a useless first element !)
  Dom_size.push_back(genes_pos[0] + (_genome_size - genes_pos.back()));

  std::cout << "\tbefore indels (genome size " << _genome_size << "): ";
  display_vector(genes_pos);

  // Perform indels
  uint n_indels = _rand_generator.indels_nb();
  uint indel_nb = 0;
  uint dom, start;
  DNApos pos;

  while (indel_nb < n_indels) {
    dom = get_rnd_dom_btwn_genes(Dom_size);
    pos = get_rnd_pos_in_domain(dom, genes_pos, Dom_size);

    if (rand() % 2) { // insertion
      std::cout << "\tinsertion in position " << pos << std::endl;
      Dom_size[2 + 2 * dom] += _params->DELTA_X;
      if ((dom == _genes->size() - 1) and (pos < genes_pos[0]))
        start = 0;
      else
        start = (2 + 2 * dom);
      for (uint i = start; i < genes_pos.size(); i++) {
        genes_pos[i] += _params->DELTA_X; // local vector of genes positions
        if (i % 2)
          (*_genes)[i / 2].shift(_params->DELTA_X, 1); // Transcript objects
      }

      _genome_size += _params->DELTA_X;
    } else { // deletion
      std::cout << "\tdeletion in position " << pos << std::endl;
      Dom_size[2 + 2 * dom] -= _params->DELTA_X;
      if ((dom == _genes->size() - 1) and (pos < genes_pos[0]))
        start = 0;
      else
        start = (2 + 2 * dom);
      for (uint i = start; i < genes_pos.size(); i++) {
        if (genes_pos[i] >= (uint)_params->DELTA_X)
          genes_pos[i] -= _params->DELTA_X;
        if (i % 2)
          (*_genes)[i / 2].shift(-_params->DELTA_X, -1);
      }
      _genome_size -= _params->DELTA_X;
    }

    indel_nb++;
  }
  std::cout << "\tafter indels (genome size " << _genome_size << "): ";
  display_vector(genes_pos);

  // Perform an inversion (or not)
  if (_rand_generator.inversion_occurs()) {
    DNApos pos1, pos2, first, sec;
    pos1 = get_rnd_pos_btwn_genes(genes_pos, Dom_size);
    pos2 = get_rnd_pos_btwn_genes(genes_pos, Dom_size);
    (pos1 < pos2) ? first = pos1, sec = pos2 : first = pos2, sec = pos1;
  }
}

void Individual::update_fitness(void) {
  // Simulate the expression process. Update the expr_count_ attribute of
  // transcripts in _genes.
  estimate_exression();

  // Sum the number of transcripts of each gene to compute a frequency
  double total_transcripts = 0;
  for (auto tr = (*_genes).begin(); tr != (*_genes).end(); tr++)
    total_transcripts += tr->expr_count_;

  // _fitness = SUM[ | x[i] - e[i] | ] with x[i] the proportion of gene i, e[i]
  // the target proportion
  _fitness = 0.0;
  for (size_t i = 0; i < _target_envir.size(); i++)
    _fitness +=
        abs(((*_genes)[i].expr_count_ / total_transcripts) - _target_envir[i]);
  std::cout << "\tcost = " << _fitness << ",";
}

Individual::~Individual(void) {
  delete _genes;
  delete _barr_fix;
}

uint Individual::get_rnd_dom_btwn_genes(const vector<size_t> &Dom_size) {
  vector<double> domain_probs, cum_probs;
  double non_gene_cum_size = 0;

  domain_probs.reserve((Dom_size.size() - 1) / 2.0);
  cum_probs.reserve((Dom_size.size() - 1) / 2.0);

  // probability of domains
  for (auto i = Dom_size.begin() + 2; i < Dom_size.end(); i += 2) {
    non_gene_cum_size += *i;
  }
  for (auto i = Dom_size.begin() + 2; i < Dom_size.end(); i += 2) {
    domain_probs.push_back(double(*i) / non_gene_cum_size);
  }
  // Compute cumulated probs
  for (auto it = domain_probs.begin(); it != domain_probs.end(); it++)
    cum_probs.push_back(std::accumulate(domain_probs.begin(), it, 0.0));
  // Get an index at random
  double k = rand() / double(RAND_MAX);
  uint dom = find_if(cum_probs.begin(), cum_probs.end(),
                     [k](double p) -> bool { return p > k; }) -
             cum_probs.begin() - 1;
  return dom;
}

DNApos Individual::get_rnd_pos_btwn_genes(vector<DNApos> &gene_pos,
                                          vector<size_t> &Dom_size) {
  // chose random domain in which to get a random DNApos
  uint dom = get_rnd_dom_btwn_genes(Dom_size);
  // chose random DNApos in domain
  // (Domain 0 is between gene 0 and gene 1 in my model)
  DNApos pos = gene_pos[2 * dom + 1] + (rand() % Dom_size[2 + 2 * dom]);
  if (pos > _genome_size) // pos is from the domain which contains base n°0
    return pos % _genome_size;
  return pos;
}

DNApos Individual::get_rnd_pos_in_domain(uint dom, vector<DNApos> &gene_pos,
                                         vector<size_t> &Dom_size) {
  // chose random DNApos in given domain
  // (Domain 0 is between gene 0 and gene 1 in my model)
  DNApos pos = gene_pos[2 * dom + 1] + (rand() % Dom_size[2 + 2 * dom]);
  if (pos > _genome_size) // pos is from the domain which contains base n°0
    return pos % _genome_size;
  return pos;
}

void Individual::estimate_exression() {
  // ====================== Topological barriers ============================
  std::vector<uint> Barr_pos(_barr_fix->size(), 0);
  copy(_barr_fix->begin(), _barr_fix->end(), Barr_pos.begin());
  std::vector<int> Dom_size(_barr_fix->size(), 0);
  std::vector<int> Barr_type(_barr_fix->size(), 0);
  std::vector<double> Barr_sigma(_barr_fix->size(), _params->SIGMA_0);
  // For later efficiency while inserting RNAPs in the barriers:
  Barr_pos.reserve(_barr_fix->size() + _params->RNAPS_NB);
  Barr_type.reserve(_barr_fix->size() + _params->RNAPS_NB);
  Barr_sigma.reserve(_barr_fix->size() + _params->RNAPS_NB);
  Dom_size.reserve(_barr_fix->size() + _params->RNAPS_NB);

  // ===================== Polymerases initialization ======================

  uint N_RNAPs_unhooked(_params->RNAPS_NB); // Available polymerases unhooked.
  std::vector<RNAP> RNAPs_hooked; // Vector of hooked polymerases, empty for now
  RNAPs_hooked.reserve(_params->RNAPS_NB);

  // ===================== Simulation ======================================
  std::vector<uint> TSS_pos_idx, rm_RNAPs_idx;
  std::vector<int> tss_and_unhooked_RNAPs, picked_tr;
  std::vector<double> init_rates, prob_init_rates, prob_unhooked_rates,
      all_prob;
  std::valarray<bool> isfinished;
  double sum_init_rates, prob_unhooked_rate;
  size_t index, j;
  DNApos pos;
  init_rates.reserve(_genes->size());
  prob_init_rates.reserve(_genes->size());
  all_prob.resize(_genes->size() + N_RNAPs_unhooked);
  picked_tr.reserve(N_RNAPs_unhooked);
  rm_RNAPs_idx.reserve(N_RNAPs_unhooked);
  int genome = int(_genome_size / _params->DELTA_X);
  for (size_t t = 0; t < (size_t)_params->ITERATIONS_NB; ++t) {
    // cout << "=========== t = " << t << "===========" << endl;
    TSS_pos_idx.clear();
    init_rates.clear();
    prob_init_rates.clear();
    // get init_rate (raw)
    for (auto it = _genes->begin(); it != _genes->end(); it++)
      init_rates.push_back(it->f_init_rate(_params->sigma_t, _params->epsilon,
                                           _params->m, Barr_pos, Barr_sigma));
    // compute the probabilities of expression of each tr
    sum_init_rates = vector_sum(init_rates);
    for (auto it = _genes->begin(); it != _genes->end(); it++)
      prob_init_rates.push_back(
          it->f_prob_init_rate(sum_init_rates, _params->DELTA_T));
    // cout << _params->RNAPS_NB - N_RNAPs_unhooked << " RNAPs working. ";
    // if there are some unhooked RNAPs, hook them
    if (N_RNAPs_unhooked) {
      // get the probability for a RNAP to stay unhooked
      prob_unhooked_rate = f_prob_unhooked_rate(
          sum_init_rates, _params->DELTA_T, N_RNAPs_unhooked);
      prob_unhooked_rates.assign(N_RNAPs_unhooked, prob_unhooked_rate);
      // create an array that will contains [ 1 2 ... nTR , -1 -1 -1 ... ]
      tss_and_unhooked_RNAPs.assign(_genes->size() + N_RNAPs_unhooked, -1);
      iota(tss_and_unhooked_RNAPs.begin(),
           tss_and_unhooked_RNAPs.begin() + _genes->size(), 0);
      // concatenation of the probabilities arrays
      copy(prob_init_rates.begin(), prob_init_rates.end(), all_prob.begin());
      copy(prob_unhooked_rates.begin(), prob_unhooked_rates.end(),
           all_prob.begin() + _genes->size());
      // pick up transcripts randomly
      random_choice(picked_tr, tss_and_unhooked_RNAPs, N_RNAPs_unhooked,
                    all_prob);
      // cout << "tir des nouveaux tr: ";
      // if all RNAPs are hooked, assign transcripts to RNAPs
      for (int i : picked_tr) {
        // cout << i << ' ';
        if (i != -1) {
          // cout << "(brin" << (*_genes)[i].s_ << ") ";
          // Add a hooked polymerase in the vector, working with the chosen
          // transcript.
          RNAPs_hooked.push_back(RNAP((*_genes)[i]));
          N_RNAPs_unhooked--;
          // Add the RNAP in the list of topological barriers:
          pos = (*_genes)[i].start_;
          index =
              std::find_if(Barr_pos.begin(), Barr_pos.end(),
                           [pos](DNApos barr) -> bool { return barr > pos; }) -
              Barr_pos.begin();
          Barr_pos.insert(Barr_pos.begin() + index, pos);
          Barr_type.insert(Barr_type.begin() + index, (*_genes)[i].s_);
          Barr_sigma.insert(Barr_sigma.begin() + index, Barr_sigma[index - 1]);
        }
      }
      // cout << endl;
      // Update the topological domains
      Dom_size.resize(Barr_pos.size());
      adjacent_difference(Barr_pos.begin(), Barr_pos.end(), Dom_size.begin());
      Dom_size.erase(Dom_size.begin()); // because unlike np.ediff1d,
                                        // adjacent_diff keeps the first
                                        // element.
      Dom_size.push_back(Barr_pos[0] + (genome - Barr_pos.back()));
    }
    // Move each polymerase on its transcript
    for (auto it = RNAPs_hooked.begin(); it != RNAPs_hooked.end(); it++)
      it->move();
    for (uint i = 0; i < Barr_pos.size(); i++)
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
    for (uint i = 0; i < RNAPs_hooked.size(); i++) {
      if (RNAPs_hooked[i].hasfinished()) {
        // cout << "a pol has finished !" << endl;
        (*_genes)[RNAPs_hooked[i].tr_id_].expr_count_++;
        rm_RNAPs_idx.push_back(i);
      }
    }
    // Unhook the polymerase(s)
    j = 0;
    for (auto it = rm_RNAPs_idx.begin(); it != rm_RNAPs_idx.end(); it++) {
      pos = RNAPs_hooked[(*it) - j].pos_;
      index = find_if(Barr_pos.begin(), Barr_pos.end(),
                      [pos](DNApos barr) -> bool { return barr > pos; }) -
              Barr_pos.begin();
      Barr_sigma[index - 2] = (Dom_size[index - 2] * Barr_sigma[index - 2] +
                               Dom_size[index - 1] * Barr_sigma[index - 1]) /
                              (Dom_size[index - 2] + Dom_size[index - 1]);
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
    std::adjacent_difference(Barr_pos.begin(), Barr_pos.end(),
                             Dom_size.begin());
    Dom_size.erase(Dom_size.begin());
    Dom_size.push_back(Barr_pos[0] + (genome - Barr_pos.back()));
    // Update Sigma:
    int a, b, d;
    for (uint i = 0; i < Barr_pos.size(); i++) {
      a = Barr_type[i];
      if (i + 1 == Barr_pos.size())
        b = Barr_type[0];
      else
        b = Barr_type[i + 1];
      d = Dom_size[i];
      switch (a) {
      case 0:
        switch (b) {
        case 1:
          Barr_sigma[i] = Barr_sigma[i] * (d - 1.0) / d - _RNAPs_genSC / d;
          break;
        case -1:
          Barr_sigma[i] = Barr_sigma[i] * (d + 1.0) / d + _RNAPs_genSC / d;
        }
        break;
      case -1:
        switch (b) {
        case 1:
          Barr_sigma[i] = Barr_sigma[i] * (d - 2.0) / d - 2 * _RNAPs_genSC / d;
          break;
        case 0:
          Barr_sigma[i] = Barr_sigma[i] * (d - 1.0) / d - _RNAPs_genSC / d;
        }
        break;
      case 1:
        switch (b) {
        case -1:
          Barr_sigma[i] = Barr_sigma[i] * (d + 2.0) / d + 2 * _RNAPs_genSC / d;
          break;
        case 0:
          Barr_sigma[i] = Barr_sigma[i] * (d + 1.0) / d + _RNAPs_genSC / d;
        }
      }
    }
    calc_sigma(Barr_sigma, _params);
    // cout << "torsion: ";
    // display_vector(Barr_sigma);
  }
  // std::cout << "Simulation completed successfully !! \nNumber of transcripts
  // : "<< endl;
  // for (uint i=0; i<_genes->size(); i++)
  //     std::cout << "Transcript " << i << " : " << (*_genes)[i].expr_count_ <<
  //     endl;
}
