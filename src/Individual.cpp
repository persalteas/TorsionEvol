#include "Individual.h"
#include "RNAP.h"
#include <random>

#define N_DELTA_X 300

double Individual::_mean = 5,
       Individual::_p = 0.5; // modifiable with command line arguments 3 and 4
Params *Individual::_params = nullptr; // set later to params.ini
float Individual::_RNAPs_genSC = 0.1;  // not modifiable (yet)
vector<double> Individual::_target_envir;
Random Individual::_rand_generator = Random();

Individual::Individual(const Individual &indiv) {
  _genome_size = indiv._genome_size;
  _genes = indiv._genes;
  _barr_fix = indiv._barr_fix;
  _fitness = indiv._fitness;
}

Individual::Individual(unsigned genome_size, vector<Transcript> &tr,
                       vector<DNApos> &Barr_fix)
    : _genome_size(genome_size) {
  _genes = tr;
  _barr_fix = Barr_fix;
  _fitness = 1.0;
}

Individual &Individual::operator=(const Individual &indiv) {
  if (this != &indiv) {
    _genome_size = indiv._genome_size;
    _genes = indiv._genes;
    _barr_fix = indiv._barr_fix;
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
  std::vector<size_t> Dom_size(1 + _genes.size(), 0);

  // get the start of each gene on a segment [1, _genome_size+1]
  for (Transcript &tr : _genes) {
    if (tr.s_ == 1) { // (+) strand
      genes_pos.push_back(tr.start_);
      genes_pos.push_back(tr.end_);
    } else { // (-) strand
      genes_pos.push_back(tr.end_);
      genes_pos.push_back(tr.start_);
    }
  }
  genes_pos.push_back(_genome_size / _params->DELTA_X + 1);

  // compute the size of domains
  Dom_size[0] = genes_pos[0] - 1;
  for (uint i = 2; i < genes_pos.size(); i += 2) {
    Dom_size[i / 2] = genes_pos[i] - genes_pos[i - 1] - 1;
  }
  // std::cout << "\tbefore mutations: \n";
  // display_state();

  // Perform indels
  uint n_indels = _rand_generator.indels_nb();
  uint indel_nb = 0;
  int dom, start;
  int sense;
  DNApos pos;

  // cout << endl;
  // display_state();
  while (indel_nb < n_indels) {
    dom = get_rnd_dom_btwn_genes(Dom_size);
    // cout << "\tin domain " << dom << " (from ";
    // if (dom != -1)
    //   cout << genes_pos[2 * dom + 1] + 1 << " to "
    //        << genes_pos[2 * dom + 2] - 1;
    // else
    //   cout << "1 to " << genes_pos[0] - 1;
    // cout << "), ";
    pos = get_rnd_pos_in_domain(dom, genes_pos, Dom_size);

    if (rand() % 2) {
      // insertion
      std::cout << "\tinsertion in " << pos;
      sense = 1;
    } else {
      // deletion
      std::cout << "\tdeletion in  " << pos;
      sense = -1;
    }
    start = (2 + 2 * dom);
    if (int(Dom_size[1 + dom]) + sense >= 0) {
      // Update domain size
      cout << endl;
      Dom_size[1 + dom] += sense;

      for (uint i = start; i < genes_pos.size() - 1; i++) {
        // Update local vector of genes positions
        genes_pos[i] += sense;
        // Update Transcript objects
        if (!(i % 2)) {
          _genes[i / 2].shift(sense);
        }
      }
      // Update topo barriers
      for (uint i = 0; i < _barr_fix.size(); i++) {
        if (_barr_fix[i] >= pos)
          _barr_fix[i] += sense;
      }
      _genome_size += sense * _params->DELTA_X;
    }
    // cout << endl;
    // display_state();
    indel_nb++;
  }

  // Perform an inversion (or not)
  if (_rand_generator.inversion_occurs()) {
    int d1, d2, dom1, dom2;
    DNApos pos1, pos2, first, sec;
    dom1 = get_rnd_dom_btwn_genes(Dom_size);
    dom2 = get_rnd_dom_btwn_genes(Dom_size);
    pos1 = get_rnd_pos_in_domain(dom1, genes_pos, Dom_size);
    pos2 = get_rnd_pos_in_domain(dom2, genes_pos, Dom_size);
    (pos1 < pos2) ? (first = pos1, sec = pos2) : (first = pos2, sec = pos1);
    (pos1 < pos2) ? (d1 = dom1, d2 = dom2) : (d1 = dom2, d2 = dom1);
    std::cout << "\tinversion of [" << first << "-" << sec << "]\n";
    // Update Transcript objects
    if (d2 > d1) {
      for (auto it = _genes.begin() + d1 + 1; it < _genes.begin() + d2 + 1;
           it++) {
        it->revert(first, sec);
      }
      std::reverse(_genes.begin() + d1 + 1, _genes.begin() + d2 + 1);
    }

    // Update topo barriers
    for (uint i = 0; i < _barr_fix.size(); i++) {
      if ((_barr_fix[i] >= first) and (_barr_fix[i] <= sec))
        _barr_fix[i] = first + sec - _barr_fix[i] - 1;
    }
    std::sort(_barr_fix.begin(), _barr_fix.end());

    // // Update gene_pos local vector (useless)
    // auto it2 = genes_pos.begin() + 2 * d1 + 2;
    // for (auto it = _genes.begin() + d1 + 1; it < _genes.begin() + d2 + 1;
    //      it++, it2 += 2) {
    //   (*it2) = (it->s_ == 1) * (it->start_) + (it->s_ == -1) * (it->end_);
    //   (*(it2 + 1)) = (it->s_ == 1) * (it->end_) + (it->s_ == -1) *
    //   (it->start_);
    // }
    // // Update size of domains (useless)
    // Dom_size[0] = genes_pos[0] - 1;
    // for (uint i = 2; i < genes_pos.size(); i += 2) {
    //   Dom_size[i / 2] = genes_pos[i] - genes_pos[i - 1] - 1;
    // }

    // display_state();
  }
  // display_state();
}

void Individual::update_fitness(void) {
  // Simulate the expression process. Update the expr_count_ attribute of
  // transcripts in _genes.
  // std::cout << "\trunning Meyer's simulation, with _barr_fix: ";
  // display_vector(*_barr_fix);
  estimate_exression();

  // Sum the number of transcripts of each gene to compute a frequency
  double total_transcripts = 0;
  for (auto tr = _genes.begin(); tr != _genes.end(); tr++)
    total_transcripts += tr->expr_count_;

  // _fitness = SUM[ | x[i] - e[i] | ] with x[i] the proportion of gene i, e[i]
  // the target proportion
  _fitness = 0.0;
  for (size_t i = 0; i < _target_envir.size(); i++)
    _fitness +=
        abs((_genes[i].expr_count_ / total_transcripts) - _target_envir[i]);
  std::cout << "\tcost = " << _fitness << std::endl;
}

Individual::~Individual(void) {}

int Individual::get_rnd_dom_btwn_genes(const vector<size_t> &Dom_size) {
  vector<double> domain_probs, cum_probs;
  double non_gene_cum_size = 0;

  domain_probs.reserve(Dom_size.size());
  cum_probs.reserve(Dom_size.size());

  // probability of domains
  non_gene_cum_size = double(vector_sum(Dom_size));
  for (auto i = Dom_size.begin(); i < Dom_size.end(); i++) {
    domain_probs.push_back(double(*i) / non_gene_cum_size);
  }
  // Compute cumulated probs
  for (auto it = domain_probs.begin() + 1; it < domain_probs.end() + 1; it++)
    cum_probs.push_back(std::accumulate(domain_probs.begin(), it, 0.0));
  // Get an index at random
  double k = rand() / double(RAND_MAX);
  return find_if(cum_probs.begin(), cum_probs.end(),
                 [k](double p) -> bool { return p > k; }) -
         cum_probs.begin() - 1;
}

DNApos Individual::get_rnd_pos_in_domain(int dom, vector<DNApos> &gene_pos,
                                         vector<size_t> &Dom_size) {
  // chose random DNApos in given domain
  // (Domain 0 is between gene 0 and gene 1 in my model, domain -1 is before
  // gene 0)
  if (dom == -1) {
    return 1 + rand() % (Dom_size[0] + 1);
  }
  return gene_pos[2 * dom + 1] + 1 + (rand() % (Dom_size[1 + dom] + 1));
}

void Individual::estimate_exression() {
  // ====================== Topological barriers ============================
  std::vector<uint> Barr_pos(_barr_fix.size());
  copy(_barr_fix.begin(), _barr_fix.end(), Barr_pos.begin());
  std::vector<int> Dom_size(_barr_fix.size(), 0);
  std::vector<int> Barr_type(_barr_fix.size(), 0);
  std::vector<double> Barr_sigma(_barr_fix.size(), _params->SIGMA_0);
  // For later efficiency while inserting RNAPs in the barriers:
  Barr_pos.reserve(_barr_fix.size() + _params->RNAPS_NB);
  Barr_type.reserve(_barr_fix.size() + _params->RNAPS_NB);
  Barr_sigma.reserve(_barr_fix.size() + _params->RNAPS_NB);
  Dom_size.reserve(_barr_fix.size() + _params->RNAPS_NB);

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
  init_rates.reserve(_genes.size());
  prob_init_rates.reserve(_genes.size());
  all_prob.resize(_genes.size() + N_RNAPs_unhooked);
  picked_tr.reserve(N_RNAPs_unhooked);
  rm_RNAPs_idx.reserve(N_RNAPs_unhooked);
  int genome = int(_genome_size / _params->DELTA_X);
  for (size_t t = 0; t < (size_t)_params->ITERATIONS_NB; ++t) {
    // std::cout << "\t\t=========== t = " << t << "===========" << endl;
    // std::cout << "\t\tt = " << t << std::endl;
    // display_state(Barr_pos, Barr_type);
    // getchar();
    TSS_pos_idx.clear();
    init_rates.clear();
    prob_init_rates.clear();
    // get init_rate (raw)
    for (auto it = _genes.begin(); it != _genes.end(); it++)
      init_rates.push_back(it->f_init_rate(_params->sigma_t, _params->epsilon,
                                           _params->m, Barr_pos, Barr_sigma));
    // compute the probabilities of expression of each tr
    sum_init_rates = vector_sum(init_rates);
    for (auto it = _genes.begin(); it != _genes.end(); it++)
      prob_init_rates.push_back(
          it->f_prob_init_rate(sum_init_rates, _params->DELTA_T));
    // std::cout << "\t\t" << _params->RNAPS_NB - N_RNAPs_unhooked
    // << " RNAPs working. ";
    // if there are some unhooked RNAPs, hook them
    if (N_RNAPs_unhooked) {
      // get the probability for a RNAP to stay unhooked
      prob_unhooked_rate = f_prob_unhooked_rate(
          sum_init_rates, _params->DELTA_T, N_RNAPs_unhooked);
      prob_unhooked_rates.assign(N_RNAPs_unhooked, prob_unhooked_rate);
      // create an array that will contains [ 1 2 ... nTR , -1 -1 -1 ... ]
      tss_and_unhooked_RNAPs.assign(_genes.size() + N_RNAPs_unhooked, -1);
      iota(tss_and_unhooked_RNAPs.begin(),
           tss_and_unhooked_RNAPs.begin() + _genes.size(), 0);
      // concatenation of the probabilities arrays
      copy(prob_init_rates.begin(), prob_init_rates.end(), all_prob.begin());
      copy(prob_unhooked_rates.begin(), prob_unhooked_rates.end(),
           all_prob.begin() + _genes.size());
      // pick up transcripts randomly
      random_choice(picked_tr, tss_and_unhooked_RNAPs, N_RNAPs_unhooked,
                    all_prob);
      // std::cout << "\t\ttir des nouveaux tr: ";
      // if all RNAPs are hooked, assign transcripts to RNAPs
      for (int i : picked_tr) {
        // std::cout << i << ' ';
        if (i != -1) {
          // std::cout << "(brin" << _genes[i].s_ << ") ";
          // Add a hooked polymerase in the vector, working with the chosen
          // transcript.
          RNAPs_hooked.push_back(RNAP(_genes[i]));
          N_RNAPs_unhooked--;
          // Add the RNAP in the list of topological barriers:
          pos = _genes[i].start_;
          index =
              std::find_if(Barr_pos.begin(), Barr_pos.end(),
                           [pos](DNApos barr) -> bool { return barr > pos; }) -
              Barr_pos.begin();
          Barr_pos.insert(Barr_pos.begin() + index, pos);
          Barr_type.insert(Barr_type.begin() + index, _genes[i].s_);
          Barr_sigma.insert(Barr_sigma.begin() + index, Barr_sigma[index - 1]);
        }
      }
      // std::cout << std::endl;
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

    // std::cout << "\t\tprogress: ";
    // for (RNAP pol : RNAPs_hooked)
    //   std::cout << abs(int(pol.last_pos_) - int(pol.pos_)) << " ";
    // std::cout << std::endl;
    // std::cout << "\t\tposition of barriers: ";
    // for (DNApos pos : Barr_pos)
    //   std::cout << pos << " ";
    // std::cout << std::endl;

    // Check if RNAPs have finished a transcript
    rm_RNAPs_idx.clear();
    for (uint i = 0; i < RNAPs_hooked.size(); i++) {
      if (RNAPs_hooked[i].hasfinished()) {
        // std::cout << "\t\ta pol has finished !" << endl;
        _genes[RNAPs_hooked[i].tr_id_].expr_count_++;
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

    // if (rm_RNAPs_idx.size()) {
    //   std::cout << "\t\tremoved " << rm_RNAPs_idx.size() << " polymerases."
    //             << std::endl;
    //   std::cout << "\t\tnew barriers: ";
    //   for (DNApos pos : Barr_pos)
    //     std::cout << pos << " ";
    //   std::cout << std::endl;
    //   std::cout << "\t\thooked pols: " << RNAPs_hooked.size() << std::endl;
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
    // std::cout << "\t\ttorsion: ";
    // display_vector(Barr_sigma);
  }
  // std::cout << "\tSimulation completed successfully !!" << std::endl;
  // for (uint i = 0; i < _genes.size(); i++)
  //   std::cout << "\tTranscript " << i << " : " << _genes[i].expr_count_
  //             << std::endl;
}

void Individual::display_state(void) {
  vector<char> sprite(_genome_size / _params->DELTA_X + 2, '.');
  // vector<char> scale(_genome_size / _params->DELTA_X + 2, ' ');
  sprite[_genome_size / _params->DELTA_X + 1] = 0;
  sprite[0] = ' ';
  // scale[_genome_size / _params->DELTA_X + 1] = 0;
  // scale
  // for (uint i = 1; i < scale.size(); i++) {
  //   if (!(i % 10))
  //     scale[i] = '|';
  //   // if (!((i - 1) % 10))
  //   // scale[i] = i / 10;
  // }
  // genes
  for (uint i = 0; i < _genes.size(); i++) {
    Transcript &tr = _genes[i];
    if (tr.s_ == 1) {
      for (DNApos_n j = tr.start_; j <= tr.end_; j++)
        sprite[j] = '>';
    } else {
      for (DNApos_n j = tr.start_; j >= tr.end_; j--)
        sprite[j] = '<';
    }
  }
  // barriers
  for (uint i = 0; i < _barr_fix.size(); i++) {
    sprite[_barr_fix[i]] = 'X';
  }
  std::cout << "\t" << reinterpret_cast<char *>(sprite.data()) << std::endl;
  // std::cout << "\t" << reinterpret_cast<char *>(scale.data()) << std::endl
  // << std::endl;
}
