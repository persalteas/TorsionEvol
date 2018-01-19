#include "Individual.h"
#include "RNAP.h"
#include <random>

double Individual::_mean = 5,
       Individual::_p = 0.5; // modifiable with command line arguments 3 and 4
Params *Individual::_params = nullptr; // set later to params.ini
float Individual::_RNAPs_genSC = 0.0;  // not modifiable (yet)
vector<double> Individual::_target_envir;
Random Individual::_rand_generator = Random();

Individual::Individual(unsigned genome_size, vector<Transcript> &tr,
                       vector<DNApos> &Barr_fix) {
  // First constructor used when loading data from param files
  _genome_size = genome_size;
  _fitness = 1.0;
  _genes = vector<Transcript>(tr);
  _barr_fix = vector<DNApos_n>(Barr_fix);
  size_t a = _barr_fix.size();
  size_t b = _params->RNAPS_NB;
  size_t c = a + b;
  size_t d = _genes.size();
  Barr_pos = vector<uint>(a, 0);
  Dom_size = vector<int>(a, 0);
  Dom_size_1 = vector<size_t>(1 + d, 0);
  Barr_type = vector<int>(a, 0);
  init_rates = vector<double>(d, 0);
  prob_init_rates = vector<double>(d, 0);
  Barr_sigma = vector<double>(a, _params->SIGMA_0);
  RNAPs_hooked = vector<RNAP>();
  rm_RNAPs_idx = vector<uint>();
  tss_and_unhooked_RNAPs = vector<int>();
  picked_tr = vector<int>();
  prob_unhooked_rates = vector<double>();
  all_prob = vector<double>();
  genes_pos = vector<DNApos_n>();
  Barr_pos.reserve(c);
  Barr_type.reserve(c);
  Barr_sigma.reserve(c);
  Dom_size.reserve(c);
  RNAPs_hooked.reserve(b);
  all_prob.reserve(d + b);
  picked_tr.reserve(b);
  rm_RNAPs_idx.reserve(b);
}

Individual::Individual(const Individual &indiv) {
  // Copy constructor used to duplicate an Indiv before mutation
  _genome_size = indiv._genome_size;
  _fitness = 1.0;
  _genes = vector<Transcript>(indiv._genes);
  _barr_fix = vector<DNApos_n>(indiv._barr_fix);
  size_t a = _barr_fix.size();
  size_t b = _params->RNAPS_NB;
  size_t c = a + b;
  size_t d = _genes.size();
  Barr_pos = vector<uint>(a, 0);
  Dom_size = vector<int>(a, 0);
  Dom_size_1 = vector<size_t>(1 + d, 0);
  Barr_type = vector<int>(a, 0);
  init_rates = vector<double>(d, 0);
  prob_init_rates = vector<double>(d, 0);
  Barr_sigma = vector<double>(a, _params->SIGMA_0);
  RNAPs_hooked = vector<RNAP>();
  rm_RNAPs_idx = vector<uint>();
  tss_and_unhooked_RNAPs = vector<int>();
  picked_tr = vector<int>();
  prob_unhooked_rates = vector<double>();
  all_prob = vector<double>();
  genes_pos = vector<DNApos_n>(2 * d, 0);
  Barr_pos.reserve(c);
  Barr_type.reserve(c);
  Barr_sigma.reserve(c);
  Dom_size.reserve(c);
  RNAPs_hooked.reserve(b);
  all_prob.reserve(d + b);
  picked_tr.reserve(b);
  rm_RNAPs_idx.reserve(b);
}

void Individual::set_mutation(double indel_nb, double inv_p) {
  _mean = indel_nb;
  _p = inv_p;
  _rand_generator = Random(_p, _mean);
}

void Individual::set_simulation_params(Params *params) { _params = params; }

void Individual::set_target_envir(vector<double> &env) { _target_envir = env; }

void Individual::mutate(void) {

  // get the start of each gene on a segment [1, _genome_size+1]
  for (uint i = 0; i < _genes.size(); i++) {
    Transcript &tr = _genes[i];
    if (tr.s_ == 1) { // (+) strand
      genes_pos[2 * i] = tr.start_;
      genes_pos[2 * i + 1] = tr.end_;
    } else { // (-) strand
      genes_pos[2 * i] = tr.end_;
      genes_pos[2 * i + 1] = tr.start_;
    }
  }
  genes_pos.push_back(_genome_size / _params->DELTA_X + 1);

  // compute the size of domains
  Dom_size_1[0] = genes_pos[0] - 1;
  for (uint i = 2; i < genes_pos.size(); i += 2) {
    Dom_size_1[i / 2] = genes_pos[i] - genes_pos[i - 1] - 1;
  }
  // cout << "\tbefore mutations: \n";
  // cout << "\tgenes positions: ";
  // display_vector(genes_pos);
  // cout << "\tbarriers positions: ";
  // display_vector(_barr_fix);
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
    dom = get_rnd_dom_btwn_genes();
    // cout << "\tin domain " << dom << " (from ";
    // if (dom != -1)
    //   cout << genes_pos[2 * dom + 1] + 1 << " to "
    //        << genes_pos[2 * dom + 2] - 1;
    // else
    //   cout << "1 to " << genes_pos[0] - 1;
    // cout << "), ";
    pos = get_rnd_pos_in_domain(dom);

    if (rand() % 2) {
      // insertion
      // cout << "\tinsertion in " << pos;
      sense = 1;
    } else {
      // deletion
      // cout << "\tdeletion in  " << pos;
      sense = -1;
    }
    start = (2 + 2 * dom);
    if (int(Dom_size_1[1 + dom]) + sense >= 0) {
      // Update domain size
      // cout << endl;
      Dom_size_1[1 + dom] += sense;
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
        if (_barr_fix[i] >= pos and sense == 1)
          _barr_fix[i]++;
        if (_barr_fix[i] > pos and sense == -1)
          _barr_fix[i]--;
      }
      _genome_size += sense * _params->DELTA_X;
    }
    // cout << "\tgenes positions: ";
    // display_vector(genes_pos);
    // cout << "\tbarriers positions: ";
    // display_vector(_barr_fix);
    // display_state();
    indel_nb++;
  }

  // Perform an inversion (or not)
  if (_rand_generator.inversion_occurs()) {
    int d1, d2, dom1, dom2;
    DNApos pos1, pos2, first, sec;
    dom1 = get_rnd_dom_btwn_genes();
    dom2 = get_rnd_dom_btwn_genes();
    pos1 = get_rnd_pos_in_domain(dom1);
    pos2 = get_rnd_pos_in_domain(dom2);
    (pos1 < pos2) ? (first = pos1, sec = pos2) : (first = pos2, sec = pos1);
    (pos1 < pos2) ? (d1 = dom1, d2 = dom2) : (d1 = dom2, d2 = dom1);
    // cout << "\tinversion of [" << first << "-" << sec << "]\n";
    // Update Transcript objects
    if (d2 > d1) {
      for (auto it = _genes.begin() + d1 + 1; it < _genes.begin() + d2 + 1;
           it++) {
        it->revert(first, sec);
      }
      reverse(_genes.begin() + d1 + 1, _genes.begin() + d2 + 1);
    }

    // Update topo barriers
    for (uint i = 0; i < _barr_fix.size(); i++) {
      if ((_barr_fix[i] >= first) and (_barr_fix[i] <= sec))
        _barr_fix[i] = first + sec - _barr_fix[i];
    }
    sort(_barr_fix.begin(), _barr_fix.end());
  }
  // cout << "\tgenes positions: ";
  // display_vector(genes_pos);
  // cout << "\tbarriers positions: ";
  // display_vector(_barr_fix);
  // display_state();
}

void Individual::update_fitness(void) {
  // Simulate the expression process. Update the expr_count_ attribute of
  // transcripts in _genes.
  // cout << "\trunning Meyer's simulation, with _barr_fix: ";
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
  // cout << "\tcost = " << _fitness << endl;
}

int Individual::get_rnd_dom_btwn_genes(void) {
  vector<double> domain_probs, cum_probs;
  double non_gene_cum_size = 0;

  domain_probs.reserve(Dom_size_1.size());
  cum_probs.reserve(Dom_size_1.size());

  // probability of domains
  non_gene_cum_size = double(vector_sum(Dom_size_1));
  for (auto i = Dom_size_1.begin(); i < Dom_size_1.end(); i++) {
    domain_probs.push_back(double(*i) / non_gene_cum_size);
  }
  // Compute cumulated probs
  for (auto it = domain_probs.begin() + 1; it < domain_probs.end() + 1; it++)
    cum_probs.push_back(accumulate(domain_probs.begin(), it, 0.0));
  // Get an index at random
  double k = rand() / double(RAND_MAX);
  return find_if(cum_probs.begin(), cum_probs.end(),
                 [k](double p) -> bool { return p > k; }) -
         cum_probs.begin() - 1;
}

DNApos Individual::get_rnd_pos_in_domain(int dom) {
  // chose random DNApos in given domain
  // (Domain 0 is between gene 0 and gene 1 in my model, domain -1 is before
  // gene 0)
  if (dom == -1) {
    return 1 + rand() % (Dom_size_1[0] + 1);
  }
  return genes_pos[2 * dom + 1] + 1 + (rand() % (Dom_size_1[1 + dom]));
}

void Individual::estimate_exression() {
  Barr_pos.resize(_barr_fix.size());
  copy(_barr_fix.begin(), _barr_fix.end(), Barr_pos.begin());
  Dom_size.clear();
  Dom_size.resize(_barr_fix.size(), 0);
  Barr_type.clear(); // IMPORTANT ! Deletes the old types of previous iteration.
  Barr_type.resize(_barr_fix.size(), 0);
  Barr_sigma.clear();
  Barr_sigma.resize(_barr_fix.size(), _params->SIGMA_0);
  tss_and_unhooked_RNAPs.clear();
  picked_tr.clear();
  prob_unhooked_rates.clear();
  RNAPs_hooked.clear();
  uint N_RNAPs_unhooked = _params->RNAPS_NB; // Available polymerases unhooked.
  all_prob.resize(_genes.size() + N_RNAPs_unhooked);

  // ===================== Simulation ======================================

  valarray<bool> isfinished;
  double sum_init_rates;
  size_t index, j;
  DNApos pos;
  int genome = int(_genome_size / _params->DELTA_X);
  // display_state();
  // cout << "\t\tposition of barriers: ";
  // display_vector(*);
  // cout << "\t\tsize of domains: ";
  // display_vector(*Dom_size);
  // cout << "\t\tactual torsions: ";
  // display_vector(*Barr_sigma);
  for (size_t t = 0; t < (size_t)_params->ITERATIONS_NB; ++t) {
    // cout << "\n\t\t=========== t = " << t << "===========" << endl;
    // cout << "\t\t" << _params->RNAPS_NB - N_RNAPs_unhooked
    //           << " RNAPs working. ";

    // ============== try to bind unbound RNAPs ===============
    if (N_RNAPs_unhooked) {
      // get the probability for a RNAP to chose each gene
      for (uint i = 0; i < _genes.size(); i++)
        init_rates[i] =
            _genes[i].f_init_rate(_params->sigma_t, _params->epsilon,
                                  _params->m, Barr_pos, Barr_sigma);
      sum_init_rates = vector_sum(init_rates);
      for (uint i = 0; i < _genes.size(); i++)
        prob_init_rates[i] =
            _genes[i].f_prob_init_rate(sum_init_rates, _params->DELTA_T);

      // get the probability for a RNAP to stay unhooked
      prob_unhooked_rates.assign(N_RNAPs_unhooked,
                                 f_prob_unhooked_rate(sum_init_rates,
                                                      _params->DELTA_T,
                                                      N_RNAPs_unhooked));
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
      // cout << "\t\ttir des nouveaux tr: ";
      for (int i : picked_tr) {
        // cout << i << ' ';
        if (i != -1) {
          // cout << "(brin" << _genes[i].s_ << ") ";
          // Add a hooked polymerase in the vector, working with the chosen
          // transcript.
          RNAPs_hooked.push_back(RNAP(_genes[i]));
          N_RNAPs_unhooked--;
          // Add the RNAP in the list of topological barriers:
          pos = _genes[i].start_;
          index = find_if(Barr_pos.begin(), Barr_pos.end(),
                          [pos](DNApos barr) -> bool { return barr > pos; }) -
                  Barr_pos.begin();
          Barr_pos.insert(Barr_pos.begin() + index, pos);
          Barr_type.insert(Barr_type.begin() + index, _genes[i].s_);
          Barr_sigma.insert(Barr_sigma.begin() + index,
                            Barr_sigma[index * (index != Barr_sigma.size())]);
          Dom_size.push_back(0); // just to increase its size
          adjacent_difference(Barr_pos.begin(), Barr_pos.end(),
                              Dom_size.begin());
          Dom_size[0] += (genome - Barr_pos.back());
        }
      }
      // cout << endl;
    }
    // cout << "\t\ttorsion avant mouvement: ";
    // display_vector(*Barr_sigma);

    // ================ Unbound RNAPs that have finished =================
    rm_RNAPs_idx.clear();
    for (uint i = 0; i < RNAPs_hooked.size(); i++) {
      if (RNAPs_hooked[i].hasfinished()) {
        _genes[RNAPs_hooked[i].tr_id_].expr_count_++;
        rm_RNAPs_idx.push_back(i);
      }
    }
    if (rm_RNAPs_idx.size()) {
      // cout << endl;
      // Unhook the polymerase(s)
      j = 0;
      for (auto it = rm_RNAPs_idx.begin(); it != rm_RNAPs_idx.end(); it++) {
        // retrieve the position of the RNAP in Barr_pos
        pos = RNAPs_hooked[(*it) - j].pos_;
        index = find_if(Barr_pos.begin(), Barr_pos.end(),
                        [pos](DNApos barr) -> bool { return barr > pos; }) -
                Barr_pos.begin();
        // cout << "\t\tunhooking " << *it << " from pos " << pos
        //      << ", index=" << index << endl;
        // remove it and update barriers, domains, and sigma
        if (index < Barr_sigma.size()) {
          Barr_sigma[index] = (Dom_size[index] * Barr_sigma[index] +
                               Dom_size[index - 1] * Barr_sigma[index - 1]) /
                              (Dom_size[index] + Dom_size[index - 1]);
          Dom_size[index] = Dom_size[index] + Dom_size[index - 1];
        } else { // we remove the last barr of the vector
          Barr_sigma[0] = (Dom_size[0] * Barr_sigma[0] +
                           Dom_size[index - 1] * Barr_sigma[index - 1]) /
                          (Dom_size[0] + Dom_size[index - 1]);
          Dom_size[0] = Dom_size[0] + Dom_size[index - 1];
        }
        Barr_pos.erase(Barr_pos.begin() + index - 1);
        Barr_type.erase(Barr_type.begin() + index - 1);
        Barr_sigma.erase(Barr_sigma.begin() + index - 1);
        Dom_size.erase(Dom_size.begin() + index - 1);

        adjacent_difference(Barr_pos.begin(), Barr_pos.end(), Dom_size.begin());
        Dom_size[0] += (genome - Barr_pos.back());
        RNAPs_hooked.erase(RNAPs_hooked.begin() + (*it) - j);
        N_RNAPs_unhooked++;
        j++;
      }
      // cout << "\t\tremoved " << rm_RNAPs_idx.size() << " polymerases."
      //           << endl;
      // cout << "\t\ttype of barriers: ";
      // display_vector(Barr_type);
      // cout << "\t\tposition of barriers: ";
      // display_vector(Barr_pos);
      // cout << "\t\tsize of domains: (" << Dom_size.size() << ") ";
      // display_vector(Dom_size);
      // cout << "\t\tactual torsions: (" << Barr_sigma.size() << ") ";
      // display_vector(Barr_sigma);
      // cout << endl;
    }

    // ========== Move each polymerase on its transcript =================
    for (auto it = RNAPs_hooked.begin(); it != RNAPs_hooked.end(); it++)
      it->move();
    for (uint i = 0; i < Barr_pos.size(); i++)
      Barr_pos[i] += Barr_type[i];

    // ================== Update domain vectors ==========================
    // Dom_size[i] and Barr_sigma[i] correspond to the domain before barrier
    // Barr_pos[i] of type Barr_type[i]
    adjacent_difference(Barr_pos.begin(), Barr_pos.end(), Dom_size.begin());
    Dom_size[0] += (genome - Barr_pos.back());
    double a, b, d;
    for (uint i = 0; i < Barr_sigma.size(); i++) {
      if (i)
        a = Barr_type[i - 1];
      else
        a = Barr_type.back();
      b = Barr_type[i];
      d = Dom_size[i];
      Barr_sigma[i] *= (d + (a - b) * (1.0 + _RNAPs_genSC)) / d;
    }
    calc_sigma(Barr_sigma, _params);

    // =================== display user infos ===========================
    // cout << "\t\tprogress: ";
    // for (RNAP pol : (*RNAPs_hooked))
    //   cout << abs(int(pol.last_pos_) - int(pol.pos_)) << " ";
    // cout << endl;
    // display_state();
    // cout << "\t\ttype of barriers    : ";
    // display_vector(Barr_type);
    // cout << "\t\tposition of barriers: ";
    // display_vector(Barr_pos);
    // cout << "\t\tsize of domains: (" << Dom_size.size() << ") ";
    // display_vector(Dom_size);
    // cout << "\t\tactual torsions: (" << Barr_sigma.size() << ") ";
    // display_vector(Barr_sigma);
    // cout << endl;
  }
  // cout << "\tSimulation completed successfully !!" << endl;
  // for (uint i = 0; i < _genes.size(); i++)
  //   cout << "\tTranscript " << i << " : " << _genes[i].expr_count_
  //             << endl;
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
  // fix barriers
  for (uint i = 0; i < _barr_fix.size(); i++) {
    sprite[_barr_fix[i]] = 'X';
  }
  // barriers
  for (uint i = 0; i < Barr_pos.size(); i++) {
    if (Barr_type[i]) {
      sprite[Barr_pos[i]] = 'O';
    } else {
      sprite[Barr_pos[i]] = 'X';
    }
  }
  cout << "\t" << reinterpret_cast<char *>(sprite.data()) << endl;
  // cout << "\t" << reinterpret_cast<char *>(scale.data()) << endl
  // << endl;
}
