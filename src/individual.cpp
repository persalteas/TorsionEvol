
//##########################################################
//         Transcription Process (Simulation)              #
//##########################################################

void Individual::start_transcribing(INI_file, output_dir)
{
    //###################### Params info ###################
    config = read_config_file(INI_file);

    // get inputs infos from the config file
    GFF_file = config.get('INPUTS', 'GFF');
    TSS_file = config.get('INPUTS', 'TSS');
    TTS_file = config.get('INPUTS', 'TTS');
    Prot_file = config.get('INPUTS', 'BARR_FIX');    

    // get values from the config file
    m = config.getfloat('GLOBAL', 'm');
    sigma_t = config.getfloat('GLOBAL', 'sigma_t');
    epsilon = config.getfloat('GLOBAL', 'epsilon');

    SIGMA_0 = config.getfloat('SIMULATION', 'SIGMA_0');
    DELTA_X = config.getfloat('SIMULATION', 'DELTA_X');
    DELTA_T = config.getfloat('SIMULATION', 'DELTA_T');
    RNAPS_NB = config.getint('SIMULATION', 'RNAPS_NB');
    ITERATIONS_NB = config.getfloat('SIMULATION', 'ITERATIONS_NB');
    OUTPUT_STEP = config.getfloat('SIMULATION', 'OUTPUT_STEP');

    GYRASE_CONC = config.getfloat('SIMULATION', 'GYRASE_CONC');
    TOPO_CONC = config.getfloat('SIMULATION', 'TOPO_CONC');
    TOPO_CTE = config.getfloat('SIMULATION', 'TOPO_CTE');
    GYRASE_CTE = config.getfloat('SIMULATION', 'GYRASE_CTE');
    TOPO_EFFICIENCY = config.getfloat('SIMULATION', 'TOPO_EFFICIENCY');
    k_GYRASE = config.getfloat('SIMULATION', 'k_GYRASE');
    x0_GYRASE = config.getfloat('SIMULATION', 'x0_GYRASE');
    k_TOPO = config.getfloat('SIMULATION', 'k_TOPO');
    x0_TOPO = config.getfloat('SIMULATION', 'x0_TOPO');
    //SIGMA_0 = 0 #((-np.log(((GYRASE_CONC*GYRASE_CTE)/TOPO_CONC*TOPO_CTE)-1))/k)+x_0
    //$print("SIGMA_0 --> ", SIGMA_0)

    // define the output directory
    os.makedirs(output_dir, exist_ok=True);
    
    // path to the input files (remove the "params.ini" from the path)
    pth = INI_file[:-10];
    gff_df_raw = load_gff(pth+GFF_file);
    tss = load_tab_file(pth+TSS_file);
    tts = load_tab_file(pth+TTS_file);
    prot = load_tab_file(pth+Prot_file);

    // TSS_pos
    TSS_pos = (tss['TSS_pos'].values/DELTA_X).astype(int);
    
    // Kon
    Kon = tss['TSS_strength'].values;
    
    // Poff
    Poff = tts['TTS_proba_off'].values;

    // Genome size
    genome_size = get_genome_size(gff_df_raw);
    gff_df = rename_gff_cols(gff_df_raw);

    // Dict of transciption units with the list of tts belonging to TU.
    TU_tts = get_TU_tts(tss, tts);

    // The RNAPs id
    RNAPs_id = np.full(RNAPS_NB, range(0, RNAPS_NB), dtype=int);

    // The position of RNAPs
    RNAPs_pos = np.full(RNAPS_NB, NaN) #np.zeros(RNAPS_NB, dtype=int);

    // RNAPs_last_pos
    RNAPs_last_pos = np.full(RNAPS_NB, NaN) #np.zeros(RNAPS_NB, dtype=int);

    // Strands orientation
    strands = str2num(gff_df['strand'].values);

    // list of all possible transcripts
    tr_id, tr_strand, tr_start, tr_end, tr_rate, tr_size, ts_beg_all_trs, ts_remain_all = get_tr_info(tss, tts, TU_tts, Kon, Poff);

    // convert all variables to numpy array
    tr_id = np.array(tr_id);
    tr_strand = np.array(tr_strand);

    tr_start = np.array(tr_start)/DELTA_X;
    tr_start = tr_start.astype(int64);

    tr_end = np.array(tr_end)/DELTA_X;
    tr_end = tr_end.astype(int64);

    tr_rate = np.array(tr_rate);

    tr_size = np.array(tr_size)/DELTA_X;
    tr_size = tr_size.astype(int64);

    ts_beg_all_trs = np.array(ts_beg_all_trs);

    ts_remain_all = np.array(ts_remain_all)/DELTA_X;
    ts_remain_all = ts_remain_all.astype(int64);
    
    // The number of times transcripts has been transcribed
    tr_nbr = np.zeros(len(tr_id), dtype=int);

    genome = int(genome_size/DELTA_X);

    Barr_fix = (prot['prot_pos'].values/DELTA_X).astype(int);
    
    // just for the echo we can assign it directely
    Barr_pos = np.copy(Barr_fix);
    Dom_size = np.ediff1d(Barr_pos);
    Dom_size = np.append(Dom_size, genome-Barr_fix[-1]+Barr_fix[0]); // !! change Barr_fix to Barr_pos case : O | |
    
    Barr_type = np.full(len(Barr_fix), 0, dtype=int);
    Barr_sigma = np.full(len(Barr_fix), SIGMA_0);
    
    // here we need to make an Barr_ts_remain
    // to track the position of each RNAPol
    // each position in Barr_ts_remain is associated with the same position in Barr_pos
    Barr_ts_remain = np.full(len(Barr_fix), NaN); // The Barr_ts_remain of fixed barr is NaN


    //######## Variables used to get the coverage ##########

    id_shift_fwd = list(range(1, genome));
    id_shift_fwd.append(0);
    id_shift_fwd = np.array(id_shift_fwd);
    id_shift_bwd = list(range(0, genome-1));
    id_shift_bwd.insert(0, genome-1);
    id_shift_bwd = np.array(id_shift_bwd);

    cov_bp = np.arange(0, genome_size, DELTA_X);
    cov_bp = np.resize(cov_bp, genome);

    //mmmm sigma = np.full(genome, SIGMA_0)

    //##########################################################
    //                 initiation of values                    #
    //##########################################################

    // save the time when RNApoly is starting trasncribing a specific transcript
    tr_times = col.defaultdict(list);

    // numpy array where will save all RNAPs info
    save_RNAPs_info = np.full([RNAPS_NB, 2, int(ITERATIONS_NB/DELTA_T)], np.nan); // nbr d'ele (cols)

    // the same for transcripts info
    save_tr_info = np.full([len(tr_id), 2, int(ITERATIONS_NB/DELTA_T)], np.nan);

    // in those variables, we will save/append info in each time step to save them as --> all_res ;-)
    save_Barr_sigma = list();
    save_Dom_size = list();
    save_mean_sig_wholeGenome = list();

    //########## Go !

    RNAPs_unhooked_id = np.copy(RNAPs_id);

    RNAPs_strand = np.full(RNAPS_NB, NaN);
    ts_beg = np.full(RNAPS_NB, NaN);
    ts_remain = np.full(RNAPS_NB, NaN);
    // RNAPs_tr will contain the id of the picked transcript
    RNAPs_tr = np.full(RNAPS_NB, -1, dtype=(int64));
    // get the TSSs ids
    tss_id = tss.index.values;

    // in the case of RNAP_NBR = 0
    RNAPs_hooked_id = [];

    for (int t{0}; t < static_cast<int>(ITERATIONS_NB/DELTA_T); ++t)
    {
        // we need to know each TSS belong to which Domaine
        TSS_pos_idx = np.searchsorted(Barr_pos, TSS_pos);

        // after knowing the domaine of each TSS we can get sigma
        sigma_tr_start = Barr_sigma[TSS_pos_idx-1];

        // get the initiation rates
        init_rate = f_init_rate(tr_rate, sigma_tr_start, sigma_t, epsilon, m);

        sum_init_rate = np.sum(init_rate);

        prob_init_rate = f_prob_init_rate(init_rate, sum_init_rate, DELTA_T);

        if (np.size(RNAPs_unhooked_id)!=0)
        {
            // get the unhooked rates
            prob_unhooked_rate = f_prob_unhooked_rate(sum_init_rate, DELTA_T, len(RNAPs_unhooked_id));
            // craete the numpy array
            prob_unhooked_rate = np.full(len(RNAPs_unhooked_id), prob_unhooked_rate);
            all_prob = np.concatenate([prob_init_rate, prob_unhooked_rate]);

            // create the numpy array that will contains [ nTSS , Unhooked RNAPS ]
            tss_and_unhooked_RNAPs = np.concatenate([tss_id, np.full(len(RNAPs_unhooked_id), -1, dtype=int)]);

            // pick up
            picked_tr = np.random.choice(tss_and_unhooked_RNAPs, len(RNAPs_unhooked_id), replace=False, p=all_prob) #RNAPs_unhooked_id;

            // This is the KEY !
            picked_tr_hooked_id = picked_tr[np.where(picked_tr!=-1)[0]];
            picked_tr_unhooked_id = picked_tr[np.where(picked_tr==-1)[0]];

            new_RNAPs_hooked_id = RNAPs_unhooked_id[np.where(picked_tr==picked_tr_hooked_id)];
    
            RNAPs_tr[new_RNAPs_hooked_id] = picked_tr[picked_tr!=-1];
            RNAPs_strand[new_RNAPs_hooked_id] = tr_strand[picked_tr[np.where(picked_tr!=-1)]];

            // The new position of each polymerase
            // if there is no RNAP already at this position
            RNAPs_pos[new_RNAPs_hooked_id] = tr_start[picked_tr[np.where(picked_tr!=-1)]].astype(int);

            // take the position and use them to get the index in which u will insert them in Barr_pos array
            Barr_pos_RNAPs_idx = np.searchsorted(Barr_pos, RNAPs_pos[new_RNAPs_hooked_id]);

            // after getting the idx, we start inserting
            Barr_pos = np.insert(Barr_pos, Barr_pos_RNAPs_idx, RNAPs_pos[new_RNAPs_hooked_id]);
            Dom_size = np.ediff1d(Barr_pos);
            Dom_size = np.append(Dom_size, genome-Barr_pos[-1]+Barr_pos[0]);
            Barr_type = np.insert(Barr_type, Barr_pos_RNAPs_idx, RNAPs_strand[new_RNAPs_hooked_id]);

            // Now Sigma
            Barr_sigma = np.insert(Barr_sigma, Barr_pos_RNAPs_idx, Barr_sigma[Barr_pos_RNAPs_idx-1]);
            
            // RNAPs_last_pos
            RNAPs_last_pos[new_RNAPs_hooked_id] = tr_end[picked_tr];
            ts_beg[new_RNAPs_hooked_id] = 0;
            ts_remain[new_RNAPs_hooked_id] = ts_remain_all[picked_tr];
            Barr_ts_remain = np.insert(Barr_ts_remain, Barr_pos_RNAPs_idx, ts_remain[new_RNAPs_hooked_id]);
            RNAPs_hooked_id = np.where(RNAPs_tr!=-1)[0];
        }

        ts_beg[RNAPs_hooked_id] += 1;
        ts_remain[RNAPs_hooked_id] -= 1;

        // save the time when RNApoly FINISHS trasncribing a specific transcript
        for (x : RNAPs_tr[np.where(ts_remain==0)])
            tr_times[x].append(t*DELTA_T); // + 0.5

        tr_nbr[RNAPs_tr[np.where(ts_remain==0)]]+=1;

        // look in the net : numpy where two conditions
        Barr_ts_remain[np.where(Barr_type == -1)]-=1;
        Barr_ts_remain[np.where(Barr_type == 1)]-=1;

        // Get the index of RNAPs to remove
        rm_RNAPs_idx = np.where(Barr_ts_remain == 0)[0];

        // recover sigma value of the removed position
        removed_sigma = Barr_sigma[rm_RNAPs_idx];
        removed_dom_size = Dom_size[rm_RNAPs_idx];

        // recover the old_dom_size : the size of the previous domaine before combination/merging
        old_dom_size = Dom_size[rm_RNAPs_idx-1];
        old_sigma = Barr_sigma[rm_RNAPs_idx-1];


        // update Dom_size 
        //Dom_size[rm_RNAPs_idx-1] += removed_dom_size 
        // or
        Dom_size = np.ediff1d(Barr_pos);
        Dom_size = np.append(Dom_size, genome-Barr_fix[-1]+Barr_fix[0]);
        
        Barr_sigma[rm_RNAPs_idx-1] = (old_dom_size*old_sigma+removed_dom_size*removed_sigma)/(old_dom_size+removed_dom_size);
        
        // and reomve them
        Barr_pos = np.delete(Barr_pos, rm_RNAPs_idx);
        Barr_type = np.delete(Barr_type, rm_RNAPs_idx);
        Barr_ts_remain = np.delete(Barr_ts_remain, rm_RNAPs_idx);
        Barr_sigma = np.delete(Barr_sigma, rm_RNAPs_idx);
        Dom_size = np.delete(Dom_size, rm_RNAPs_idx);

        // update the RNAPs_tr array
        RNAPs_tr[np.where(ts_remain==0)] = -1;
        // update the RNAPs_unhooked_id based on RNAPs_tr
        RNAPs_unhooked_id = np.where(RNAPs_tr==-1)[0];

        // reset the arrays
        RNAPs_strand[RNAPs_unhooked_id] = NaN;
        RNAPs_pos[RNAPs_unhooked_id] = NaN;
        RNAPs_last_pos[RNAPs_unhooked_id] = NaN;
        ts_beg[RNAPs_unhooked_id] = NaN;
        ts_remain[RNAPs_unhooked_id] = NaN;

        Barr_pos[np.where(Barr_type == -1)]-=1;
        Barr_pos[np.where(Barr_type == 1)]+=1;

        // Update the position of polymerases still transcribing
        RNAPs_pos[np.where(RNAPs_strand == 1)]+=1;
        RNAPs_pos[np.where(RNAPs_strand == -1)]-=1;

        //plt.scatter(t*DELTA_T, tr_nbr[0])

        // Update the Dom_size (+1 or -1)
        Dom_size = np.ediff1d(Barr_pos);
        Dom_size = np.append(Dom_size, genome-Barr_pos[-1]+Barr_pos[0]);
        
        // UPDATE SIGMA
        // R_plus_pos : the ids of RNA pol in the + strand
        R_plus_pos = np.where(Barr_type == 1)[0].astype(int);
        // R_minus_pos : the ids of RNA pol in the - strand
        R_minus_pos = np.where(Barr_type == -1)[0].astype(int);

        //### Extract all types of domaines (Those are ids of domaines)
        // Barr_type_ahead to make the extraction circular ;)
        Barr_type_ahead = np.roll(Barr_type, -1);
        // O +
        Barr_Dom_RPlus = np.where((Barr_type==0) & (Barr_type_ahead==1));
        // O -
        Barr_Dom_RMinus = np.where((Barr_type==0) & (Barr_type_ahead==-1));
        // O O
        Barr_Dom_Barr = np.where((Barr_type==0) & (Barr_type_ahead==0));
        // + +
        RPlus_Dom_RPlus = np.where((Barr_type==1) & (Barr_type_ahead==1));
        // - -
        RMinus_Dom_RMinus = np.where((Barr_type==-1) & (Barr_type_ahead==-1));
        // + -
        RPlus_Dom_RMinus = np.where((Barr_type==1) & (Barr_type_ahead==-1));
        // - +
        RMinus_Dom_RPlus = np.where((Barr_type==-1) & (Barr_type_ahead==+1));
        // - O
        RMinus_Dom_Barr = np.where((Barr_type==-1) & (Barr_type_ahead==0));
        // + O
        RPlus_Dom_Barr = np.where((Barr_type_ahead==0) & (Barr_type==+1));

        //### And then correct the value of Sigma in each case (before/after)
        corr_sig_Barr_Dom_RPlus = (Dom_size[Barr_Dom_RPlus]-1)/(Dom_size[Barr_Dom_RPlus]); // Sigma decrease x1
        corr_sig_Barr_Dom_RMinus = (Dom_size[Barr_Dom_RMinus]+1)/(Dom_size[Barr_Dom_RMinus]); // Sigma increase x1
        corr_sig_Barr_Dom_Barr = (Dom_size[Barr_Dom_Barr])/(Dom_size[Barr_Dom_Barr]); // Sigma FIX
        corr_sig_RPlus_Dom_RPlus = (Dom_size[RPlus_Dom_RPlus])/(Dom_size[RPlus_Dom_RPlus]); // Sigma FIX
        corr_sig_RMinus_Dom_RMinus = (Dom_size[RMinus_Dom_RMinus])/(Dom_size[RMinus_Dom_RMinus]); // Sigma FIX
        corr_sig_RPlus_Dom_RMinus = (Dom_size[RPlus_Dom_RMinus]+2)/(Dom_size[RPlus_Dom_RMinus]); // Sigma increase x2
        corr_sig_RMinus_Dom_RPlus = (Dom_size[RMinus_Dom_RPlus]-2)/(Dom_size[RMinus_Dom_RPlus]); // Sigma decrease x2
        corr_sig_RMinus_Dom_Barr = (Dom_size[RMinus_Dom_Barr]-1)/(Dom_size[RMinus_Dom_Barr]); // Sigma decrease x1
        corr_sig_RPlus_Dom_Barr = (Dom_size[RPlus_Dom_Barr]+1)/(Dom_size[RPlus_Dom_Barr]); // Sigma increase x1

        //## Multiply Sigma *= Corr (Each sigma value correspond to an specific domaine)
        Barr_sigma[Barr_Dom_RPlus] *= corr_sig_Barr_Dom_RPlus;
        Barr_sigma[Barr_Dom_RMinus] *= corr_sig_Barr_Dom_RMinus;
        Barr_sigma[Barr_Dom_Barr] *= corr_sig_Barr_Dom_Barr;
        Barr_sigma[RPlus_Dom_RPlus] *= corr_sig_RPlus_Dom_RPlus;
        Barr_sigma[RMinus_Dom_RMinus] *= corr_sig_RMinus_Dom_RMinus;
        Barr_sigma[RPlus_Dom_RMinus] *= corr_sig_RPlus_Dom_RMinus;
        Barr_sigma[RMinus_Dom_RPlus] *= corr_sig_RMinus_Dom_RPlus;
        Barr_sigma[RMinus_Dom_Barr] *= corr_sig_RMinus_Dom_Barr;
        Barr_sigma[RPlus_Dom_Barr] *= corr_sig_RPlus_Dom_Barr;

        //## Now calculate the SC generated in each domaine
        // RNAPs_genSC_all : contains an array of RNAPs_genSC that should be added or substracted from each domaine
        RNAPs_genSC_all = RNAPs_genSC/Dom_size;

        // Now update the value of sigma
        Barr_sigma[Barr_Dom_RPlus] -= RNAPs_genSC_all[Barr_Dom_RPlus];
        Barr_sigma[Barr_Dom_RMinus] += RNAPs_genSC_all[Barr_Dom_RMinus];
        Barr_sigma[RPlus_Dom_RMinus] += 2*RNAPs_genSC_all[RPlus_Dom_RMinus];
        Barr_sigma[RMinus_Dom_RPlus] -= 2*RNAPs_genSC_all[RMinus_Dom_RPlus];
        Barr_sigma[RMinus_Dom_Barr] -= RNAPs_genSC_all[RMinus_Dom_Barr];
        Barr_sigma[RPlus_Dom_Barr] += RNAPs_genSC_all[RPlus_Dom_Barr];
        

        // Now calc_sigma
        Barr_sigma = calc_sigma(Barr_sigma, GYRASE_CONC, k_GYRASE, x0_GYRASE, GYRASE_CTE, TOPO_CONC, k_TOPO, x0_TOPO, TOPO_CTE, DELTA_T);

        mean_sig_wholeGenome = np.sum(Barr_sigma*Dom_size)/genome;
        
        // Update the initiation rate        
        init_rate = f_init_rate(tr_rate, sigma_tr_start, sigma_t, epsilon, m);
        
        if (t%OUTPUT_STEP == 0)
        {
            // save all informations to npz file    
            // RNAPs_info
            save_RNAPs_info[:, 0, t] = RNAPs_tr; 
            save_RNAPs_info[:, 1, t] = RNAPs_pos;
            // tr_info
            save_tr_info[:, 0, t] = tr_nbr;
            save_tr_info[:, 1, t] = init_rate;
        }
        save_Barr_sigma.append(Barr_sigma);
        save_Dom_size.append(Dom_size);
        save_mean_sig_wholeGenome.append(mean_sig_wholeGenome);
    }
    save_Barr_sigma = np.array(save_Barr_sigma);
    save_Dom_size = np.array(save_Dom_size);
    save_mean_sig_wholeGenome = np.array(save_mean_sig_wholeGenome);
    save_files(output_dir, Barr_pos, Barr_type, Dom_size, Barr_ts_remain, Barr_sigma, tr_nbr, tr_times, save_RNAPs_info, save_tr_info, save_Barr_sigma, save_Dom_size, save_mean_sig_wholeGenome, DELTA_X, RNAPs_genSC, RNAPs_tr, RNAPs_pos, RNAPs_unhooked_id, init_rate, Kon, RNAPS_NB, SIGMA_0, GYRASE_CONC, TOPO_CONC);

    printf("Simulation completed successfully !! \nNumber of transcripts : \n");
    for i, v in enumerate(tr_nbr):
        print("Transcript{} : {}".format(i, v))

    return (GFF_file, TSS_file, TTS_file,
            ITERATIONS_NB, RNAPS_NB,
            tr_nbr, tr_times, init_rate, 
            RNAPs_tr, RNAPs_pos, RNAPs_unhooked_id,
            save_RNAPs_info, save_tr_info, save_Barr_sigma, save_Dom_size,
            cov_bp, tr_end);
}

