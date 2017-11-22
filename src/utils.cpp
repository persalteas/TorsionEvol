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
