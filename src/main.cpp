#include <cstdlib>
#include "utils.h"
#include <boost/filesystem.hpp>

using namespace std;

int main(int argc, char** argv) {
	if (!argc) cout << "torsionEvol path/to/params.ini";
	Params* params = readIni(argv[1]);

	cout << "\tusing " << params->GFF << endl;
	cout << "\tusing " << params->TSS << endl;
	cout << "\tusing " << params->TTS << endl;

	// define the output directory
	boost::filesystem::create_directories(argv[2]);
	
	// define the input directory
	argv[1][strlen(argv[1])-10] = 0; 	//removes "params.ini" in string "path/to/params.ini"
	const char* pth = argv[1];	// by setting "p" to 0 (end of string)
	cout << "\tand input files in " << pth << endl;

	// read files
	vector<prot_t>* prot = readProt(pth + params->BARR_FIX);
	vector<TTS_t>* tts = readTTS(pth + params->TTS);
	vector<TSS_t>* tss = readTSS(pth + params->TSS);

	delete params;
	delete prot;
	delete tss;
	delete tts;
	return(EXIT_SUCCESS);
}
