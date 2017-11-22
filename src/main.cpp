#include <cstdlib>
#include "utils.h"

using namespace std;

int main(int argc, char** argv) {
	if (!argc) cout << "torsionEvol path/to/params.ini";
	Params* params = readIni(argv[1]);

	cout << "Welcome to TorsionEvol simulation.\n";
	cout << "using " << params->GFF << endl;

	return(EXIT_SUCCESS);
}
