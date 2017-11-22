#include <cstdlib>
#include "utils.h"

int main(int argc, char** argv)
{
    Params* params_ini = readIni(argv[1]);
    return EXIT_SUCCESS;
}